
% ============== analyze_gpsense_ultimate.m (v5.8 - SNR对比图改为子图样式) ==============
%
% 功能:
%   本脚本是GPSense分析流程的终极整合版本。它集成了信号幅度(归一化与原始值)、
%   伪距和两种核心相位去趋势方法（滑动平均法 和 带通滤波法）的计算与可视化。
%   【V5.6更新】: 同步添加包含原始SNR的三合一绘图模块（时间轴版本）。
%   【V5.7更新】: 新增 SNR 滤波前后对比图，用于评估去噪效果。
%   【V5.8更新】: SNR 对比图改为上下子图样式，使用默认颜色。
%   【多系统更新】: 扩展以支持 Galileo (E) 和 QZSS (J) 卫星。
%
% 输出:
%   time_vector      - 有效数据点的时间向量 (datetime)
%   amplitude_norm   - 归一化后的信号幅度 (0-1范围)
%   pseudorange_detrend - 去趋势后的伪距 (米)
%   phase_detrend_movmean - 通过滑动平均法去趋势后的相位 (弧度)
%   phase_detrend_bpf     - 通过带通滤波法处理后的相位 (弧度)
%
function [time_vector, amplitude_norm, pseudorange_detrend, phase_detrend_movmean, phase_detrend_bpf] = analyze_gpsense_ultimate(obs_data, nav_data, target_satellite_id, epoch_range)

% --- 1. 参数检查与初始化 ---
if nargin < 3, error('请提供目标卫星ID (例如 ''G01'')。'); end
if nargin < 4 || isempty(epoch_range), epoch_range = 1:length(obs_data); end

% --- 物理常数 ---
c = 299792458.0; % 光速 (m/s)
k_boltzmann = 1.380649e-23; % 玻尔兹曼常数 (J/K)
T_noise_kelvin = 290.0; % 系统噪声温度 (K), 标准室温假设
Pn_noise_power = k_boltzmann * T_noise_kelvin; % 噪声功率

% ====================【核心修改：扩展卫星系统支持】====================
% --- 动态适配卫星系统 (G/C/E/J) ---
sys_char = upper(target_satellite_id(1));
if sys_char == 'G' || sys_char == 'E' || sys_char == 'J'
    % GPS (G), Galileo (E), QZSS (J) 均使用 L1/E1 频段
    pr_code = 'C1C'; dop_code = 'D1C'; snr_code = 'S1C';
    f_carrier = 1575.42e6; % GPS L1 / GAL E1 / QZS L1 频率
elseif sys_char == 'C'
    % BeiDou (C) 使用 B1I 频段
    pr_code = 'C2I'; dop_code = 'D2I'; snr_code = 'S2I';
    f_carrier = 1561.098e6; % BDS B1I 频率
else
    % 更新错误消息以反映新支持的系统
    error('不支持的卫星系统: %s。目前只支持 GPS (G), BeiDou (C), Galileo (E) 和 QZSS (J)。', target_satellite_id);
end
lambda_carrier = c / f_carrier; % 计算载波波长
% =========================【修改结束】=========================

% --- 初始化结果存储向量 ---
num_epochs_to_process = length(epoch_range);
time_vector = NaT(1, num_epochs_to_process);
phase_vector_rad = NaN(1, num_epochs_to_process);
processed_pr_vector = NaN(1, num_epochs_to_process);
amplitude_vector = NaN(1, num_epochs_to_process);
raw_snr_vector = NaN(1, num_epochs_to_process); % <--- 新增: 用于存储原始SNR
valid_results_count = 0;
initial_phase_anchor_cycles = NaN;
accumulated_residual_phase = 0;
last_epoch_time = NaT;

% --- 平滑器状态变量 ---
smooth_receiver_pos = []; % 平滑后的位置
smooth_receiver_clk_err = NaN; % 平滑后的钟差
alpha_pos = 0.1; % 位置平滑因子
alpha_clk = 0.1; % 钟差平滑因子

fprintf('--> 开始GPSense终极分析流程，目标卫星 %s ...\n', target_satellite_id);

% --- 2. 核心处理循环 (逐历元进行信号重构) ---
for loop_idx = 1:num_epochs_to_process
    epoch_idx = epoch_range(loop_idx);
    current_time = obs_data(epoch_idx).time;
    try
        % 步骤2.1: 获取该时刻的接收机位置和原始钟差 (这是所有计算的基础)
        % !! 依赖: 必须使用支持 G/C/E/J 的 calculate_receiver_position.m !!
        [receiver_pos_raw, receiver_clk_err_raw, sat_states] = ...
            calculate_receiver_position(obs_data, nav_data, epoch_idx);
        
        % 步骤2.2: 数据有效性检查
        if isempty(receiver_pos_raw) || any(isnan(receiver_pos_raw)) || isnan(receiver_clk_err_raw), continue; end
        if ~isfield(sat_states, target_satellite_id), continue; end
        sat_obs = obs_data(epoch_idx).data.(target_satellite_id);
        
        % 检查 (G/E/J 使用 C1C/D1C/S1C, C 使用 C2I/D2I/S2I)
        if ~(isfield(sat_obs, 'pseudorange') && isfield(sat_obs.pseudorange, pr_code) && ...
             isfield(sat_obs, 'doppler') && isfield(sat_obs.doppler, dop_code) && ...
             isfield(sat_obs, 'snr') && isfield(sat_obs.snr, snr_code))
            continue;
        end
        P_raw = sat_obs.pseudorange.(pr_code);
        measured_doppler = sat_obs.doppler.(dop_code);
        cn0_dbhz = sat_obs.snr.(snr_code);
        if isnan(P_raw) || isnan(measured_doppler) || isnan(cn0_dbhz), continue; end
        
        % --- 在线平滑位置与钟差 ---
        if isempty(smooth_receiver_pos)
            smooth_receiver_pos = receiver_pos_raw;
            smooth_receiver_clk_err = receiver_clk_err_raw;
        else
            smooth_receiver_pos = alpha_pos * receiver_pos_raw + (1 - alpha_pos) * smooth_receiver_pos;
            smooth_receiver_clk_err = alpha_clk * receiver_clk_err_raw + (1 - alpha_clk) * smooth_receiver_clk_err;
        end

        % 步骤2.3: 重构信号幅度 (依据论文公式4)
        current_amplitude = sqrt(Pn_noise_power * 10^(cn0_dbhz / 10));
        
        % 步骤2.4: 重构伪距 (使用平滑后的位置和钟差)
        sat_state = sat_states.(target_satellite_id);
        sat_pos = sat_state.position;
        geom_range = norm(sat_pos - smooth_receiver_pos);
        clock_corrected_pr = P_raw - c * (smooth_receiver_clk_err - sat_state.clock_error);
        processed_pr = clock_corrected_pr - geom_range;
        
        % 步骤2.5: 重构相位 (使用平滑值进行首次锚定)
        los_vec = (sat_pos - smooth_receiver_pos) / norm(sat_pos - smooth_receiver_pos);
        fD_geometric = -(sat_state.velocity' * los_vec) / lambda_carrier;
        residual_doppler = measured_doppler - fD_geometric;
        
        if isnan(initial_phase_anchor_cycles)
            initial_phase_anchor_cycles = (P_raw - c * (smooth_receiver_clk_err - sat_state.clock_error)) / lambda_carrier;
            accumulated_residual_phase = 0;
        end
        
        if ~isnat(last_epoch_time)
            dt = seconds(current_time - last_epoch_time);
            accumulated_residual_phase = accumulated_residual_phase + residual_doppler * dt;
        end
        total_phase_cycles = initial_phase_anchor_cycles + accumulated_residual_phase;
        
        % 步骤2.6: 存储当前历元的有效结果
        valid_results_count = valid_results_count + 1;
        time_vector(valid_results_count) = current_time;
        phase_vector_rad(valid_results_count) = total_phase_cycles * (2 * pi);
        processed_pr_vector(valid_results_count) = processed_pr;
        amplitude_vector(valid_results_count) = current_amplitude;
        raw_snr_vector(valid_results_count) = cn0_dbhz; % <--- 新增: 存储SNR值
        
        last_epoch_time = current_time;
    catch ME
        % 可以在此添加错误处理，例如: fprintf('处理历元 %d (%s) 时出错: %s\n', epoch_idx, target_satellite_id, ME.message);
    end
end
fprintf('--> 核心数据处理循环完成！\n\n');

% --- 3. 数据修剪与后处理 ---
if valid_results_count == 0, error('未能计算出任何有效结果，程序终止。'); end
time_vector = time_vector(1:valid_results_count);
phase_vector_rad = phase_vector_rad(1:valid_results_count);
processed_pr_vector = processed_pr_vector(1:valid_results_count);
amplitude_vector = amplitude_vector(1:valid_results_count);
raw_snr_vector = raw_snr_vector(1:valid_results_count); % <--- 新增: 修剪SNR向量

% 转换时间坐标为北京时间 (UTC+8)，并减去GPS时与UTC时的18秒闰秒差
time_vector_cst = time_vector + hours(8) - seconds(20);
% time_vector_cst = time_vector - seconds(18) + hours(8);
unwrapped_phase = unwrap(phase_vector_rad);

% --- 3.5 核心计算分离 (确保变量总是被定义) ---
% 归一化幅度计算
min_amp = min(amplitude_vector);
max_amp = max(amplitude_vector);
if (max_amp - min_amp) > eps
    amplitude_norm = (amplitude_vector - min_amp) / (max_amp - min_amp);
else
    amplitude_norm = ones(size(amplitude_vector)) * 0.5; % 避免除以零
end

% 平滑幅度计算 (使用Savitzky-Golay滤波器)
if length(amplitude_norm) > 11
    order = 2;
    framelen = 11;
    amplitude_norm_smooth = sgolayfilt(amplitude_norm, order, framelen);
else
    amplitude_norm_smooth = amplitude_norm; % 如果数据不足，则不平滑
end

% === [新增] SNR 平滑计算 ===
if length(raw_snr_vector) > 11
    % 使用 Savitzky-Golay 滤波器对 SNR 进行去噪 (参数与 step1 建议保持一致)
    % Order=2, FrameLen=7 (约0.7s窗口)
    smooth_snr_vector = sgolayfilt(raw_snr_vector, 2, 7);
else
    smooth_snr_vector = raw_snr_vector;
end
% =========================

% 去趋势伪距计算
window_size_pr = 25;
pr_trend = movmean(processed_pr_vector, window_size_pr, 'omitnan');
pseudorange_detrend = processed_pr_vector - pr_trend;

% 去趋势相位计算
window_size_phase = 25;
phase_trend_movmean = movmean(unwrapped_phase, window_size_phase, 'omitnan');
phase_detrend_movmean = unwrapped_phase - phase_trend_movmean;

% =========================================================================
% ======================= 4. 可视化模块 (按需注释) =========================
% =========================================================================

% --- 绘图模块 1: 未归一化的原始幅度 ---
% 作用: 显示信号强度的绝对物理量，单位为 sqrt(Watts)。
% 如何禁用: 在下面 "figure" 到 "datetick" 的所有行前加上 '%' 符号。
% -------------------------------------------------------------------------
fprintf('--> 正在处理和绘制: 未归一化的原始幅度...\n');
% figure('Name', sprintf('卫星 %s - 原始幅度', target_satellite_id), 'NumberTitle', 'off');
% plot(time_vector_cst, amplitude_vector, 'Color', [0.8500 0.3250 0.0980]); % 使用橙色
% title(sprintf('卫星 %s - 原始幅度 (未归一化)', target_satellite_id));
% xlabel('时间 (北京时间)');
% ylabel('原始幅度 (sqrt(W))');
% grid on;
% datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');


% ========================= 原始SNR绘图模块 =========================
% --- 绘图模块 1.5: 原始SNR观测值 ---
% 作用: 显示从RINEX文件中直接解析出的原始载噪比(C/N0)值。
% 如何禁用: 在下面 "figure" 到 "datetick" 的所有行前加上 '%' 符号。
% -------------------------------------------------------------------------
% fprintf('--> 正在处理和绘制: 原始SNR观测值...\n');
% figure('Name', sprintf('卫星 %s - 原始SNR', target_satellite_id), 'NumberTitle', 'off');
% plot(time_vector_cst, raw_snr_vector, '.-');
% ylim([30 50]);
% title(sprintf('卫星 %s - 原始SNR观测值 (C/N0)', target_satellite_id));
% xlabel('时间 (北京时间)');
% ylabel('SNR (dB-Hz)');
% grid on;
% datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');
% ========================= 原始SNR绘图模块结束 =========================


% ========================= [新增] SNR滤波对比绘图模块 =========================
% --- 绘图模块 1.6: SNR 滤波前后对比图 ---
% 作用: 对比显示原始SNR和经过SG滤波后的SNR，直观评估去噪效果。
% 如何禁用: 在下面 "figure" 到 "datetick" 的所有行前加上 '%' 符号。
% -------------------------------------------------------------------------
fprintf('--> 正在处理和绘制: SNR 滤波前后对比图...\n');
figure('Name', sprintf('卫星 %s - SNR 滤波对比', target_satellite_id), 'NumberTitle', 'off');

% 子图 1: 原始 SNR (默认颜色)
subplot(2, 1, 1);
plot(time_vector_cst, raw_snr_vector); 
title(sprintf('原始 SNR (C/N0) - 卫星 %s', target_satellite_id));
ylabel('SNR (dB-Hz)');
grid on;
datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');

% 子图 2: 滤波后 SNR (默认颜色)
subplot(2, 1, 2);
plot(time_vector_cst, smooth_snr_vector);
title(sprintf('滤波后 SNR (Savitzky-Golay) - 卫星 %s', target_satellite_id));
xlabel('时间 (北京时间)');
ylabel('SNR (dB-Hz)');
grid on;
datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');
% ========================= SNR滤波对比绘图模块结束 =========================


% --- 绘图模块 2: 归一化幅度 ---
% 作用: 显示信号强度的相对变化，对应论文中的幅度分析。
% 如何禁用: 在下面 "figure" 到 "datetick" 的所有行前加上 '%' 符号。
% -------------------------------------------------------------------------
fprintf('--> 正在处理和绘制: 归一化幅度...\n');
% figure('Name', sprintf('卫星 %s - 归一化幅度', target_satellite_id), 'NumberTitle', 'off');
% plot(time_vector_cst, amplitude_norm);
% title(sprintf('卫星 %s - 归一化幅度 (原始离散值)', target_satellite_id));
% xlabel('时间 (北京时间)');
% ylabel('归一化幅度 (无单位)');
% grid on;
% datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');


% --- 绘图模块 2.5: 平滑处理后的归一化幅度 ---
% 作用: 使用Savitzky-Golay滤波器对离散的归一化幅度进行平滑处理，使其更适合进行算法分析。
% 如何禁用: 在下面 "figure" 到 "datetick" 的所有行前加上 '%' 符号。
% -------------------------------------------------------------------------
fprintf('--> 正在处理和绘制: 平滑处理后的归一化幅度...\n');
% figure('Name', sprintf('卫星 %s - 平滑后的归一化幅度', target_satellite_id), 'NumberTitle', 'off');
% plot(time_vector_cst, amplitude_norm_smooth);
% title(sprintf('卫星 %s - 平滑后的归一化幅度 (Savitzky-Golay)', target_satellite_id));
% xlabel('时间 (北京时间)');
% ylabel('归一化幅度 (无单位)');
% grid on;
% datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');


% --- 绘图模块 3.5: 修正后但【未去趋势】的伪距 ---
% 作用: 独立绘制经过位置和钟差双重平滑修正后的伪距残差，即 processed_pr_vector。
% 如何禁用: 在下面 "figure" 到 "datetick" 的所有行前加上 '%' 符号。
% -------------------------------------------------------------------------
fprintf('--> 正在绘制修正后但【未去趋势】的伪距图...\n');
% figure('Name', sprintf('卫星 %s - 修正后未去趋势伪距', target_satellite_id), 'NumberTitle', 'off');
% plot(time_vector_cst, processed_pr_vector, 'b-');
% title(sprintf('卫星 %s - 修正后未去趋势伪距 (位置与钟差平滑)', target_satellite_id));
% xlabel('时间 (北京时间)');
% ylabel('残余伪距 (米)');
% grid on;
% datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');


% --- 绘图模块 3.6: 去趋势前后伪距对比图 ---
% 作用: 在同一张图中对比显示修正后的伪距残差，以及对其进行滑动平均去趋势后的结果。
% 如何禁用: 在下面 "figure" 到 "datetick" 的所有行前加上 '%' 符号。
% -------------------------------------------------------------------------
fprintf('--> 正在绘制去趋势前后伪距对比图...\n');
% figure('Name', sprintf('卫星 %s - 去趋势前后伪距对比', target_satellite_id), 'NumberTitle', 'off');
% plot(time_vector_cst, processed_pr_vector, 'Color', [0.7 0.7 0.7], 'DisplayName', '修正后、未去趋势伪距');
% hold on;
% plot(time_vector_cst, pseudorange_detrend, 'b-', 'LineWidth', 1.5, 'DisplayName', '最终去趋势伪距');
% hold off;
% title(sprintf('卫星 %s - 去趋势前后伪距对比', target_satellite_id));
% xlabel('时间 (北京时间)');
% ylabel('残余伪距 (米)');
% legend('show');
% grid on;
% datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');


% ========================= 原始三合一绘图模块开始 =========================
fprintf('--> 正在生成三合一对比图 (幅度/伪距/相位)...\n');
% figure('Name', sprintf('卫星 %s - 三合一分析图 (幅度/伪距/相位)', target_satellite_id), 'NumberTitle', 'off', 'Position', [100, 100, 800, 900]);
% % 子图 1: 平滑处理后的归一化幅度
% subplot(3, 1, 1);
% plot(time_vector_cst, amplitude_norm_smooth);
% title(sprintf('平滑后的归一化幅度 (Savitzky-Golay) - 卫星 %s', target_satellite_id));
% xlabel('时间 (北京时间)');
% ylabel('归一化幅度');
% grid on;
% datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');
% % 子图 2: 去趋势后的伪距
% subplot(3, 1, 2);
% plot(time_vector_cst, pseudorange_detrend, 'b-');
% title('去趋势后的伪距 (已进行在线钟差与位置平滑)');
% xlabel('时间 (北京时间)');
% ylabel('残余伪距 (米)');
% grid on;
% datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');
% % 子图 3: 去趋势相位 (滑动平均法)
% subplot(3, 1, 3);
% plot(time_vector_cst, phase_detrend_movmean, 'r-');
% title('去趋势GPSense相位 (滑动平均法)');
% xlabel('时间 (北京时间)');
% ylabel('去趋势相位 (弧度)');
% grid on;
% datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');
% ========================= 原始三合一绘图模块结束 =========================


% --- 绘图模块 5: 相位去趋势 (方法B - 带通滤波法) ---
% 作用: 使用带通滤波器精确提取人体活动频段(0.2-5.0Hz)的相位变化，是更专业的处理方法。
% 如何禁用: 在下面 "if" 到 "datetick" 的所有行前加上 '%' 符号。
% -------------------------------------------------------------------------
fprintf('--> 正在处理和绘制: 相位去趋势 (带通滤波法)...\n');
% 查找并替换 unwrapped_phase 中的 NaN 和 Inf
valid_indices = isfinite(unwrapped_phase);
if length(time_vector) > 10 && sum(valid_indices) > 20 % 确保有足够的有效数据点
    mean_dt = mean(seconds(diff(time_vector)));
    Fs = 1 / mean_dt;
    fprintf('    计算出的平均采样率为: %.2f Hz\n', Fs);
    f_low = 0.2;
    f_high = 3;
    order = 4;
    [b, a] = butter(order, [f_low f_high] / (Fs / 2), 'bandpass');
    % 在滤波前，对 NaN 进行线性插值处理
    phase_to_filter = unwrapped_phase;
    phase_to_filter(~valid_indices) = interp1(find(valid_indices), unwrapped_phase(valid_indices), find(~valid_indices), 'linear', 'extrap');
    phase_detrend_bpf = filtfilt(b, a, phase_to_filter);
%     figure('Name', sprintf('卫星 %s - 带通滤波后相位', target_satellite_id), 'NumberTitle', 'off');
%     plot(time_vector_cst, phase_detrend_bpf, 'm-');
%     title(sprintf('卫星 %s - 带通滤波后的GPSense相位 (0.2-3.0 Hz)', target_satellite_id));
%     xlabel('时间 (北京时间)');
%     ylabel('滤波后相位 (弧度)');
%     grid on;
%     datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');
else
    phase_detrend_bpf = NaN(size(time_vector));
    fprintf('    数据点不足或含有过多无效值，跳过带通滤波处理。\n');
end


% ========================= 新增SNR三合一绘图模块开始 =========================
% --- 绘图模块 6: 三合一对比图 (SNR/幅度/相位) ---
% 作用: 在同一个图窗中对比显示原始SNR、平滑后的归一化幅度、以及去趋势后的相位。
% 如何禁用: 在下面 "figure" 到 "datetick" 的所有行前加上 '%' 符号。
% -------------------------------------------------------------------------
% fprintf('--> 正在生成三合一对比图 (SNR/幅度/相位)...\n');
% figure('Name', sprintf('卫星 %s - 三合一分析图 (SNR/幅度/相位)', target_satellite_id), 'NumberTitle', 'off', 'Position', [200, -10, 800, 900]); % 调整位置避免重叠
% % 子图 1: 原始SNR
% subplot(3, 1, 1);
% plot(time_vector_cst, raw_snr_vector, '.-');
% ylim([30 50]); % 固定Y轴范围
% title(sprintf('原始SNR观测值 (C/N0) - 卫星 %s', target_satellite_id));
% xlabel('时间 (北京时间)');
% ylabel('SNR (dB-Hz)');
% grid on;
% datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits'); % <-- 使用正确的时间格式
% % 子图 2: 平滑处理后的归一化幅度
% subplot(3, 1, 2);
% plot(time_vector_cst, amplitude_norm_smooth);
% title('归一化幅度');
% xlabel('时间 (北京时间)');
% ylabel('归一化幅度');
% grid on;
% datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits'); % <-- 使用正确的时间格式
% % 子图 3: 去趋势相位 (滑动平均法)
% subplot(3, 1, 3);
% plot(time_vector_cst, phase_detrend_movmean, 'r-'); % <-- 已将颜色改回 'r-'
% title('去趋势GPSense相位');
% xlabel('时间 (北京时间)');
% ylabel('去趋势相位 (弧度)');
% grid on;
% datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits'); % <-- 使用正确的时间格式
% ========================= 新增SNR三合一绘图模块结束 =========================

fprintf('\n✅ 所有分析和绘图任务完成！\n');
end







% % ============== analyze_gpsense_ultimate.m (v5.6 - 已扩展支持 G/C/E/J) ==============
% %
% % 功能:
% %   本脚本是GPSense分析流程的终极整合版本。它集成了信号幅度(归一化与原始值)、
% %   伪距和两种核心相位去趋势方法（滑动平均法 和 带通滤波法）的计算与可视化。
% %   【V5.6更新】: 同步添加包含原始SNR的三合一绘图模块（时间轴版本）。
% %   【多系统更新】: 扩展以支持 Galileo (E) 和 QZSS (J) 卫星。
% %
% % 输出:
% %   time_vector      - 有效数据点的时间向量 (datetime)
% %   amplitude_norm   - 归一化后的信号幅度 (0-1范围)
% %   pseudorange_detrend - 去趋势后的伪距 (米)
% %   phase_detrend_movmean - 通过滑动平均法去趋势后的相位 (弧度)
% %   phase_detrend_bpf     - 通过带通滤波法处理后的相位 (弧度)
% %
% function [time_vector, amplitude_norm, pseudorange_detrend, phase_detrend_movmean, phase_detrend_bpf] = analyze_gpsense_ultimate(obs_data, nav_data, target_satellite_id, epoch_range)
% 
% % --- 1. 参数检查与初始化 ---
% if nargin < 3, error('请提供目标卫星ID (例如 ''G01'')。'); end
% if nargin < 4 || isempty(epoch_range), epoch_range = 1:length(obs_data); end
% 
% % --- 物理常数 ---
% c = 299792458.0; % 光速 (m/s)
% k_boltzmann = 1.380649e-23; % 玻尔兹曼常数 (J/K)
% T_noise_kelvin = 290.0; % 系统噪声温度 (K), 标准室温假设
% Pn_noise_power = k_boltzmann * T_noise_kelvin; % 噪声功率
% 
% % ====================【核心修改：扩展卫星系统支持】====================
% % --- 动态适配卫星系统 (G/C/E/J) ---
% sys_char = upper(target_satellite_id(1));
% if sys_char == 'G' || sys_char == 'E' || sys_char == 'J'
%     % GPS (G), Galileo (E), QZSS (J) 均使用 L1/E1 频段
%     pr_code = 'C1C'; dop_code = 'D1C'; snr_code = 'S1C';
%     f_carrier = 1575.42e6; % GPS L1 / GAL E1 / QZS L1 频率
% elseif sys_char == 'C'
%     % BeiDou (C) 使用 B1I 频段
%     pr_code = 'C2I'; dop_code = 'D2I'; snr_code = 'S2I';
%     f_carrier = 1561.098e6; % BDS B1I 频率
% else
%     % 更新错误消息以反映新支持的系统
%     error('不支持的卫星系统: %s。目前只支持 GPS (G), BeiDou (C), Galileo (E) 和 QZSS (J)。', target_satellite_id);
% end
% lambda_carrier = c / f_carrier; % 计算载波波长
% % =========================【修改结束】=========================
% 
% % --- 初始化结果存储向量 ---
% num_epochs_to_process = length(epoch_range);
% time_vector = NaT(1, num_epochs_to_process);
% phase_vector_rad = NaN(1, num_epochs_to_process);
% processed_pr_vector = NaN(1, num_epochs_to_process);
% amplitude_vector = NaN(1, num_epochs_to_process);
% raw_snr_vector = NaN(1, num_epochs_to_process); % <--- 新增: 用于存储原始SNR
% valid_results_count = 0;
% initial_phase_anchor_cycles = NaN;
% accumulated_residual_phase = 0;
% last_epoch_time = NaT;
% 
% % --- 平滑器状态变量 ---
% smooth_receiver_pos = []; % 平滑后的位置
% smooth_receiver_clk_err = NaN; % 平滑后的钟差
% alpha_pos = 0.1; % 位置平滑因子
% alpha_clk = 0.1; % 钟差平滑因子
% 
% fprintf('--> 开始GPSense终极分析流程，目标卫星 %s ...\n', target_satellite_id);
% 
% % --- 2. 核心处理循环 (逐历元进行信号重构) ---
% for loop_idx = 1:num_epochs_to_process
%     epoch_idx = epoch_range(loop_idx);
%     current_time = obs_data(epoch_idx).time;
%     try
%         % 步骤2.1: 获取该时刻的接收机位置和原始钟差 (这是所有计算的基础)
%         % !! 依赖: 必须使用支持 G/C/E/J 的 calculate_receiver_position.m !!
%         [receiver_pos_raw, receiver_clk_err_raw, sat_states] = ...
%             calculate_receiver_position(obs_data, nav_data, epoch_idx);
%         
%         % 步骤2.2: 数据有效性检查
%         if isempty(receiver_pos_raw) || any(isnan(receiver_pos_raw)) || isnan(receiver_clk_err_raw), continue; end
%         if ~isfield(sat_states, target_satellite_id), continue; end
%         sat_obs = obs_data(epoch_idx).data.(target_satellite_id);
%         
%         % 检查 (G/E/J 使用 C1C/D1C/S1C, C 使用 C2I/D2I/S2I)
%         if ~(isfield(sat_obs, 'pseudorange') && isfield(sat_obs.pseudorange, pr_code) && ...
%              isfield(sat_obs, 'doppler') && isfield(sat_obs.doppler, dop_code) && ...
%              isfield(sat_obs, 'snr') && isfield(sat_obs.snr, snr_code))
%             continue;
%         end
%         P_raw = sat_obs.pseudorange.(pr_code);
%         measured_doppler = sat_obs.doppler.(dop_code);
%         cn0_dbhz = sat_obs.snr.(snr_code);
%         if isnan(P_raw) || isnan(measured_doppler) || isnan(cn0_dbhz), continue; end
%         
%         % --- 在线平滑位置与钟差 ---
%         if isempty(smooth_receiver_pos)
%             smooth_receiver_pos = receiver_pos_raw;
%             smooth_receiver_clk_err = receiver_clk_err_raw;
%         else
%             smooth_receiver_pos = alpha_pos * receiver_pos_raw + (1 - alpha_pos) * smooth_receiver_pos;
%             smooth_receiver_clk_err = alpha_clk * receiver_clk_err_raw + (1 - alpha_clk) * smooth_receiver_clk_err;
%         end
% 
%         % 步骤2.3: 重构信号幅度 (依据论文公式4)
%         current_amplitude = sqrt(Pn_noise_power * 10^(cn0_dbhz / 10));
%         
%         % 步骤2.4: 重构伪距 (使用平滑后的位置和钟差)
%         sat_state = sat_states.(target_satellite_id);
%         sat_pos = sat_state.position;
%         geom_range = norm(sat_pos - smooth_receiver_pos);
%         clock_corrected_pr = P_raw - c * (smooth_receiver_clk_err - sat_state.clock_error);
%         processed_pr = clock_corrected_pr - geom_range;
%         
%         % 步骤2.5: 重构相位 (使用平滑值进行首次锚定)
%         los_vec = (sat_pos - smooth_receiver_pos) / norm(sat_pos - smooth_receiver_pos);
%         fD_geometric = -(sat_state.velocity' * los_vec) / lambda_carrier;
%         residual_doppler = measured_doppler - fD_geometric;
%         
%         if isnan(initial_phase_anchor_cycles)
%             initial_phase_anchor_cycles = (P_raw - c * (smooth_receiver_clk_err - sat_state.clock_error)) / lambda_carrier;
%             accumulated_residual_phase = 0;
%         end
%         
%         if ~isnat(last_epoch_time)
%             dt = seconds(current_time - last_epoch_time);
%             accumulated_residual_phase = accumulated_residual_phase + residual_doppler * dt;
%         end
%         total_phase_cycles = initial_phase_anchor_cycles + accumulated_residual_phase;
%         
%         % 步骤2.6: 存储当前历元的有效结果
%         valid_results_count = valid_results_count + 1;
%         time_vector(valid_results_count) = current_time;
%         phase_vector_rad(valid_results_count) = total_phase_cycles * (2 * pi);
%         processed_pr_vector(valid_results_count) = processed_pr;
%         amplitude_vector(valid_results_count) = current_amplitude;
%         raw_snr_vector(valid_results_count) = cn0_dbhz; % <--- 新增: 存储SNR值
%         
%         last_epoch_time = current_time;
%     catch ME
%         % 可以在此添加错误处理，例如: fprintf('处理历元 %d (%s) 时出错: %s\n', epoch_idx, target_satellite_id, ME.message);
%     end
% end
% fprintf('--> 核心数据处理循环完成！\n\n');
% 
% % --- 3. 数据修剪与后处理 ---
% if valid_results_count == 0, error('未能计算出任何有效结果，程序终止。'); end
% time_vector = time_vector(1:valid_results_count);
% phase_vector_rad = phase_vector_rad(1:valid_results_count);
% processed_pr_vector = processed_pr_vector(1:valid_results_count);
% amplitude_vector = amplitude_vector(1:valid_results_count);
% raw_snr_vector = raw_snr_vector(1:valid_results_count); % <--- 新增: 修剪SNR向量
% 
% % 转换时间坐标为北京时间 (UTC+8)，并减去GPS时与UTC时的18秒闰秒差
% time_vector_cst = time_vector + hours(8) - seconds(20);
% % time_vector_cst = time_vector - seconds(18) + hours(8);
% unwrapped_phase = unwrap(phase_vector_rad);
% 
% % --- 3.5 核心计算分离 (确保变量总是被定义) ---
% % 归一化幅度计算
% min_amp = min(amplitude_vector);
% max_amp = max(amplitude_vector);
% if (max_amp - min_amp) > eps
%     amplitude_norm = (amplitude_vector - min_amp) / (max_amp - min_amp);
% else
%     amplitude_norm = ones(size(amplitude_vector)) * 0.5; % 避免除以零
% end
% 
% % 平滑幅度计算 (使用Savitzky-Golay滤波器)
% if length(amplitude_norm) > 11
%     order = 2;
%     framelen = 11;
%     amplitude_norm_smooth = sgolayfilt(amplitude_norm, order, framelen);
% else
%     amplitude_norm_smooth = amplitude_norm; % 如果数据不足，则不平滑
% end
% 
% % 去趋势伪距计算
% window_size_pr = 25;
% pr_trend = movmean(processed_pr_vector, window_size_pr, 'omitnan');
% pseudorange_detrend = processed_pr_vector - pr_trend;
% 
% % 去趋势相位计算
% window_size_phase = 25;
% phase_trend_movmean = movmean(unwrapped_phase, window_size_phase, 'omitnan');
% phase_detrend_movmean = unwrapped_phase - phase_trend_movmean;
% 
% % =========================================================================
% % ======================= 4. 可视化模块 (按需注释) =========================
% % =========================================================================
% 
% % --- 绘图模块 1: 未归一化的原始幅度 ---
% % 作用: 显示信号强度的绝对物理量，单位为 sqrt(Watts)。
% % 如何禁用: 在下面 "figure" 到 "datetick" 的所有行前加上 '%' 符号。
% % -------------------------------------------------------------------------
% fprintf('--> 正在处理和绘制: 未归一化的原始幅度...\n');
% % figure('Name', sprintf('卫星 %s - 原始幅度', target_satellite_id), 'NumberTitle', 'off');
% % plot(time_vector_cst, amplitude_vector, 'Color', [0.8500 0.3250 0.0980]); % 使用橙色
% % title(sprintf('卫星 %s - 原始幅度 (未归一化)', target_satellite_id));
% % xlabel('时间 (北京时间)');
% % ylabel('原始幅度 (sqrt(W))');
% % grid on;
% % datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');
% 
% 
% % ========================= 原始SNR绘图模块 =========================
% % --- 绘图模块 1.5: 原始SNR观测值 ---
% % 作用: 显示从RINEX文件中直接解析出的原始载噪比(C/N0)值。
% % 如何禁用: 在下面 "figure" 到 "datetick" 的所有行前加上 '%' 符号。
% % -------------------------------------------------------------------------
% fprintf('--> 正在处理和绘制: 原始SNR观测值...\n');
% figure('Name', sprintf('卫星 %s - 原始SNR', target_satellite_id), 'NumberTitle', 'off');
% plot(time_vector_cst, raw_snr_vector, '.-');
% ylim([30 50]);
% title(sprintf('卫星 %s - 原始SNR观测值 (C/N0)', target_satellite_id));
% xlabel('时间 (北京时间)');
% ylabel('SNR (dB-Hz)');
% grid on;
% datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');
% % ========================= 原始SNR绘图模块结束 =========================
% 
% 
% % --- 绘图模块 2: 归一化幅度 ---
% % 作用: 显示信号强度的相对变化，对应论文中的幅度分析。
% % 如何禁用: 在下面 "figure" 到 "datetick" 的所有行前加上 '%' 符号。
% % -------------------------------------------------------------------------
% fprintf('--> 正在处理和绘制: 归一化幅度...\n');
% % figure('Name', sprintf('卫星 %s - 归一化幅度', target_satellite_id), 'NumberTitle', 'off');
% % plot(time_vector_cst, amplitude_norm);
% % title(sprintf('卫星 %s - 归一化幅度 (原始离散值)', target_satellite_id));
% % xlabel('时间 (北京时间)');
% % ylabel('归一化幅度 (无单位)');
% % grid on;
% % datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');
% 
% 
% % --- 绘图模块 2.5: 平滑处理后的归一化幅度 ---
% % 作用: 使用Savitzky-Golay滤波器对离散的归一化幅度进行平滑处理，使其更适合进行算法分析。
% % 如何禁用: 在下面 "figure" 到 "datetick" 的所有行前加上 '%' 符号。
% % -------------------------------------------------------------------------
% fprintf('--> 正在处理和绘制: 平滑处理后的归一化幅度...\n');
% % figure('Name', sprintf('卫星 %s - 平滑后的归一化幅度', target_satellite_id), 'NumberTitle', 'off');
% % plot(time_vector_cst, amplitude_norm_smooth);
% % title(sprintf('卫星 %s - 平滑后的归一化幅度 (Savitzky-Golay)', target_satellite_id));
% % xlabel('时间 (北京时间)');
% % ylabel('归一化幅度 (无单位)');
% % grid on;
% % datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');
% 
% 
% % --- 绘图模块 3.5: 修正后但【未去趋势】的伪距 ---
% % 作用: 独立绘制经过位置和钟差双重平滑修正后的伪距残差，即 processed_pr_vector。
% % 如何禁用: 在下面 "figure" 到 "datetick" 的所有行前加上 '%' 符号。
% % -------------------------------------------------------------------------
% fprintf('--> 正在绘制修正后但【未去趋势】的伪距图...\n');
% % figure('Name', sprintf('卫星 %s - 修正后未去趋势伪距', target_satellite_id), 'NumberTitle', 'off');
% % plot(time_vector_cst, processed_pr_vector, 'b-');
% % title(sprintf('卫星 %s - 修正后未去趋势伪距 (位置与钟差平滑)', target_satellite_id));
% % xlabel('时间 (北京时间)');
% % ylabel('残余伪距 (米)');
% % grid on;
% % datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');
% 
% 
% % --- 绘图模块 3.6: 去趋势前后伪距对比图 ---
% % 作用: 在同一张图中对比显示修正后的伪距残差，以及对其进行滑动平均去趋势后的结果。
% % 如何禁用: 在下面 "figure" 到 "datetick" 的所有行前加上 '%' 符号。
% % -------------------------------------------------------------------------
% fprintf('--> 正在绘制去趋势前后伪距对比图...\n');
% % figure('Name', sprintf('卫星 %s - 去趋势前后伪距对比', target_satellite_id), 'NumberTitle', 'off');
% % plot(time_vector_cst, processed_pr_vector, 'Color', [0.7 0.7 0.7], 'DisplayName', '修正后、未去趋势伪距');
% % hold on;
% % plot(time_vector_cst, pseudorange_detrend, 'b-', 'LineWidth', 1.5, 'DisplayName', '最终去趋势伪距');
% % hold off;
% % title(sprintf('卫星 %s - 去趋势前后伪距对比', target_satellite_id));
% % xlabel('时间 (北京时间)');
% % ylabel('残余伪距 (米)');
% % legend('show');
% % grid on;
% % datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');
% 
% 
% % ========================= 原始三合一绘图模块开始 =========================
% fprintf('--> 正在生成三合一对比图 (幅度/伪距/相位)...\n');
% % figure('Name', sprintf('卫星 %s - 三合一分析图 (幅度/伪距/相位)', target_satellite_id), 'NumberTitle', 'off', 'Position', [100, 100, 800, 900]);
% % % 子图 1: 平滑处理后的归一化幅度
% % subplot(3, 1, 1);
% % plot(time_vector_cst, amplitude_norm_smooth);
% % title(sprintf('平滑后的归一化幅度 (Savitzky-Golay) - 卫星 %s', target_satellite_id));
% % xlabel('时间 (北京时间)');
% % ylabel('归一化幅度');
% % grid on;
% % datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');
% % % 子图 2: 去趋势后的伪距
% % subplot(3, 1, 2);
% % plot(time_vector_cst, pseudorange_detrend, 'b-');
% % title('去趋势后的伪距 (已进行在线钟差与位置平滑)');
% % xlabel('时间 (北京时间)');
% % ylabel('残余伪距 (米)');
% % grid on;
% % datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');
% % % 子图 3: 去趋势相位 (滑动平均法)
% % subplot(3, 1, 3);
% % plot(time_vector_cst, phase_detrend_movmean, 'r-');
% % title('去趋势GPSense相位 (滑动平均法)');
% % xlabel('时间 (北京时间)');
% % ylabel('去趋势相位 (弧度)');
% % grid on;
% % datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');
% % ========================= 原始三合一绘图模块结束 =========================
% 
% 
% % --- 绘图模块 5: 相位去趋势 (方法B - 带通滤波法) ---
% % 作用: 使用带通滤波器精确提取人体活动频段(0.2-5.0Hz)的相位变化，是更专业的处理方法。
% % 如何禁用: 在下面 "if" 到 "datetick" 的所有行前加上 '%' 符号。
% % -------------------------------------------------------------------------
% fprintf('--> 正在处理和绘制: 相位去趋势 (带通滤波法)...\n');
% % 查找并替换 unwrapped_phase 中的 NaN 和 Inf
% valid_indices = isfinite(unwrapped_phase);
% if length(time_vector) > 10 && sum(valid_indices) > 20 % 确保有足够的有效数据点
%     mean_dt = mean(seconds(diff(time_vector)));
%     Fs = 1 / mean_dt;
%     fprintf('    计算出的平均采样率为: %.2f Hz\n', Fs);
%     f_low = 0.2;
%     f_high = 3;
%     order = 4;
%     [b, a] = butter(order, [f_low f_high] / (Fs / 2), 'bandpass');
%     % 在滤波前，对 NaN 进行线性插值处理
%     phase_to_filter = unwrapped_phase;
%     phase_to_filter(~valid_indices) = interp1(find(valid_indices), unwrapped_phase(valid_indices), find(~valid_indices), 'linear', 'extrap');
%     phase_detrend_bpf = filtfilt(b, a, phase_to_filter);
% %     figure('Name', sprintf('卫星 %s - 带通滤波后相位', target_satellite_id), 'NumberTitle', 'off');
% %     plot(time_vector_cst, phase_detrend_bpf, 'm-');
% %     title(sprintf('卫星 %s - 带通滤波后的GPSense相位 (0.2-3.0 Hz)', target_satellite_id));
% %     xlabel('时间 (北京时间)');
% %     ylabel('滤波后相位 (弧度)');
% %     grid on;
% %     datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');
% else
%     phase_detrend_bpf = NaN(size(time_vector));
%     fprintf('    数据点不足或含有过多无效值，跳过带通滤波处理。\n');
% end
% 
% 
% % ========================= 新增SNR三合一绘图模块开始 =========================
% % --- 绘图模块 6: 三合一对比图 (SNR/幅度/相位) ---
% % 作用: 在同一个图窗中对比显示原始SNR、平滑后的归一化幅度、以及去趋势后的相位。
% % 如何禁用: 在下面 "figure" 到 "datetick" 的所有行前加上 '%' 符号。
% % -------------------------------------------------------------------------
% % fprintf('--> 正在生成三合一对比图 (SNR/幅度/相位)...\n');
% % figure('Name', sprintf('卫星 %s - 三合一分析图 (SNR/幅度/相位)', target_satellite_id), 'NumberTitle', 'off', 'Position', [200, -10, 800, 900]); % 调整位置避免重叠
% % % 子图 1: 原始SNR
% % subplot(3, 1, 1);
% % plot(time_vector_cst, raw_snr_vector, '.-');
% % ylim([30 50]); % 固定Y轴范围
% % title(sprintf('原始SNR观测值 (C/N0) - 卫星 %s', target_satellite_id));
% % xlabel('时间 (北京时间)');
% % ylabel('SNR (dB-Hz)');
% % grid on;
% % datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits'); % <-- 使用正确的时间格式
% % % 子图 2: 平滑处理后的归一化幅度
% % subplot(3, 1, 2);
% % plot(time_vector_cst, amplitude_norm_smooth);
% % title('归一化幅度');
% % xlabel('时间 (北京时间)');
% % ylabel('归一化幅度');
% % grid on;
% % datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits'); % <-- 使用正确的时间格式
% % % 子图 3: 去趋势相位 (滑动平均法)
% % subplot(3, 1, 3);
% % plot(time_vector_cst, phase_detrend_movmean, 'r-'); % <-- 已将颜色改回 'r-'
% % title('去趋势GPSense相位');
% % xlabel('时间 (北京时间)');
% % ylabel('去趋势相位 (弧度)');
% % grid on;
% % datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits'); % <-- 使用正确的时间格式
% % ========================= 新增SNR三合一绘图模块结束 =========================
% 
% fprintf('\n✅ 所有分析和绘图任务完成！\n');
% end
