% 
% 
% % ============== analyze_and_plot_gpsense.m (新增幅度归一化，严格遵循论文) ==============
% function [time_vector, detrended_phase_rad, processed_pr_vector, amplitude_vector] = analyze_and_plot_gpsense(obs_data, nav_data, target_satellite_id, epoch_range)
% % ANALYZE_AND_PLOT_GPSENSE - 新增信号幅度重建和归一化，并使用北京时间绘图。
% 
% % --- 1. 参数检查与初始化 ---
% if nargin < 3, error('请提供目标卫星ID。'); end
% if nargin < 4 || isempty(epoch_range), epoch_range = 1:length(obs_data); end
% 
% % --- 物理常数 ---
% c = 299792458.0;
% f_L1 = 1575.42e6;
% lambda_L1 = c / f_L1;
% 
% %{
%  定义幅度计算所需常数 (据论文公式4)
%  Pn = k * T
%  k: 玻尔兹曼常数 (J/K)
%  T: 系统噪声温度 (K)。论文未指定，此处采用标准室温290K作为合理假设。
% %}
% k_boltzmann = 1.380649e-23;
% T_noise_kelvin = 290.0;
% Pn_noise_power = k_boltzmann * T_noise_kelvin;
% 
% % --- 初始化结果存储向量 ---
% num_epochs_to_process = length(epoch_range);
% time_vector = NaT(1, num_epochs_to_process);
% phase_vector_rad = NaN(1, num_epochs_to_process);
% processed_pr_vector = NaN(1, num_epochs_to_process);
% amplitude_vector = NaN(1, num_epochs_to_process);
% 
% valid_results_count = 0;
% initial_phase_anchor_cycles = NaN;
% accumulated_residual_phase = 0;
% last_epoch_time = NaT;
% 
% fprintf('--> 开始GPSense分析流程，历元范围 %d-%d，目标卫星 %s ...\n', epoch_range(1), epoch_range(end), target_satellite_id);
% 
% % --- 2. 遍历历元 ---
% for loop_idx = 1:num_epochs_to_process
%     epoch_idx = epoch_range(loop_idx);
%     current_time = obs_data(epoch_idx).time;
%     try
%         [receiver_pos, receiver_clk_err, sat_states] = ...
%             calculate_receiver_position(obs_data, nav_data, epoch_idx);
%         
%         if ~isfield(sat_states, target_satellite_id), continue; end
%         
%         % --- 提取当前历元所需的所有观测量 ---
%         sat_obs = obs_data(epoch_idx).data.(target_satellite_id);
%         
%         if ~(isfield(sat_obs, 'pseudorange') && isfield(sat_obs.pseudorange, 'C1C') && ...
%              isfield(sat_obs, 'doppler') && isfield(sat_obs.doppler, 'D1C') && ...
%              isfield(sat_obs, 'snr') && isfield(sat_obs.snr, 'S1C'))
%             continue;
%         end
%         
%         P_raw = sat_obs.pseudorange.C1C;
%         measured_doppler = sat_obs.doppler.D1C;
%         cn0_dbhz = sat_obs.snr.S1C;
%         
%         if isnan(P_raw) || isnan(measured_doppler) || isnan(cn0_dbhz), continue; end
% 
%         % 严格按照论文公式(4)计算幅度
%         current_amplitude = sqrt(Pn_noise_power * 10^(cn0_dbhz / 10));
%         
%         sat_state = sat_states.(target_satellite_id);
%         sat_pos = sat_state.position;
%         sat_vel = sat_state.velocity;
%         sat_clk_err = sat_state.clock_error;
%         
%         % --- 计算“处理后的伪距” ---
%         geom_range = norm(sat_pos - receiver_pos);
%         clock_corrected_pr = P_raw - c * (receiver_clk_err - sat_clk_err);
%         processed_pr = clock_corrected_pr - geom_range;
%         
%         % --- 计算“去趋势相位” ---
%         los_vec = (sat_pos - receiver_pos) / norm(sat_pos - receiver_pos);
%         fD_geometric = -(sat_vel' * los_vec) / lambda_L1;
%         residual_doppler = measured_doppler - fD_geometric;
%         
%         if isnan(initial_phase_anchor_cycles)
%             initial_phase_anchor_cycles = (P_raw - c * (receiver_clk_err - sat_clk_err)) / lambda_L1;
%             accumulated_residual_phase = 0;
%         end
%         
%         if ~isnat(last_epoch_time)
%             dt = seconds(current_time - last_epoch_time);
%             accumulated_residual_phase = accumulated_residual_phase + residual_doppler * dt;
%         end
%         total_phase_cycles = initial_phase_anchor_cycles + accumulated_residual_phase;
%         
%         % --- 存储所有结果 ---
%         valid_results_count = valid_results_count + 1;
%         time_vector(valid_results_count) = current_time;
%         phase_vector_rad(valid_results_count) = total_phase_cycles * (2 * pi);
%         processed_pr_vector(valid_results_count) = processed_pr;
%         amplitude_vector(valid_results_count) = current_amplitude;
%         
%         last_epoch_time = current_time;
%     catch ME
%         fprintf('在处理历元 %d 时发生错误: %s\n', epoch_idx, ME.message);
%     end
% end
% fprintf('--> 数据处理循环完成！\n\n');
% 
% % --- 3. 数据修剪与绘图 ---
% time_vector = time_vector(1:valid_results_count);
% phase_vector_rad = phase_vector_rad(1:valid_results_count);
% processed_pr_vector = processed_pr_vector(1:valid_results_count);
% amplitude_vector = amplitude_vector(1:valid_results_count);
% if isempty(time_vector), error('未能计算出任何有效结果。'); end
% 
% % 转换时间坐标为北京时间
% fprintf('--> 转换时间坐标为北京时间...\n');
% time_vector_cst = time_vector - seconds(18) + hours(8);
% 
% 
% % ====================【新增/修改：对幅度进行归一化并绘图】====================
% % 对幅度数据进行最小-最大归一化，将其缩放到 [0, 1] 区间
% fprintf('--> 对幅度数据进行归一化处理...\n');
% min_amp = min(amplitude_vector);
% max_amp = max(amplitude_vector);
% % 检查以避免除以零（当所有幅值都相同时）
% if (max_amp - min_amp) > eps 
%     normalized_amplitude_vector = (amplitude_vector - min_amp) / (max_amp - min_amp);
% else
%     % 如果所有值都相同，则将归一化结果设为0.5
%     normalized_amplitude_vector = ones(size(amplitude_vector)) * 0.5;
% end
% 
% figure('Name', sprintf('卫星 %s - 归一化幅度', target_satellite_id), 'NumberTitle', 'off');
% plot(time_vector_cst, normalized_amplitude_vector, 'g-'); % 使用归一化后的数据绘图
% title(sprintf('卫星 %s - 归一化幅度', target_satellite_id));
% xlabel('时间 (北京时间)');
% ylabel('归一化幅度 (无单位)');
% grid on;
% datetick('x', 'HH:MM:SS', 'keepticks', 'keeplimits');
% % ======================================================================
% 
% 
% % 绘制处理后的伪距
% window_size = 25; 
% pr_trend = movmean(processed_pr_vector, window_size, 'omitnan');
% detrended_pr_vector = processed_pr_vector - pr_trend;
% figure('Name', sprintf('卫星 %s - 去趋势后的伪距残差', target_satellite_id), 'NumberTitle', 'off');
% plot(time_vector_cst, detrended_pr_vector, 'b-');
% title(sprintf('卫星 %s - 去趋势后的伪距残差', target_satellite_id));
% xlabel('时间 (北京时间)');
% ylabel('残余伪距 (米)');
% grid on;
% datetick('x', 'HH:MM:SS', 'keepticks', 'keeplimits');
% 
% % 绘制去趋势相位
% unwrapped_phase = unwrap(phase_vector_rad);
% trend = movmean(unwrapped_phase, window_size, 'omitnan');
% detrended_phase_rad = unwrapped_phase - trend;
% figure('Name', sprintf('卫星 %s - 最终去趋势相位', target_satellite_id), 'NumberTitle', 'off');
% plot(time_vector_cst, detrended_phase_rad, 'r-');
% title(sprintf('卫星 %s - 最终去趋势GPSense相位', target_satellite_id));
% xlabel('时间 (北京时间)');
% ylabel('去趋势相位 (弧度)');
% grid on;
% datetick('x', 'HH:MM:SS', 'keepticks', 'keeplimits');
% 
% fprintf('✅ 图形绘制完成！\n');
% end











% ============== analyze_and_plot_gpsense.m (修正采样率计算错误) ==============
function [time_vector, detrended_phase_rad, processed_pr_vector, amplitude_vector] = analyze_and_plot_gpsense(obs_data, nav_data, target_satellite_id, epoch_range)
% ANALYZE_AND_PLOT_GPSENSE - 修正了采样率计算的错误。

% --- 1. 参数检查与初始化 ---
if nargin < 3, error('请提供目标卫星ID。'); end
if nargin < 4 || isempty(epoch_range), epoch_range = 1:length(obs_data); end

% --- 物理常数 ---
c = 299792458.0;
f_L1 = 1575.42e6;
lambda_L1 = c / f_L1;

% 幅度计算所需常数
k_boltzmann = 1.380649e-23;
T_noise_kelvin = 290.0;
Pn_noise_power = k_boltzmann * T_noise_kelvin;

% --- 初始化结果存储向量 ---
num_epochs_to_process = length(epoch_range);
time_vector = NaT(1, num_epochs_to_process);
phase_vector_rad = NaN(1, num_epochs_to_process);
processed_pr_vector = NaN(1, num_epochs_to_process);
amplitude_vector = NaN(1, num_epochs_to_process);
valid_results_count = 0;
initial_phase_anchor_cycles = NaN;
accumulated_residual_phase = 0;
last_epoch_time = NaT;

fprintf('--> 开始GPSense分析流程，历元范围 %d-%d，目标卫星 %s ...\n', epoch_range(1), epoch_range(end), target_satellite_id);

% --- 2. 遍历历元 ---
for loop_idx = 1:num_epochs_to_process
    epoch_idx = epoch_range(loop_idx);
    current_time = obs_data(epoch_idx).time;
    try
        [receiver_pos, receiver_clk_err, sat_states] = ...
            calculate_receiver_position(obs_data, nav_data, epoch_idx);
        
        if ~isfield(sat_states, target_satellite_id), continue; end
        
        sat_obs = obs_data(epoch_idx).data.(target_satellite_id);
        
        if ~(isfield(sat_obs, 'pseudorange') && isfield(sat_obs.pseudorange, 'C1C') && ...
             isfield(sat_obs, 'doppler') && isfield(sat_obs.doppler, 'D1C') && ...
             isfield(sat_obs, 'snr') && isfield(sat_obs.snr, 'S1C'))
            continue;
        end
        
        P_raw = sat_obs.pseudorange.C1C;
        measured_doppler = sat_obs.doppler.D1C;
        cn0_dbhz = sat_obs.snr.S1C;
        
        if isnan(P_raw) || isnan(measured_doppler) || isnan(cn0_dbhz), continue; end
        
        current_amplitude = sqrt(Pn_noise_power * 10^(cn0_dbhz / 10));
        
        sat_state = sat_states.(target_satellite_id);
        sat_pos = sat_state.position;
        sat_vel = sat_state.velocity;
        sat_clk_err = sat_state.clock_error;
        
        geom_range = norm(sat_pos - receiver_pos);
        clock_corrected_pr = P_raw - c * (receiver_clk_err - sat_clk_err);
        processed_pr = clock_corrected_pr - geom_range;
        
        los_vec = (sat_pos - receiver_pos) / norm(sat_pos - receiver_pos);
        fD_geometric = -(sat_vel' * los_vec) / lambda_L1;
        residual_doppler = measured_doppler - fD_geometric;
        
        if isnan(initial_phase_anchor_cycles)
            initial_phase_anchor_cycles = (P_raw - c * (receiver_clk_err - sat_clk_err)) / lambda_L1;
            accumulated_residual_phase = 0;
        end
        
        if ~isnat(last_epoch_time)
            dt = seconds(current_time - last_epoch_time);
            accumulated_residual_phase = accumulated_residual_phase + residual_doppler * dt;
        end
        total_phase_cycles = initial_phase_anchor_cycles + accumulated_residual_phase;
        
        valid_results_count = valid_results_count + 1;
        time_vector(valid_results_count) = current_time;
        phase_vector_rad(valid_results_count) = total_phase_cycles * (2 * pi);
        processed_pr_vector(valid_results_count) = processed_pr;
        amplitude_vector(valid_results_count) = current_amplitude;
        
        last_epoch_time = current_time;
    catch ME
        fprintf('在处理历元 %d 时发生错误: %s\n', epoch_idx, ME.message);
    end
end
fprintf('--> 数据处理循环完成！\n\n');

% --- 3. 数据修剪与绘图 ---
time_vector = time_vector(1:valid_results_count);
phase_vector_rad = phase_vector_rad(1:valid_results_count);
processed_pr_vector = processed_pr_vector(1:valid_results_count);
amplitude_vector = amplitude_vector(1:valid_results_count);
if isempty(time_vector), error('未能计算出任何有效结果。'); end

fprintf('--> 转换时间坐标为北京时间...\n');
time_vector_cst = time_vector - seconds(18) + hours(8);

% 绘制归一化幅度
fprintf('--> 对幅度数据进行归一化处理...\n');
min_amp = min(amplitude_vector);
max_amp = max(amplitude_vector);
if (max_amp - min_amp) > eps 
    normalized_amplitude_vector = (amplitude_vector - min_amp) / (max_amp - min_amp);
else
    normalized_amplitude_vector = ones(size(amplitude_vector)) * 0.5;
end
figure('Name', sprintf('卫星 %s - 归一化幅度', target_satellite_id), 'NumberTitle', 'off');
plot(time_vector_cst, normalized_amplitude_vector, 'g-');
title(sprintf('卫星 %s - 归一化幅度', target_satellite_id));
xlabel('时间 (北京时间)');
ylabel('归一化幅度 (无单位)');
grid on;
datetick('x', 'HH:MM:SS', 'keepticks', 'keeplimits');

% 绘制去趋势后的伪距
window_size = 25; 
pr_trend = movmean(processed_pr_vector, window_size, 'omitnan');
detrended_pr_vector = processed_pr_vector - pr_trend;
figure('Name', sprintf('卫星 %s - 去趋势后的伪距残差', target_satellite_id), 'NumberTitle', 'off');
plot(time_vector_cst, detrended_pr_vector, 'b-');
title(sprintf('卫星 %s - 去趋势后的伪距残差', target_satellite_id));
xlabel('时间 (北京时间)');
ylabel('残余伪距 (米)');
grid on;
datetick('x', 'HH:MM:SS', 'keepticks', 'keeplimits');

% --- 绘制带通滤波后的相位 ---
unwrapped_phase = unwrap(phase_vector_rad);

% --- 计算采样率 ---
if length(time_vector) < 2
    error('时间向量数据点不足，无法计算采样率。');
end

% ====================【核心修正】====================
% 旧的错误代码: mean_dt = mean(diff(seconds(time_vector)));
% 新的正确代码: 先用 diff 计算 duration 数组，再用 seconds 转为数值
mean_dt = mean(seconds(diff(time_vector)));
% ===================================================

Fs = 1 / mean_dt;
fprintf('    计算出的平均采样率为: %.2f Hz\n', Fs);

% --- 设计并应用带通滤波器 ---
f_low = 0.2;
f_high = 5.0;
order = 4;
[b, a] = butter(order, [f_low f_high] / (Fs / 2), 'bandpass');
filtered_phase = filtfilt(b, a, unwrapped_phase);
detrended_phase_rad = filtered_phase;

figure('Name', sprintf('卫星 %s - 带通滤波后的GPSense相位', target_satellite_id), 'NumberTitle', 'off');
plot(time_vector_cst, detrended_phase_rad, 'r-');
title(sprintf('卫星 %s - 带通滤波后的GPSense相位', target_satellite_id));
xlabel('时间 (北京时间)');
ylabel('滤波后相位 (弧度)');
grid on;
datetick('x', 'HH:MM:SS', 'keepticks', 'keeplimits');

fprintf('✅ 图形绘制完成！\n');
end

