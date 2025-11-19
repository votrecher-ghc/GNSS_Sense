% =========================================================================
% step1_final_visual_v4.m
% 功能: 手势检测 Step 1 (连通域合并 + 图3增加阈值线 + 数据导出)
% =========================================================================

%% 1. 准备工作与参数设置
clearvars -except obs_data nav_data; 
if ~exist('obs_data', 'var'), error('请先加载 obs_data!'); end

% --- 核心参数 ---
PARA.smooth_window_sec = 1.5;  % 基线平滑窗口(秒)
PARA.gvi_threshold     = 10;    % GVI 阈值 (dB)
PARA.sampling_rate     = 10;   % 初始采样率(会自动校准)

% 连通域参数 (解决 M 型波断裂)
PARA.merge_gap_sec     = 0.5;  % 中间断开小于 0.5 秒视为不断开
PARA.min_duration_sec  = 0.4;  % 小于 0.4 秒视为噪声

fprintf('--> [Step 1] 开始手势检测 (连通域合并模式)...\n');

%% 2. 数据提取与对齐 (保持不变)
all_sat_ids = {};
scan_range = unique([1:min(100, length(obs_data)), max(1, length(obs_data)-100):length(obs_data)]);
for i = scan_range
    if ~isempty(obs_data(i).data), all_sat_ids = [all_sat_ids, fieldnames(obs_data(i).data)']; end
end
unique_sat_ids = unique(all_sat_ids);

valid_sats = {};
for i = 1:length(unique_sat_ids)
    sid = unique_sat_ids{i};
    if ismember(sid(1), ['G','C','E','J']), valid_sats{end+1} = sid; end
end

raw_times = [obs_data.time];
mean_dt = mean(diff(raw_times));
if mean_dt < seconds(0.05), PARA.sampling_rate = 20;
elseif mean_dt < seconds(0.2), PARA.sampling_rate = 10;
else, PARA.sampling_rate = 1; end

t_grid = (min(raw_times) : seconds(1/PARA.sampling_rate) : max(raw_times))'; 
num_samples = length(t_grid);
num_sats = length(valid_sats);

cn0_matrix = NaN(num_samples, num_sats);
for s_idx = 1:num_sats
    sat_id = valid_sats{s_idx};
    target_snr_code = '';
    for k = 1:min(50, length(obs_data)) 
        if isfield(obs_data(k).data, sat_id) && isfield(obs_data(k).data.(sat_id), 'snr')
            fields = fieldnames(obs_data(k).data.(sat_id).snr);
            if ~isempty(fields)
                if ismember('S1C', fields), target_snr_code = 'S1C'; 
                elseif ismember('S2I', fields), target_snr_code = 'S2I';
                else, target_snr_code = fields{1}; 
                end
                break;
            end
        end
    end
    if isempty(target_snr_code), continue; end
    
    s_times = []; s_cn0 = [];
    for k = 1:length(obs_data)
        if isfield(obs_data(k).data, sat_id) && isfield(obs_data(k).data.(sat_id).snr, target_snr_code)
            val = obs_data(k).data.(sat_id).snr.(target_snr_code);
            if ~isnan(val) && val > 10
                s_times = [s_times; obs_data(k).time]; %#ok<AGROW>
                s_cn0   = [s_cn0; val]; %#ok<AGROW>
            end
        end
    end
    if length(s_times) > 20
        [u_times, u_idx] = unique(s_times);
        cn0_matrix(:, s_idx) = interp1(u_times, s_cn0(u_idx), t_grid, 'linear', NaN);
    end
end

%% 3. 计算波动与分段 (连通域逻辑)
cn0_smooth = movmean(cn0_matrix, round(PARA.smooth_window_sec * PARA.sampling_rate), 1, 'omitnan');
volatility_matrix = abs(cn0_matrix - cn0_smooth);
gvi_curve = sum(volatility_matrix, 2, 'omitnan');
gvi_curve_clean = movmean(gvi_curve, 5);

fprintf('--> 正在执行连通域分段...\n');

is_active = gvi_curve_clean > PARA.gvi_threshold;

% 填补缝隙
min_gap_idx = round(PARA.merge_gap_sec * PARA.sampling_rate);
padded_active = [1; is_active; 1];
gap_starts = find(diff(padded_active) == -1); 
gap_ends   = find(diff(padded_active) == 1) - 1; 

for i = 1:length(gap_starts)
    len = gap_ends(i) - gap_starts(i) + 1;
    if len < min_gap_idx && gap_starts(i) > 1 && gap_ends(i) < length(padded_active)
        idx_s = gap_starts(i); 
        idx_e = gap_ends(i) - 1; 
        if idx_s <= idx_e, is_active(idx_s:idx_e) = 1; end
    end
end

edges = diff([0; is_active; 0]);
s_indices = find(edges == 1);
e_indices = find(edges == -1) - 1;

segments = struct('start_idx', {}, 'end_idx', {}, 'peak_time', {}, 'peak_gvi', {});
num_segs = 0;
min_dur_idx = round(PARA.min_duration_sec * PARA.sampling_rate);

for i = 1:length(s_indices)
    s = s_indices(i);
    e = e_indices(i);
    if (e - s) >= min_dur_idx
        num_segs = num_segs + 1;
        [max_val, max_loc] = max(gvi_curve_clean(s:e));
        peak_idx = s + max_loc - 1;
        segments(num_segs).id = num_segs;
        segments(num_segs).start_idx = s;
        segments(num_segs).end_idx = e;
        segments(num_segs).peak_idx = peak_idx;
        segments(num_segs).peak_time = t_grid(peak_idx);
        segments(num_segs).peak_gvi = max_val;
    end
end

fprintf('✅ 识别到 %d 个手势片段。\n', num_segs);

%% 4. 结果可视化 (图3增加阈值虚线)
figure('Name', 'Gesture Detection Analysis', 'Position', [50, 50, 1000, 800]);

% --- 图1 ---
subplot(3, 1, 1);
plot(t_grid, cn0_matrix);
title(sprintf('1. 全星座 C/N0 原始数据 (%d 颗卫星)', num_sats));
ylabel('SNR'); axis tight; grid on;
yl = ylim;
for i=1:num_segs
    patch([t_grid(segments(i).start_idx) t_grid(segments(i).end_idx) t_grid(segments(i).end_idx) t_grid(segments(i).start_idx)], ...
          [yl(1) yl(1) yl(2) yl(2)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end

% --- 图2 ---
subplot(3, 1, 2);
plot(t_grid, gvi_curve_clean, 'k-', 'LineWidth', 1); hold on;
yline(PARA.gvi_threshold, 'b--', '阈值'); % 图2也有阈值线
title(sprintf('2. 波动指数,阈值:%d',PARA.gvi_threshold));
ylabel('GVI'); axis tight; grid on;
yl2 = ylim;
for i=1:num_segs
    x_s = t_grid(segments(i).start_idx); x_e = t_grid(segments(i).end_idx);
    patch([x_s x_e x_e x_s], [yl2(1) yl2(1) yl2(2) yl2(2)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end

% --- 图3：检测到的手势动作片段 (详细高亮 + 标签 + 红色背景 + 阈值线) ---
subplot(3, 1, 3);
plot(t_grid, gvi_curve_clean, 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5); hold on;

% [新增] 在图3中添加蓝色的阈值虚线
yline(PARA.gvi_threshold, 'b--', '阈值'); 

title('3. 检测到的手势动作片段详情 (红色高亮)');
ylabel('GVI Detail'); xlabel('时间');
grid on; axis tight;

yl3 = ylim; 
for i = 1:num_segs
    idx_range = segments(i).start_idx : segments(i).end_idx;
    t_s = t_grid(segments(i).start_idx);
    t_e = t_grid(segments(i).end_idx);
    
    % 1. 红色背景块
    patch([t_s t_e t_e t_s], [yl3(1) yl3(1) yl3(2) yl3(2)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
      
    % 2. 红色波形线
    plot(t_grid(idx_range), gvi_curve_clean(idx_range), 'r-', 'LineWidth', 2);
    
    % 3. 标签
    text(segments(i).peak_time, segments(i).peak_gvi, sprintf('  #%d', i), ...
        'Color', 'r', 'FontWeight', 'bold', 'FontSize', 11, ...
        'VerticalAlignment', 'bottom');
end

%% 5. 计算并导出起止时间
fprintf('\n--> 正在生成结果表格...\n');

ids = [segments.id]';
start_idxs = [segments.start_idx]';
end_idxs   = [segments.end_idx]';
dur_samples = end_idxs - start_idxs + 1;

% 假设原始 obs_data 时间是 UTC，转为北京时间 (+8h)
start_times_bjt = t_grid(start_idxs) + hours(8);
end_times_bjt   = t_grid(end_idxs)   + hours(8);
dur_seconds = seconds(end_times_bjt - start_times_bjt);

% 表格 1: 采样点信息
T_Index = table(ids, start_idxs, end_idxs, dur_samples, ...
    'VariableNames', {'GestureID', 'Start_Index', 'End_Index', 'Duration_Points'});

% 表格 2: 北京时间信息
T_Time = table(ids, start_times_bjt, end_times_bjt, dur_seconds, ...
    'VariableNames', {'GestureID', 'Start_Time_BJT', 'End_Time_BJT', 'Duration_Sec'});

fprintf('\n=== 表格 1: 采样点信息 ===\n');
disp(T_Index);

fprintf('\n=== 表格 2: 北京时间信息 (UTC+8) ===\n');
disp(T_Time);

% % 保存到 CSV 文件
% writetable(T_Index, 'Gesture_Result_Indices.csv');
% writetable(T_Time,  'Gesture_Result_Times_BJT.csv');
% 
% fprintf('✅ 结果已保存为: \n   - Gesture_Result_Indices.csv\n   - Gesture_Result_Times_BJT.csv\n');





% 
% 
% 
% %% 1. 准备工作
% clearvars -except obs_data nav_data; % 保留已加载的数据，避免重复解析
% if ~exist('obs_data', 'var')
%     error('请先加载 obs_data! (运行 parse_rinex_obs.m)');
% end
% 
% % 参数设置 (根据您的实验环境微调)
% PARA.smooth_window_sec = 1.5;  % 滑动平均窗口大小 (秒)，用于计算基线
% PARA.gvi_threshold     = 10;   % GVI触发阈值 (dB)，低于此值的波动被视为噪声
% PARA.min_peak_dist     = 1.0;  % 最小峰值间隔 (秒)，防止同一个手势被切成两段
% PARA.sampling_rate     = 10;   % 假设采样率为10Hz (或者根据数据自动计算)
% 
% fprintf('--> 正在提取全星座 C/N0 数据...\n');
% 
% %% 2. 提取并对齐所有卫星的 C/N0 数据
% % 2.1 获取所有出现的卫星ID
% all_sat_ids = {};
% for i = 1:length(obs_data)
%     if ~isempty(obs_data(i).data)
%         all_sat_ids = [all_sat_ids, fieldnames(obs_data(i).data)'];
%     end
% end
% unique_sat_ids = unique(all_sat_ids);
% % 过滤只支持的系统 (G/C/E/J) - 根据您的analyze脚本逻辑
% valid_sats = {};
% for i = 1:length(unique_sat_ids)
%     sid = unique_sat_ids{i};
%     if ismember(sid(1), ['G','C','E','J']) 
%         valid_sats{end+1} = sid; 
%     end
% end
% fprintf('    共识别到 %d 颗有效卫星。\n', length(valid_sats));
% 
% % 2.2 构建统一时间轴
% raw_times = [obs_data.time];
% % 转换为北京时间用于绘图
% t_grid = (min(raw_times) : seconds(1/PARA.sampling_rate) : max(raw_times))'; 
% num_samples = length(t_grid);
% num_sats = length(valid_sats);
% 
% % 2.3 构建 C/N0 矩阵 [时间 x 卫星]
% cn0_matrix = NaN(num_samples, num_sats);
% 
% for s_idx = 1:num_sats
%     sat_id = valid_sats{s_idx};
%     % 确定该卫星使用的 SNR 代码 (参考您的 analyze 脚本)
%     sys = sat_id(1);
%     if sys == 'C', snr_code = 'S2I'; else, snr_code = 'S1C'; end 
%     
%     % 提取该卫星的时间和SNR
%     s_times = [];
%     s_cn0   = [];
%     for k = 1:length(obs_data)
%         if isfield(obs_data(k).data, sat_id) && ...
%            isfield(obs_data(k).data.(sat_id), 'snr') && ...
%            isfield(obs_data(k).data.(sat_id).snr, snr_code)
%             
%             val = obs_data(k).data.(sat_id).snr.(snr_code);
%             if ~isnan(val) && val > 0
%                 s_times = [s_times; obs_data(k).time]; %#ok<AGROW>
%                 s_cn0   = [s_cn0; val]; %#ok<AGROW>
%             end
%         end
%     end
%     
%     % 对齐到统一时间轴 (线性插值)
%     if length(s_times) > 10
%         % 去重
%         [u_times, u_idx] = unique(s_times);
%         u_cn0 = s_cn0(u_idx);
%         % 插值填充
%         cn0_matrix(:, s_idx) = interp1(u_times, u_cn0, t_grid, 'linear', NaN);
%     end
% end
% % 将NaN (无信号) 替换为0，避免计算出错，或者保持NaN不参与计算
% % 这里我们保持NaN，但在计算GVI时忽略NaN
% 
% %% 3. 计算全局波动指数 (GVI)
% fprintf('--> 正在计算 GVI...\n');
% 
% % 3.1 计算每颗星的波动量 V_k(t)
% % 方法：原始值 - 平滑基线
% window_points = round(PARA.smooth_window_sec * PARA.sampling_rate);
% cn0_smooth = movmean(cn0_matrix, window_points, 1, 'omitnan');
% volatility_matrix = abs(cn0_matrix - cn0_smooth);
% 
% % 3.2 求和得到 GVI (忽略NaN)
% gvi_curve = sum(volatility_matrix, 2, 'omitnan');
% 
% % 3.3 对 GVI 自身再做一次轻微平滑，去除尖刺噪声
% gvi_curve_clean = movmean(gvi_curve, 5); 
% 
% %% 4. 核心：切割波动段 (Segmentation)
% fprintf('--> 正在执行分段识别...\n');
% 
% % 使用 findpeaks 寻找显著的波动事件
% % MinPeakHeight: GVI必须超过此值才算事件
% % MinPeakDistance: 两个手势之间的最小间隔
% [pks, locs, w, p] = findpeaks(gvi_curve_clean, ...
%     'MinPeakHeight', PARA.gvi_threshold, ...
%     'MinPeakDistance', PARA.min_peak_dist * PARA.sampling_rate, ...
%     'WidthReference', 'halfheight');
% 
% % 定义段的开始和结束
% % 简单策略：以峰值为中心，左右各扩展一定的宽度，或者寻找下降到阈值的时间点
% segments = struct('start_idx', {}, 'end_idx', {}, 'peak_time', {}, 'peak_gvi', {});
% 
% for i = 1:length(locs)
%     peak_idx = locs(i);
%     
%     % 向左搜索开始点 (GVI 降到 阈值的 20% 以下)
%     start_idx = peak_idx;
%     while start_idx > 1 && gvi_curve_clean(start_idx) > (0.2 * pks(i))
%         start_idx = start_idx - 1;
%     end
%     
%     % 向右搜索结束点
%     end_idx = peak_idx;
%     while end_idx < num_samples && gvi_curve_clean(end_idx) > (0.2 * pks(i))
%         end_idx = end_idx + 1;
%     end
%     
%     segments(i).start_idx = start_idx;
%     segments(i).end_idx = end_idx;
%     segments(i).peak_time = t_grid(peak_idx);
%     segments(i).peak_gvi = pks(i);
% end
% 
% fprintf('✅ 成功识别到 %d 个手势波动段！\n', length(segments));
% 
% %% 5. 可视化结果
% figure('Name', 'GVI Segmentation', 'Position', [100, 100, 1000, 600]);
% 
% % 上图：所有卫星的 C/N0 叠画
% subplot(2,1,1);
% plot(t_grid + hours(8), cn0_matrix); % 转北京时间
% title(sprintf('全星座 C/N0 原始数据 (%d 颗卫星)', num_sats));
% ylabel('SNR (dB-Hz)');
% grid on;
% axis tight;
% datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');
% 
% % 下图：GVI 曲线与识别结果
% subplot(2,1,2);
% plot(t_grid + hours(8), gvi_curve_clean, 'k', 'LineWidth', 1.5); hold on;
% ylabel('全局波动指数 (GVI)');
% xlabel('时间 (北京时间)');
% title(sprintf('GVI 指数与分段结果 (阈值: %.1f)', PARA.gvi_threshold));
% 
% % 绘制识别出的段
% y_limits = ylim;
% for i = 1:length(segments)
%     t_start = t_grid(segments(i).start_idx) + hours(8);
%     t_end   = t_grid(segments(i).end_idx) + hours(8);
%     
%     % 绘制高亮区域
%     patch([t_start t_end t_end t_start], ...
%           [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
%           'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
%       
%     % 标记峰值
%     text(t_grid(segments(i).peak_time == t_grid) + hours(8), ...
%          segments(i).peak_gvi, ...
%          sprintf('Seg #%d', i), 'Vert', 'bottom', 'Horiz', 'center');
% end
% grid on;
% datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');
% 
% %% 6. (可选) 输出第一段的数据供下一步分析
% if ~isempty(segments)
%     fprintf('\n准备进入难点2分析...\n');
%     fprintf('Seg #1 范围: %s 到 %s\n', ...
%         datestr(t_grid(segments(1).start_idx)+hours(8), 'HH:MM:SS.fff'), ...
%         datestr(t_grid(segments(1).end_idx)+hours(8), 'HH:MM:SS.fff'));
%     
%     % 将波动量矩阵传递给下一步
%     % segment_volatility = volatility_matrix(segments(1).start_idx:segments(1).end_idx, :);
% end
% 
% 
% 
% 
% 
