% =========================================================================
% run_gesture_analysis_3d_separate_v3.m
% 功能: 手势感知全流程 (独立窗口 + 视觉通透版)
% 改进:
%   1. [空间感] 射线不再终结于手势平面，而是延伸至高空(Z=2m)，模拟“通天”视线。
%   2. [精致] 手势箭头变细 (LineWidth=2)。
%   3. [完整] Miss 卫星也绘制虚线射线，展示完整几何环境。
% =========================================================================

%% 1. 环境检查与参数设置
clearvars -except obs_data nav_data; 
if ~exist('obs_data', 'var'), error('请先加载 obs_data!'); end

% --- [Step 1] 分段参数 ---
PARA.smooth_window_sec = 1.5;   
PARA.gvi_threshold     = 6;     
PARA.sampling_rate     = 10;    
PARA.merge_gap_sec     = 0.5;   
PARA.min_duration_sec  = 0.4;   

% --- [Step 2] 3D 轨迹参数 ---
TRAJ.gesture_height    = 0.30;  % 手势平面高度 (0.3m)
TRAJ.sky_height        = 1.0;   % [新增] 射线绘制的顶部高度 (让视觉更开阔)
TRAJ.energy_threshold  = 0.4;   
TRAJ.min_hit_sats      = 2;     
TRAJ.miss_conflict_dist= 0.15;  
TRAJ.min_elevation     = 15;    

addpath(genpath('sky_plot')); 
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

fprintf('--> 启动手势感知全流程分析 (视觉通透版 v3)...\n');

%% ================= [Step 1] 数据提取、滤波与分段 =================
% (核心逻辑保持不变，确保数据处理的一致性)

fprintf('--> [Step 1] 提取全星座数据...\n');
all_sat_ids = {};
for i = 1:min(100, length(obs_data)), if ~isempty(obs_data(i).data), all_sat_ids = [all_sat_ids, fieldnames(obs_data(i).data)']; end; end
unique_sat_ids = unique(all_sat_ids);
valid_sats = {};
for i = 1:length(unique_sat_ids)
    sid = unique_sat_ids{i};
    if ismember(sid(1), ['G','C','E','J']), valid_sats{end+1} = sid; end
end

raw_times = [obs_data.time];
mean_dt = mean(diff(raw_times));
if mean_dt < seconds(0.05), PARA.sampling_rate = 20; elseif mean_dt < seconds(0.2), PARA.sampling_rate = 10; else, PARA.sampling_rate = 1; end

t_grid = (min(raw_times) : seconds(1/PARA.sampling_rate) : max(raw_times))'; 
t_grid_plot = t_grid + hours(8) - seconds(20); 
num_samples = length(t_grid);
num_sats = length(valid_sats);

cn0_matrix = NaN(num_samples, num_sats);
for s_idx = 1:num_sats
    sat_id = valid_sats{s_idx};
    target_snr_code = '';
    sys = sat_id(1);
    if sys=='C', codes={'S2I','S2X','S1I','S6I','S7I'}; else, codes={'S1C','S1X','S2C','S2X','S5X'}; end
    for k = 1:min(50, length(obs_data))
        if isfield(obs_data(k).data, sat_id) && isfield(obs_data(k).data.(sat_id), 'snr')
            avail_codes = fieldnames(obs_data(k).data.(sat_id).snr);
            for c = 1:length(codes), if ismember(codes{c}, avail_codes), target_snr_code = codes{c}; break; end; end
            if ~isempty(target_snr_code), break; end
        end
    end
    if isempty(target_snr_code), continue; end
    s_times = []; s_cn0 = [];
    for k = 1:length(obs_data)
        if isfield(obs_data(k).data, sat_id) && isfield(obs_data(k).data.(sat_id).snr, target_snr_code)
            val = obs_data(k).data.(sat_id).snr.(target_snr_code);
            if ~isnan(val) && val > 10, s_times = [s_times; obs_data(k).time]; s_cn0 = [s_cn0; val]; end
        end
    end
    if length(s_times) > 20
        [u_times, u_idx] = unique(s_times);
        cn0_matrix(:, s_idx) = interp1(u_times, s_cn0(u_idx), t_grid, 'linear', NaN);
    end
end

% 1.2 预滤波
fprintf('--> [Step 1] 执行 Savitzky-Golay 预滤波去噪...\n');
sg_order = 2; sg_len = 7;
for s = 1:num_sats
    col = cn0_matrix(:, s);
    valid = ~isnan(col);
    if sum(valid) > sg_len*2
        idx = 1:length(col);
        filled = interp1(idx(valid), col(valid), idx, 'linear', 'extrap')';
        cn0_matrix(:, s) = sgolayfilt(filled, sg_order, sg_len);
    end
end

% 1.3 分段
cn0_smooth = movmean(cn0_matrix, round(PARA.smooth_window_sec * PARA.sampling_rate), 1, 'omitnan');
volatility_matrix = abs(cn0_matrix - cn0_smooth);
gvi_curve = sum(volatility_matrix, 2, 'omitnan');
gvi_curve_clean = movmean(gvi_curve, 5);

is_active = gvi_curve_clean > PARA.gvi_threshold;
gap_pts = round(PARA.merge_gap_sec * PARA.sampling_rate);
pad_act = [1; is_active; 1];
g_starts = find(diff(pad_act)==-1); g_ends = find(diff(pad_act)==1)-1;
for i=1:length(g_starts)
    if (g_ends(i)-g_starts(i)+1) < gap_pts && g_starts(i)>1
        is_active(g_starts(i):g_ends(i)-1) = 1;
    end
end

edges = diff([0; is_active; 0]);
s_idxs = find(edges==1); e_idxs = find(edges==-1)-1;
min_dur = round(PARA.min_duration_sec * PARA.sampling_rate);

segments = struct('id', {}, 'start_idx', {}, 'end_idx', {}, 'peak_time', {}, 'peak_gvi', {});
cnt = 0;
for i=1:length(s_idxs)
    if (e_idxs(i)-s_idxs(i)) >= min_dur
        cnt = cnt + 1;
        [m_v, m_i] = max(gvi_curve_clean(s_idxs(i):e_idxs(i)));
        peak_idx = s_idxs(i) + m_i - 1;
        segments(cnt).id = cnt;
        segments(cnt).start_idx = s_idxs(i);
        segments(cnt).end_idx = e_idxs(i);
        segments(cnt).peak_idx = peak_idx;
        segments(cnt).peak_time = t_grid(peak_idx);
        segments(cnt).peak_gvi = m_v;
    end
end
fprintf('✅ [Step 1] 完成。识别到 %d 个有效手势片段。\n', cnt);

% 1.4 可视化分段
figure('Name', 'Step 1: Segmentation Summary', 'Position', [50, 500, 800, 300]);
plot(t_grid_plot, gvi_curve_clean, 'k', 'LineWidth', 1); hold on;
yline(PARA.gvi_threshold, 'b--', 'Threshold');
yl = ylim;
for i=1:length(segments)
    idx = segments(i).start_idx : segments(i).end_idx;
    plot(t_grid_plot(idx), gvi_curve_clean(idx), 'r', 'LineWidth', 2);
end
title(sprintf('手势分段概览 (阈值=%d)', PARA.gvi_threshold));
xlabel('Time (BJT)'); ylabel('GVI');
datetick('x','HH:MM:ss','keepticks','keeplimits'); grid on;


%% ================= [Step 2] 3D 轨迹推演 (独立窗口 - 视觉通透版) =================

fprintf('--> [Step 2] 开始 3D 空间轨迹推演 (通透射线版)...\n');

traj_colors = lines(length(segments));

for i = 1:length(segments)
    seg = segments(i);
    idx_range = seg.start_idx : seg.end_idx;
    seg_times = t_grid(idx_range);
    
    sub_vol = volatility_matrix(idx_range, :);
    
    [~, ep_idx] = min(abs([obs_data.time] - seg.peak_time));
    [rec_pos, ~, sat_states] = calculate_receiver_position(obs_data, nav_data, ep_idx);
    if isempty(rec_pos), continue; end
    [lat0, lon0, alt0] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));
    
    % 存储 3D 点及其方向向量 (vec_u)
    sat_pts = struct('pos', {}, 'vec_u', {}, 'energy', {}, 't_off', {});
    pt_cnt = 0;
    
    for s = 1:length(valid_sats)
        sid = valid_sats{s};
        if ~isfield(sat_states, sid), continue; end
        
        sat_p = sat_states.(sid).position;
        [e, n, u] = ecef2enu(sat_p(1)-rec_pos(1), sat_p(2)-rec_pos(2), sat_p(3)-rec_pos(3), lat0, lon0, alt0);
        vec = [e, n, u]; dist = norm(vec); vec_u = vec/dist;
        
        if asind(vec_u(3)) < TRAJ.min_elevation, continue; end
        
        if vec_u(3) > 0
            t_int = TRAJ.gesture_height / vec_u(3);
            pt_int = t_int * vec_u;
            
            % 稍微放宽范围，避免边缘卫星被切掉
            if norm(pt_int(1:2)) < 3.0 
                pt_cnt = pt_cnt + 1;
                sat_pts(pt_cnt).pos = pt_int;
                sat_pts(pt_cnt).vec_u = vec_u; % 关键：保存单位向量用于画射线
                sat_pts(pt_cnt).energy = sum(sub_vol(:, s), 'omitnan');
                [~, mx_i] = max(sub_vol(:, s));
                sat_pts(pt_cnt).t_off = seconds(seg_times(mx_i) - seg_times(1));
            end
        end
    end
    
    if pt_cnt < 2, continue; end
    
    eners = [sat_pts.energy];
    th_e = max(eners) * TRAJ.energy_threshold;
    hits = sat_pts(eners > th_e);
    misses = sat_pts(eners <= th_e);
    
    if length(hits) < TRAJ.min_hit_sats
        fprintf('   Seg #%d: 有效卫星不足 (%d)，跳过。\n', i, length(hits));
        continue;
    end
    
    % PCA & 时序
    coords = vertcat(hits.pos);
    pts_xy = coords(:, 1:2);
    mean_xy = mean(pts_xy);
    centered = pts_xy - mean_xy;
    [coeff, ~, ~] = pca(centered);
    dir_xy = coeff(:, 1)';
    
    proj = centered * dir_xy';
    times = [hits.t_off]';
    corr_v = corr(proj, times);
    if ~isnan(corr_v) && corr_v < 0, dir_xy = -dir_xy; end
    
    proj_vals = (pts_xy - mean_xy) * dir_xy';
    p_start_2d = mean_xy + (min(proj_vals)-0.1) * dir_xy;
    p_end_2d   = mean_xy + (max(proj_vals)+0.1) * dir_xy;
    p_start = [p_start_2d, TRAJ.gesture_height];
    p_end   = [p_end_2d,   TRAJ.gesture_height];
    
    conflict = false;
    vec_seg = p_end - p_start; len_sq = dot(vec_seg, vec_seg);
    for m=1:length(misses)
        mp = misses(m).pos;
        if len_sq>0
            t = dot(mp-p_start, vec_seg)/len_sq;
            if t>0 && t<1 && norm(mp-(p_start+t*vec_seg)) < TRAJ.miss_conflict_dist
                conflict = true; break;
            end
        end
    end
    
    % ================= [绘图：视觉通透版] =================
    fig_name = sprintf('Gesture #%d 3D View (T=%s)', i, datestr(seg.peak_time + hours(8)-seconds(20), 'HH:MM:SS'));
    f = figure('Name', fig_name, 'Position', [100 + (i-1)*30, 100 + (i-1)*30, 900, 700], 'Color', 'w');
    ax3d = axes('Parent', f);
    hold(ax3d, 'on'); grid(ax3d, 'on'); axis(ax3d, 'equal'); view(ax3d, 3);
    xlabel(ax3d, 'East (m)'); ylabel(ax3d, 'North (m)'); zlabel(ax3d, 'Up (m)');
    
    az = atan2d(dir_xy(1), dir_xy(2)); if az<0, az=az+360; end
    str_conflict = ""; if conflict, str_conflict = "[冲突]"; end
    title(ax3d, {sprintf('手势 #%d: 方向 %.1f° (相关性 %.2f) %s', i, az, abs(corr_v), str_conflict), ...
                 '彩色实线: Hit卫星 | 灰色虚线: Miss卫星'});

    % 1. 设置更开阔的坐标系
    view_range = 4.0; % [修改] 扩大到 4米
    xlim(ax3d, [-view_range, view_range]);
    ylim(ax3d, [-view_range, view_range]);
    zlim(ax3d, [0, TRAJ.sky_height + 0.5]); % [修改] 显示到2.5米高

    % 2. 画环境 (手势平面)
    plot3(ax3d, 0,0,0, '^', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'DisplayName', 'Receiver');
    plane_size = view_range; 
    patch(ax3d, [-plane_size plane_size plane_size -plane_size], ...
                [-plane_size -plane_size plane_size plane_size], ...
          [TRAJ.gesture_height TRAJ.gesture_height TRAJ.gesture_height TRAJ.gesture_height], ...
          [0.95 0.95 1], 'FaceAlpha', 0.1, 'EdgeColor', [0.8 0.8 1], 'LineStyle', ':');
    
    col = traj_colors(i, :);
    
    % 3. [改进] 画 Miss 卫星的“通天”射线
    if ~isempty(misses)
        for m = 1:length(misses)
            mp = misses(m).pos;
            vec = misses(m).vec_u;
            
            % 计算高空点 (Z = sky_height)
            if vec(3) > 0
                t_top = TRAJ.sky_height / vec(3);
                p_top = t_top * vec;
                
                % 画射线: 原点 -> 高空点 (虚线)
                plot3(ax3d, [0 p_top(1)], [0 p_top(2)], [0 p_top(3)], '--', ...
                      'Color', [0.85 0.85 0.85], 'LineWidth', 0.5); 
                
                % 画交点 (手势平面上的点)
                plot3(ax3d, mp(1), mp(2), mp(3), 'o', ...
                      'Color', [0.7 0.7 0.7], 'MarkerSize', 3);
            end
        end
    end

    % 4. [改进] 画 Hit 卫星的“通天”射线
    for k=1:length(hits)
        hp = hits(k).pos;
        vec = hits(k).vec_u;
        
        if vec(3) > 0
            t_top = TRAJ.sky_height / vec(3);
            p_top = t_top * vec;
            
            % 画射线: 原点 -> 高空点 (实线, 半透明)
            plot3(ax3d, [0 p_top(1)], [0 p_top(2)], [0 p_top(3)], '-', ...
                  'Color', [col, 0.5], 'LineWidth', 1.5);
            
            % 画交点 (实心大点)
            plot3(ax3d, hp(1), hp(2), hp(3), 'o', ...
                  'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
        end
    end
    
    % 5. 画轨迹箭头 (变细)
    quiver3(ax3d, p_start(1), p_start(2), p_start(3), ...
            p_end(1)-p_start(1), p_end(2)-p_start(2), 0, ...
            0, 'Color', col, 'LineWidth', 2, 'MaxHeadSize', 0.5); % [修改] LineWidth=2
    
    fprintf('   Seg #%d: 方向 %.1f 度 (相关性 %.2f) %s\n', i, az, abs(corr_v), str_conflict);
end

fprintf('✅ 全流程分析完成！\n');