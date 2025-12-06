% =========================================================================
% run_gesture_analysis_robust.m
% 功能: 鲁棒手势感知全流程 (最终定稿版)
% 核心逻辑:
%   1. [抗干扰] 天顶角加权 + 密度惩罚 -> 压制手臂噪声。
%   2. [几何层] 加权 RANSAC -> 筛选出几何一致的手势内点 (Inliers)。
%   3. [拟合层] 加权 PCA -> 确定线段的最佳几何姿态 (不改变线段形状)。
%   4. [方向层] 投影相关性 -> 根据时间顺序确定箭头指向 (不改变线段斜率)。
% =========================================================================

%% 1. 环境检查与参数设置
clearvars -except obs_data nav_data; 
if ~exist('obs_data', 'var'), error('请先加载 obs_data!'); end

% --- 路径加载 ---
% 尝试加载常用路径，确保工具函数可用
addpath(genpath('../sky_plot')); 
addpath(genpath('../calculate_clock_bias_and_positon'));
addpath(genpath('../nav_parse'));
addpath(genpath('sky_plot')); 
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

fprintf('--> 启动鲁棒手势感知分析 (RANSAC + PCA + Correlation)...\n');

% --- [Step 1] 分段参数 ---
PARA.smooth_window_sec = 1.5;   
PARA.gvi_threshold     = 6;     
PARA.sampling_rate     = 10;    
PARA.merge_gap_sec     = 0.5;   
PARA.min_duration_sec  = 0.4;   

% --- [Step 2] 鲁棒轨迹参数 ---
TRAJ.gesture_height    = 0.30;  % 手势平面高度 (米)
TRAJ.min_elevation     = 15;    % 最低仰角过滤
TRAJ.energy_th_ratio   = 0.4;   % 能量阈值比例

% --- [Step 2] 抗干扰算法参数 ---
ALG.zenith_safe_deg    = 30;    % 天顶角 < 30度 (头顶) 不受密度惩罚
ALG.az_neighbor_dist   = 20;    % 邻居定义：方位角相差 < 20度
ALG.density_penalty_k  = 1.0;   % 密度惩罚系数
ALG.ransac_iter        = 500;   % RANSAC 迭代次数
ALG.ransac_dist_th     = 0.20;  % RANSAC 内点距离阈值 (米)

%% ================= [Step 1] 数据提取、滤波与分段 =================
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
num_samples = length(t_grid);
num_sats = length(valid_sats);

cn0_matrix = NaN(num_samples, num_sats);
for s_idx = 1:num_sats
    sat_id = valid_sats{s_idx};
    target_snr_code = '';
    for k = 1:min(50, length(obs_data))
        if isfield(obs_data(k).data, sat_id) && isfield(obs_data(k).data.(sat_id), 'snr')
            fields = fieldnames(obs_data(k).data.(sat_id).snr);
            if ismember('S1C', fields), target_snr_code = 'S1C'; elseif ismember('S2I', fields), target_snr_code = 'S2I'; else, if ~isempty(fields), target_snr_code = fields{1}; end; end
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

% SG 预滤波
for s = 1:num_sats
    col = cn0_matrix(:, s); valid = ~isnan(col);
    if sum(valid) > 14
        idx = 1:length(col); filled = interp1(idx(valid), col(valid), idx, 'linear', 'extrap')';
        cn0_matrix(:, s) = sgolayfilt(filled, 2, 7);
    end
end

% 分段逻辑
cn0_smooth = movmean(cn0_matrix, round(PARA.smooth_window_sec * PARA.sampling_rate), 1, 'omitnan');
volatility_matrix = abs(cn0_matrix - cn0_smooth);
gvi_curve_clean = movmean(sum(volatility_matrix, 2, 'omitnan'), 5);
is_active = gvi_curve_clean > PARA.gvi_threshold;
gap_pts = round(PARA.merge_gap_sec * PARA.sampling_rate);
pad_act = [1; is_active; 1];
g_starts = find(diff(pad_act)==-1); g_ends = find(diff(pad_act)==1)-1;
for i=1:length(g_starts), if (g_ends(i)-g_starts(i)+1) < gap_pts && g_starts(i)>1, is_active(g_starts(i):g_ends(i)-1) = 1; end; end
edges = diff([0; is_active; 0]); s_idxs = find(edges==1); e_idxs = find(edges==-1)-1;
min_dur = round(PARA.min_duration_sec * PARA.sampling_rate);
segments = struct('id', {}, 'start_idx', {}, 'end_idx', {}, 'peak_time', {}, 'peak_gvi', {});
cnt = 0;
for i=1:length(s_idxs)
    if (e_idxs(i)-s_idxs(i)) >= min_dur
        cnt = cnt + 1;
        [m_v, m_i] = max(gvi_curve_clean(s_idxs(i):e_idxs(i)));
        segments(cnt).id = cnt; segments(cnt).start_idx = s_idxs(i); segments(cnt).end_idx = e_idxs(i);
        segments(cnt).peak_time = t_grid(s_idxs(i)+m_i-1); segments(cnt).peak_gvi = m_v;
    end
end
fprintf('✅ [Step 1] 分段完成，共 %d 个片段。\n', cnt);

%% ================= [Step 2] 鲁棒 3D 轨迹推演 (核心逻辑) =================
fprintf('--> [Step 2] 开始鲁棒轨迹推演 (加权 RANSAC + PCA拟合 + 时序相关性定方向)...\n');

for i = 1:length(segments)
    seg = segments(i);
    idx_range = seg.start_idx : seg.end_idx;
    seg_times = t_grid(idx_range);
    sub_vol = volatility_matrix(idx_range, :);
    
    % 1. 获取接收机位置
    [~, ep_idx] = min(abs([obs_data.time] - seg.peak_time));
    try
        [rec_pos, ~, sat_states] = calculate_receiver_position(obs_data, nav_data, ep_idx);
    catch
        continue;
    end
    if isempty(rec_pos) || all(isnan(rec_pos)), continue; end
    [lat0, lon0, alt0] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));
    
    % 2. 收集 Hit 卫星及其物理投影点
    hits_data = struct('pos', {}, 'w_final', {}, 't_off', {}, 'id', {}, 'zen_deg', {});
    raw_sats = struct('pos', {}, 'az', {}, 'zen', {}, 'energy', {}, 't_off', {}, 'id', {});
    raw_cnt = 0; hit_cnt = 0;
    
    for s = 1:length(valid_sats)
        sid = valid_sats{s};
        if ~isfield(sat_states, sid), continue; end
        
        sat_p = sat_states.(sid).position;
        [e, n, u] = ecef2enu(sat_p(1)-rec_pos(1), sat_p(2)-rec_pos(2), sat_p(3)-rec_pos(3), lat0, lon0, alt0);
        dist = norm([e, n, u]); 
        vec_u = [e, n, u]/dist;
        
        zen_deg = acosd(vec_u(3)); % 天顶角
        az_deg  = atan2d(e, n); if az_deg<0, az_deg=az_deg+360; end
        
        % 基础过滤
        if vec_u(3) <= 0 || (90-zen_deg) < TRAJ.min_elevation, continue; end
        
        % 投影到 Z = h 平面
        t_int = TRAJ.gesture_height / vec_u(3);
        pt_int = t_int * vec_u;
        
        if norm(pt_int(1:2)) > 3.0, continue; end 
        
        energy = sum(sub_vol(:, s), 'omitnan');
        [~, mx_i] = max(sub_vol(:, s));
        
        raw_cnt = raw_cnt + 1;
        raw_sats(raw_cnt).pos = pt_int;
        raw_sats(raw_cnt).az  = az_deg;
        raw_sats(raw_cnt).zen = zen_deg;
        raw_sats(raw_cnt).energy = energy;
        raw_sats(raw_cnt).t_off = seconds(seg_times(mx_i) - seg_times(1));
        raw_sats(raw_cnt).id = sid;
    end
    
    if raw_cnt < 2, continue; end
    
    % 能量筛选 Hit
    all_e = [raw_sats.energy];
    th_e = max(all_e) * TRAJ.energy_th_ratio;
    hit_candidates = raw_sats(all_e > th_e);
    
    if length(hit_candidates) < 3
        fprintf('   Seg #%d: Hit卫星过少 (%d)，无法进行 RANSAC\n', i, length(hit_candidates));
        continue;
    end
    
    % 3. === [算法核心 A+B] 计算抗干扰权重 ===
    for k = 1:length(hit_candidates)
        cand = hit_candidates(k);
        
        % A. 天顶角加权 (头顶=1, 地平=0)
        w_base = cosd(cand.zen); 
        
        % B. 密度惩罚 (仅针对侧面卫星)
        penalty = 1.0;
        if cand.zen > ALG.zenith_safe_deg 
            n_neighbors = 0;
            for j = 1:length(hit_candidates)
                if k == j, continue; end
                d_az = abs(cand.az - hit_candidates(j).az);
                if d_az > 180, d_az = 360 - d_az; end
                if d_az < ALG.az_neighbor_dist
                    n_neighbors = n_neighbors + 1;
                end
            end
            penalty = 1.0 / (1.0 + ALG.density_penalty_k * n_neighbors);
        end
        
        hit_cnt = hit_cnt + 1;
        % [修复] 显式字段赋值，防止结构体不一致错误
        hits_data(hit_cnt).pos     = cand.pos;
        hits_data(hit_cnt).t_off   = cand.t_off;
        hits_data(hit_cnt).id      = cand.id;
        hits_data(hit_cnt).zen_deg = cand.zen;
        hits_data(hit_cnt).w_final = w_base * penalty;
    end
    
    % 4. === [算法核心 C] 加权 RANSAC 拟合几何直线 ===
    pts_xy = vertcat(hits_data.pos); pts_xy = pts_xy(:, 1:2);
    weights = [hits_data.w_final]';
    
    best_score = -1;
    best_inliers = false(hit_cnt, 1);
    num_pts = size(pts_xy, 1);
    
    for iter = 1:ALG.ransac_iter
        if num_pts >= 2
            sample_idx = randsample(num_pts, 2, true, weights);
        else
            sample_idx = [1, 2];
        end
        p1 = pts_xy(sample_idx(1), :);
        p2 = pts_xy(sample_idx(2), :);
        
        vec = p2 - p1;
        if norm(vec) < 1e-3, continue; end 
        vec = vec / norm(vec);
        normal = [-vec(2), vec(1)]; 
        
        diffs = pts_xy - p1;
        dists = abs(diffs * normal');
        
        is_inlier = dists < ALG.ransac_dist_th;
        current_score = sum(weights(is_inlier));
        
        if current_score > best_score
            best_score = current_score;
            best_inliers = is_inlier;
        end
    end
    
    % 5. === [算法核心] 最终精修 (PCA拟合 + 投影相关性) ===
    % 逻辑: 
    % 1. 用 PCA 确定 Inliers 的最佳拟合直线 (Axis)，确保线段几何最正确。
    % 2. 用 Inliers 的 (位置 vs 时间) 相关性确定方向 (Sign)，确保箭头不指反。
    
    if sum(best_inliers) < 2
        fprintf('   Seg #%d: RANSAC 未找到足够内点\n', i);
        continue;
    end
    
    inlier_pts = pts_xy(best_inliers, :);
    inlier_w   = weights(best_inliers);
    
    % 加权均值
    w_sum = sum(inlier_w);
    mean_w = sum(inlier_pts .* inlier_w) / w_sum;
    
    % 加权 PCA (确定主轴方向，不管正负)
    centered = inlier_pts - mean_w;
    weighted_centered = centered .* sqrt(inlier_w);
    [U_svd, ~, ~] = svd(weighted_centered' * weighted_centered);
    dir_final = U_svd(:, 1)'; % 这是最佳拟合的几何线段方向
    
    % 时序判定 (Correlation)
    % 将所有 Inliers 投影到主轴上，看它们随时间是变大还是变小
    t_inliers = [hits_data(best_inliers).t_off]';
    proj_vals = centered * dir_final';
    
    corr_v = corr(proj_vals, t_inliers);
    
    % 如果相关性为负，说明点随时间推移在轴上是倒退的 -> 翻转方向
    if ~isnan(corr_v) && corr_v < 0
        dir_final = -dir_final; 
    end
    
    traj_az = atan2d(dir_final(1), dir_final(2));
    if traj_az < 0, traj_az = traj_az + 360; end
    
    fprintf('   Seg #%d: 鲁棒方向 %.1f° (Inliers权重占比: %.2f/%.2f) Corr: %.2f\n', ...
        i, traj_az, best_score, sum(weights), corr_v);

    %% 6. 可视化 (2D 投影平面 + 权重展示)
    fig_name = sprintf('Robust Analysis #%d (Az=%.1f)', i, traj_az);
    f = figure('Name', fig_name, 'Position', [100+(i-1)*20, 200, 600, 600], 'Color', 'w');
    ax = axes('Parent', f); hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
    xlabel('East (m)'); ylabel('North (m)');
    title({sprintf('手势 #%d 鲁棒分析结果', i), '大小=权重 | 红色=Inliers(手) | 灰色=Outliers(干扰)'});
    
    % 画接收机
    plot(ax, 0, 0, '^', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    
    % 画所有 Hit 点 (大小代表权重，颜色代表是否 Inlier)
    for k = 1:length(hits_data)
        pt = hits_data(k).pos;
        w  = hits_data(k).w_final;
        is_in = best_inliers(k);
        
        if is_in
            col = 'r'; % 手势点
        else
            col = [0.7 0.7 0.7]; % 干扰点 (灰色)
        end
        
        % Marker 大小与权重成正比 (min 5, max 15)
        ms = 5 + w * 15; 
        plot(ax, pt(1), pt(2), 'o', 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'MarkerSize', ms);
    end
    
    % 画最终轨迹 (基于 Inliers 的 PCA 轴)
    line_len = 2.0;
    p_start = mean_w - line_len * dir_final;
    p_end   = mean_w + line_len * dir_final;
    plot(ax, [p_start(1), p_end(1)], [p_start(2), p_end(2)], 'b--', 'LineWidth', 1.5);
    
    % 画向量箭头 (表示相关性确定的方向)
    quiver(ax, mean_w(1)-0.5*dir_final(1), mean_w(2)-0.5*dir_final(2), ...
           dir_final(1), dir_final(2), 'Color', 'b', 'LineWidth', 3, 'MaxHeadSize', 0.5, 'AutoScale', 'off');
       
    xlim([-2.5 2.5]); ylim([-2.5 2.5]);
end
fprintf('✅ 所有分析完成。\n');