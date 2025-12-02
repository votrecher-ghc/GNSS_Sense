% =========================================================================
% step2_direction_estimation_v9_3D_Rigorous.m
% 功能: 基于三维视距切割模型 (3D LoS Cutting Model) 的手势轨迹推断
% 改进点:
%   1. 【3D建模】建立局部 ENU 坐标系，真实还原 "接收机-手-卫星" 的空间几何。
%   2. 【物理平面】假设手在高度 h 处运动，计算视线与该平面的 3D 交点。
%   3. 【视觉升级】输出带有空间感的 3D 轨迹图，解决“模型粗糙”的问题。
% =========================================================================

%% 1. 检查环境
clearvars -except obs_data nav_data segments volatility_matrix valid_sats t_grid PARA;
if ~exist('segments', 'var') || isempty(segments)
    error('错误: 未找到 segments 变量。请先运行 step1_segmentation_GVI.m');
end

addpath(genpath('sky_plot')); 
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

%% 2. 初始化 3D 画布
fprintf('--> 正在初始化 3D 空间轨迹分析...\n');
fig_handle = figure('Name', '3D Gesture Trajectory Model', 'Position', [100, 100, 1200, 900], 'Color', 'w');
ax = axes('Parent', fig_handle);
hold(ax, 'on');
grid(ax, 'on');
axis(ax, 'equal');
view(ax, 3); % 开启3D视角

% 设置坐标轴标签
xlabel(ax, 'East (m)');
ylabel(ax, 'North (m)');
zlabel(ax, 'Up (m)');
title(ax, '基于三维视距路径切割的手势轨迹推演');

% --- 参数设置 ---
TRAJ_PARA.gesture_height       = 0.30;  % 【关键】假设手势发生的物理高度 (米)
TRAJ_PARA.energy_threshold_ratio = 0.4; % 能量阈值
TRAJ_PARA.min_hit_sats         = 2;     % 最少卫星数
TRAJ_PARA.miss_conflict_dist   = 0.15;  % 3D空间中的冲突距离 (米)
TRAJ_PARA.min_elevation        = 15;    % 最低仰角

% 绘制接收机 (原点)
plot3(ax, 0, 0, 0, '^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', 'Receiver');

% 绘制手势平面 (半透明参考面)
plane_range = 1.5; % 平面大小 +/- 1.5米
patch(ax, [-plane_range plane_range plane_range -plane_range], ...
          [-plane_range -plane_range plane_range plane_range], ...
          [TRAJ_PARA.gesture_height TRAJ_PARA.gesture_height TRAJ_PARA.gesture_height TRAJ_PARA.gesture_height], ...
          [0.9 0.9 1.0], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', 'Gesture Plane');

%% 3. 核心循环
colors = lines(length(segments)); 

fprintf('\n--> 开始 3D 轨迹估算...\n');

for i = 1:length(segments)
    seg = segments(i);
    idx_range = seg.start_idx : seg.end_idx;
    seg_times = t_grid(idx_range);
    
    fprintf('\n=== 手势片段 #%d (Peak: %s) ===\n', i, datestr(seg.peak_time, 'HH:MM:SS'));
    
    sub_volatility = volatility_matrix(idx_range, :);
    
    % 存储 3D 交点信息: [x, y, z]
    sat_points = struct('id', {}, 'pos_3d', {}, 'energy', {}, 'time_offset', {}, 'vec_u', {});
    num_pts = 0;
    
    % 获取接收机位置 (用于计算 ENU 向量)
    [~, epoch_idx] = min(abs([obs_data.time] - seg.peak_time));
    [rec_pos, ~, sat_states] = calculate_receiver_position(obs_data, nav_data, epoch_idx);
    if isempty(rec_pos), continue; end
    [rec_lat, rec_lon, rec_alt] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));

    % --- 3.1 计算每颗卫星的视线向量与平面的交点 ---
    for s = 1:length(valid_sats)
        s_id = valid_sats{s};
        if ~isfield(sat_states, s_id), continue; end
        
        % 能量 & 时间
        s_energy = sum(sub_volatility(:, s), 'omitnan');
        [~, local_max_idx] = max(sub_volatility(:, s));
        t_offset = seconds(seg_times(local_max_idx) - seg_times(1));
        
        % 计算 ENU 视线向量 (Unit Vector)
        sat_pos = sat_states.(s_id).position;
        [e, n, u] = ecef2enu(sat_pos(1)-rec_pos(1), sat_pos(2)-rec_pos(2), sat_pos(3)-rec_pos(3), rec_lat, rec_lon, rec_alt);
        
        dist = norm([e, n, u]);
        vec_u = [e, n, u] / dist; % 单位向量 [Ux, Uy, Uz]
        
        el_deg = asind(vec_u(3));
        if el_deg < TRAJ_PARA.min_elevation, continue; end 
        
        % --- 【核心物理建模】计算射线与平面 Z = h 的交点 ---
        % 射线方程: P = t * vec_u
        % 我们需要 P.z = h  =>  t * vec_u(3) = h  =>  t = h / vec_u(3)
        if vec_u(3) <= 0, continue; end % 忽略地平线以下的
        
        t_intersect = TRAJ_PARA.gesture_height / vec_u(3);
        intersect_point = t_intersect * vec_u; % [x, y, z] 交点坐标
        
        % 限制范围 (只考虑接收机附近的交点)
        if norm(intersect_point(1:2)) > 2.0, continue; end
        
        num_pts = num_pts + 1;
        sat_points(num_pts).id = s_id;
        sat_points(num_pts).pos_3d = intersect_point;
        sat_points(num_pts).vec_u = vec_u; % 保存视线方向用于画射线
        sat_points(num_pts).energy = s_energy;
        sat_points(num_pts).time_offset = t_offset;
    end
    
    if num_pts < 2, continue; end
    
    % --- 3.2 Hit / Miss 分类 ---
    all_energies = [sat_points.energy];
    threshold = max(all_energies) * TRAJ_PARA.energy_threshold_ratio;
    
    hits_mask = all_energies > threshold;
    hits = sat_points(hits_mask);
    misses = sat_points(~hits_mask);
    
    if length(hits) < TRAJ_PARA.min_hit_sats
        fprintf('   [跳过] 有效卫星不足\n');
        continue; 
    end
    
    % --- 3.3 在 3D 空间拟合直线 (实际是在 Z=h 平面上) ---
    % 提取 Hit 点的 (x, y) 坐标进行 PCA
    coords_3d = vertcat(hits.pos_3d); % N x 3
    P_xy = coords_3d(:, 1:2);         % N x 2 (只用 xy 拟合，z是固定的)
    
    mean_P = mean(P_xy);
    P_centered = P_xy - mean_P;
    
    [coeff, ~, ~] = pca(P_centered);
    dir_vec_2d = coeff(:, 1)'; % [dx, dy]
    
    % --- 3.4 时序定向 ---
    projections = P_centered * dir_vec_2d';
    times = [hits.time_offset]';
    corr_val = corr(projections, times);
    
    if ~isnan(corr_val) && corr_val < 0
        dir_vec_2d = -dir_vec_2d; 
    end
    
    % --- 3.5 计算 3D 轨迹端点 ---
    proj_final = (P_xy - mean_P) * dir_vec_2d';
    start_xy = mean_P + (min(proj_final) - 0.1) * dir_vec_2d;
    end_xy   = mean_P + (max(proj_final) + 0.1) * dir_vec_2d;
    
    % 恢复为 3D 坐标 (z = h)
    start_3d = [start_xy, TRAJ_PARA.gesture_height];
    end_3d   = [end_xy,   TRAJ_PARA.gesture_height];
    
    % --- 3.6 冲突检测 (3D 距离) ---
    % 检测 Miss 点到轨迹线段的距离
    conflict_count = 0;
    seg_vec = end_3d - start_3d;
    len_sq = dot(seg_vec, seg_vec);
    
    for m = 1:length(misses)
        m_pt = misses(m).pos_3d;
        if len_sq > 0
            t = dot(m_pt - start_3d, seg_vec) / len_sq;
            if t > 0 && t < 1
                dist = norm(m_pt - (start_3d + t * seg_vec));
                if dist < TRAJ_PARA.miss_conflict_dist
                    conflict_count = conflict_count + 1;
                    % 在 3D 图中标记冲突
                    plot3(ax, m_pt(1), m_pt(2), m_pt(3), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
                end
            end
        end
    end
    
    % --- 3.7 结果输出 ---
    traj_az = atan2d(dir_vec_2d(1), dir_vec_2d(2));
    if traj_az < 0, traj_az = traj_az + 360; end
    fprintf('   >>> 3D推演方向: %.1f 度 (相关性: %.2f)\n', traj_az, abs(corr_val));

    % --- 3.8 3D 绘图 ---
    draw_color = colors(mod(i-1, size(colors,1)) + 1, :);
    
    % 1. 画射线 (LoS Paths)
    % 仅对 Hit 卫星画彩色射线，Miss 卫星画淡灰色射线
    for k = 1:length(hits)
        pt = hits(k).pos_3d;
        % 从原点连线到交点
        plot3(ax, [0, pt(1)], [0, pt(2)], [0, pt(3)], '-', 'Color', [draw_color, 0.3], 'LineWidth', 1);
        % 画交点 (实心球)
        plot3(ax, pt(1), pt(2), pt(3), 'o', 'MarkerFaceColor', draw_color, 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
    end
    
    % 2. 画 Miss 点 (灰色小点，不画射线以免太乱)
    if ~isempty(misses)
        m_pts = vertcat(misses.pos_3d);
        plot3(ax, m_pts(:,1), m_pts(:,2), m_pts(:,3), '.', 'Color', [0.8 0.8 0.8], 'MarkerSize', 5);
    end
    
    % 3. 画手势轨迹 (粗箭头)
    quiver3(ax, start_3d(1), start_3d(2), start_3d(3), ...
            end_3d(1)-start_3d(1), end_3d(2)-start_3d(2), end_3d(3)-start_3d(3), ...
            0, 'Color', draw_color, 'LineWidth', 3, 'MaxHeadSize', 0.5);
        
    % 4. 文本标签 (浮在轨迹上方)
    text(ax, end_3d(1), end_3d(2), end_3d(3) + 0.05, sprintf('#%d: %.0f^o', i, traj_az), ...
        'Color', draw_color, 'FontWeight', 'bold', 'FontSize', 10, 'BackgroundColor', 'w');
end

hold(ax, 'off');
fprintf('\n✅ 3D 视距路径建模与轨迹拟合完成！\n');
fprintf('提示: 请在弹出的图形窗口中使用鼠标旋转视角，以观察空间关系。\n');













% % =========================================================================
% % step2_direction_estimation_v8_final.m
% % 功能: 手势方向估计 (完全匹配用户天空图风格版)
% % 特性:
% %   1. 【底图】使用标准 Zenith-Center 天空图 (与通常习惯一致，天顶在中心)。
% %   2. 【配色】完全复刻 calculate_and_plot_all_skyplot.m 中的 G/C/E/J 系统配色。
% %   3. 【内核】保持 Gnomonic 投影拟合直线的严谨数学逻辑。
% % =========================================================================
% 
% %% 1. 检查环境
% clearvars -except obs_data nav_data segments volatility_matrix valid_sats t_grid PARA;
% if ~exist('segments', 'var') || isempty(segments)
%     error('错误: 未找到 segments 变量。请先运行 step1_segmentation_GVI.m');
% end
% 
% addpath(genpath('sky_plot')); 
% addpath(genpath('calculate_clock_bias_and_positon'));
% addpath(genpath('nav_parse'));
% 
% %% 2. 初始化天空图画布 (风格对齐)
% fprintf('--> 正在初始化天空图轨迹分析...\n');
% fig_handle = figure('Name', 'Gesture Direction Analysis', 'Position', [100, 100, 1000, 800]);
% 
% pax = polaraxes('Parent', fig_handle);
% hold(pax, 'on');
% 
% % --- 样式设置 (对齐 calculate_and_plot_all_skyplot.m) ---
% set(pax, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise');
% % 注意：这里使用标准的正向坐标 (天顶在中心 r=0)，而不是 reverse
% set(pax, 'RLim', [0 90]); 
% set(pax, 'RTick', [0 30 60 90]);
% % 标签对齐：中心是90度，边缘是0度
% set(pax, 'RTickLabel', {'90', '60', '30', '0'}); 
% pax.RAxis.Label.String = '仰角 (°)';
% title(pax, '手势方向推断 (基于视距切割模型)');
% 
% %% 3. 参数设置
% TRAJ_PARA.energy_threshold_ratio = 0.4; 
% TRAJ_PARA.min_hit_sats         = 2;    
% TRAJ_PARA.miss_conflict_dist   = 0.25;  
% TRAJ_PARA.min_elevation        = 15;   
% 
% %% 4. 系统配色定义 (复刻用户代码)
% get_sys_color = @(sid) ...
%     (startsWith(sid,'G') * [0 0.4470 0.7410]) + ... % GPS: Blue
%     (startsWith(sid,'C') * [0.8500 0.3250 0.0980]) + ... % BDS: Orange
%     (startsWith(sid,'E') * [0.4660 0.6740 0.1880]) + ... % GAL: Green
%     (startsWith(sid,'J') * [0.4940 0.1840 0.5560]) + ... % QZSS: Purple
%     (~ismember(sid(1), 'GCEJ') * [0.5 0.5 0.5]);         % Others: Gray
% 
% %% 5. 核心循环
% traj_colors = lines(length(segments)); % 轨迹线的颜色 (区分不同手势)
% 
% fprintf('\n--> 开始方向估算...\n');
% 
% for i = 1:length(segments)
%     seg = segments(i);
%     idx_range = seg.start_idx : seg.end_idx;
%     seg_times = t_grid(idx_range);
%     
%     fprintf('\n=== 手势片段 #%d (Peak: %s) ===\n', i, datestr(seg.peak_time, 'HH:MM:SS'));
%     
%     sub_volatility = volatility_matrix(idx_range, :);
%     
%     sat_points = struct('id', {}, 'az_deg', {}, 'el_deg', {}, 'gx', {}, 'gy', {}, 'energy', {}, 'time_offset', {});
%     num_pts = 0;
%     
%     [~, epoch_idx] = min(abs([obs_data.time] - seg.peak_time));
%     [rec_pos, ~, sat_states] = calculate_receiver_position(obs_data, nav_data, epoch_idx);
%     if isempty(rec_pos), continue; end
%     [rec_lat, rec_lon, rec_alt] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));
% 
%     for s = 1:length(valid_sats)
%         s_id = valid_sats{s};
%         if ~isfield(sat_states, s_id), continue; end
%         
%         s_energy = sum(sub_volatility(:, s), 'omitnan');
%         [~, local_max_idx] = max(sub_volatility(:, s));
%         t_offset = seconds(seg_times(local_max_idx) - seg_times(1));
%         
%         sat_pos = sat_states.(s_id).position;
%         [e, n, u] = ecef2enu(sat_pos(1)-rec_pos(1), sat_pos(2)-rec_pos(2), sat_pos(3)-rec_pos(3), rec_lat, rec_lon, rec_alt);
%         az_rad = atan2(e, n); 
%         el_deg = asind(u / norm([e, n, u]));
%         az_deg = rad2deg(az_rad); if az_deg<0, az_deg=az_deg+360; end
%         
%         if el_deg < TRAJ_PARA.min_elevation, continue; end 
%         
%         % Gnomonic 投影 (数学拟合用)
%         R_math = 1.0 / tand(el_deg);
%         if R_math > 5.0, R_math = 5.0; end
%         g_x = R_math * sin(az_rad);
%         g_y = R_math * cos(az_rad);
%         
%         num_pts = num_pts + 1;
%         sat_points(num_pts).id = s_id;
%         sat_points(num_pts).az_deg = az_deg;
%         sat_points(num_pts).el_deg = el_deg;
%         sat_points(num_pts).gx = g_x;
%         sat_points(num_pts).gy = g_y;
%         sat_points(num_pts).energy = s_energy;
%         sat_points(num_pts).time_offset = t_offset;
%     end
%     
%     if num_pts < 2, continue; end
%     
%     % --- Hit / Miss 分类 ---
%     all_energies = [sat_points.energy];
%     threshold = max(all_energies) * TRAJ_PARA.energy_threshold_ratio;
%     
%     hits = sat_points(all_energies > threshold);
%     misses = sat_points(all_energies <= threshold);
%     
%     if length(hits) < TRAJ_PARA.min_hit_sats
%         fprintf('   [跳过] 有效波动卫星不足\n');
%         continue; 
%     end
%     
%     % --- 直线拟合 (PCA in Gnomonic Space) ---
%     P_hit = [ [hits.gx]', [hits.gy]' ]; 
%     mean_P = mean(P_hit);
%     P_centered = P_hit - mean_P;
%     
%     [coeff, ~, ~] = pca(P_centered);
%     dir_vec = coeff(:, 1)'; 
%     
%     % --- 时序定向 ---
%     projections = P_centered * dir_vec';
%     times = [hits.time_offset]';
%     corr_val = corr(projections, times);
%     
%     if ~isnan(corr_val) && corr_val < 0
%         dir_vec = -dir_vec; 
%     end
%     
%     % --- 计算起点终点 (Gnomonic) ---
%     proj_final = (P_hit - mean_P) * dir_vec';
%     g_start = mean_P + (min(proj_final) - 0.2) * dir_vec;
%     g_end   = mean_P + (max(proj_final) + 0.2) * dir_vec;
%     
%     % --- 转换回 Skyplot 坐标 ---
%     convert_to_skyplot = @(x, y) deal(atan2(x, y), 90 - atand(1.0 / sqrt(x^2 + y^2)));
%     [th_start, rho_start] = convert_to_skyplot(g_start(1), g_start(2));
%     [th_end, rho_end]     = convert_to_skyplot(g_end(1), g_end(2));
%     
%     % --- 冲突检测 ---
%     line_vec = g_end - g_start;
%     line_len_sq = dot(line_vec, line_vec);
%     conflict_found = false;
%     for m = 1:length(misses)
%         m_pt = [misses(m).gx, misses(m).gy];
%         if line_len_sq > 0
%             t = dot(m_pt - g_start, line_vec) / line_len_sq;
%             if t > 0 && t < 1
%                 dist = norm(m_pt - (g_start + t * line_vec));
%                 if dist < TRAJ_PARA.miss_conflict_dist
%                     % 画冲突标记 (红色叉叉)
%                     polarplot(pax, deg2rad(misses(m).az_deg), 90 - misses(m).el_deg, 'rx', 'MarkerSize', 10, 'LineWidth', 1.5);
%                     conflict_found = true;
%                 end
%             end
%         end
%     end
%     if conflict_found, fprintf('   [警告] 存在路径冲突\n'); end
% 
%     % --- 结果输出 ---
%     traj_az = rad2deg(atan2(dir_vec(1), dir_vec(2)));
%     if traj_az < 0, traj_az = traj_az + 360; end
%     fprintf('   >>> 推测方向: %.1f 度 (相关性: %.2f)\n', traj_az, abs(corr_val));
% 
%     % --- 绘图 ---
%     draw_color = traj_colors(mod(i-1, size(traj_colors,1)) + 1, :);
%     
%     % 1. 画 Hit 点 (使用系统专属颜色)
%     for k = 1:length(hits)
%         sys_col = get_sys_color(hits(k).id);
%         % 用实心圆点表示 Hit
%         polarplot(pax, deg2rad(hits(k).az_deg), 90 - hits(k).el_deg, 'o', ...
%             'MarkerFaceColor', sys_col, 'MarkerEdgeColor', 'k', 'MarkerSize', 9);
%     end
%     
%     % 2. 画轨迹箭头 (颜色区分手势 ID)
%     polarplot(pax, [th_start, th_end], [rho_start, rho_end], '-', 'LineWidth', 3, 'Color', draw_color);
%     polarplot(pax, th_end, rho_end, '^', 'MarkerSize', 8, 'MarkerFaceColor', draw_color, 'MarkerEdgeColor', 'k');
%         
%     % 3. 标签
%     text(pax, th_end, rho_end, sprintf('  #%d', i), 'Color', draw_color, 'FontWeight', 'bold', 'FontSize', 12, 'BackgroundColor', 'w');
% end
% 
% hold(pax, 'off');
% fprintf('\n✅ 天空图分析完成 (Standard View)\n');






