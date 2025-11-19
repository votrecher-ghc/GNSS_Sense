% =========================================================================
% step2_direction_estimation_v4_trajectory.m
% 功能: 攻克难点2 - 基于“视距路径切割模型”的时空轨迹拟合
% 核心: 
%   1. 区分 Hit (波动) 和 Miss (静默) 卫星
%   2. 使用 PCA (主成分分析) 拟合 Hit 卫星的空间分布轴线
%   3. 利用每颗卫星的 Peak Time 差异确定运动矢量方向
% =========================================================================

%% 1. 检查环境
clearvars -except obs_data nav_data segments volatility_matrix valid_sats t_grid PARA;
if ~exist('segments', 'var') || isempty(segments)
    error('错误: 未找到 segments 变量。请先运行 step1_segmentation_GVI.m');
end

% 确保必要路径
addpath(genpath('sky_plot')); 
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

%% 2. 绘制背景天空图 (作为画布)
fprintf('--> 正在初始化天空图画布...\n');
fig_handle = figure('Name', 'Gesture Trajectory Analysis', 'Position', [100, 100, 1000, 800]);
% 调用您的天空图函数 (仅画背景网格，暂不画卫星点，后面我们自己画)
% 这里为了灵活控制，我们手动建立极坐标轴
pax = polaraxes('Parent', fig_handle);
title('手势时空轨迹推断 (基于视距切割模型)');
set(pax, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise'); % 北为0，顺时针
set(pax, 'RLim', [0 90]); % 半径表示天顶角 (0=天顶, 90=地平线)
hold(pax, 'on');

%% 3. 高级参数设置
TRAJ_PARA.energy_threshold_ratio = 0.4; % 相对阈值: 能量 > 最大星能量的 40% 算 Hit
TRAJ_PARA.min_hit_sats           = 2;   % 至少需要 2 颗星才能拟合直线
TRAJ_PARA.projection_radius      = 90;  % 投影平面的归一化半径

%% 4. 核心循环：逐段分析
colors = lines(length(segments)); % 生成不同颜色

fprintf('\n--> 开始进行时空轨迹拟合...\n');

for i = 1:length(segments)
    seg = segments(i);
    idx_range = seg.start_idx : seg.end_idx;
    seg_times = t_grid(idx_range);
    
    fprintf('\n=== 手势片段 #%d (%s) ===\n', i, datestr(seg.peak_time, 'HH:MM:SS'));
    
    % --- 4.1 计算每颗卫星的能量和具体的峰值时刻 ---
    sub_volatility = volatility_matrix(idx_range, :);
    
    sat_stats = struct('id', {}, 'energy', {}, 'peak_time_offset', {}, 'az', {}, 'el', {}, 'x', {}, 'y', {});
    num_stats = 0;
    
    % 获取参考时刻的接收机位置 (用于计算 Az/El)
    [~, epoch_idx] = min(abs([obs_data.time] - seg.peak_time));
    [rec_pos, ~, sat_states] = calculate_receiver_position(obs_data, nav_data, epoch_idx);
    if isempty(rec_pos), continue; end
    [rec_lat, rec_lon, rec_alt] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));

    % 遍历所有卫星
    for s = 1:length(valid_sats)
        s_id = valid_sats{s};
        if ~isfield(sat_states, s_id), continue; end
        
        % 1. 能量积分
        s_energy = sum(sub_volatility(:, s), 'omitnan');
        
        % 2. 寻找该卫星自己的峰值时刻 (Local Peak)
        [~, local_max_idx] = max(sub_volatility(:, s));
        % 转换为相对于段开始的时间偏移 (秒)
        t_offset = seconds(seg_times(local_max_idx) - seg_times(1));
        
        % 3. 计算方位角/仰角
        sat_pos = sat_states.(s_id).position;
        [e, n, u] = ecef2enu(sat_pos(1)-rec_pos(1), sat_pos(2)-rec_pos(2), sat_pos(3)-rec_pos(3), rec_lat, rec_lon, rec_alt);
        az_deg = atan2d(e, n); if az_deg < 0, az_deg = az_deg + 360; end
        el_deg = asind(u / norm([e, n, u]));
        
        if el_deg < 0, continue; end % 忽略地平线以下的
        
        % 4. 投影到 2D 平面 (用于 PCA 拟合)
        % 极坐标转笛卡尔坐标: 北为Y轴，东为X轴
        % 半径 r = 90 - Elevation (天顶在中心)
        r = 90 - el_deg;
        theta_rad = deg2rad(az_deg);
        
        % 注意: 数学坐标系通常 0度是右(东)，逆时针。GPS坐标系 0度是上(北)，顺时针。
        % GPS(Az, r) -> Cartesian(x, y)
        % x (East)  = r * sin(Az)
        % y (North) = r * cos(Az)
        pos_x = r * sind(az_deg);
        pos_y = r * cosd(az_deg);
        
        num_stats = num_stats + 1;
        sat_stats(num_stats).id = s_id;
        sat_stats(num_stats).energy = s_energy;
        sat_stats(num_stats).peak_time_offset = t_offset;
        sat_stats(num_stats).az = az_deg;
        sat_stats(num_stats).el = el_deg;
        sat_stats(num_stats).x  = pos_x;
        sat_stats(num_stats).y  = pos_y;
    end
    
    if num_stats < 2, continue; end
    
    % --- 4.2 区分 Hit (波动) 和 Miss (静默) ---
    all_energies = [sat_stats.energy];
    max_energy = max(all_energies);
    threshold = max_energy * TRAJ_PARA.energy_threshold_ratio;
    
    hit_mask = all_energies > threshold;
    miss_mask = ~hit_mask;
    
    hits = sat_stats(hit_mask);
    misses = sat_stats(miss_mask);
    
    fprintf('   Hit卫星数: %d, Miss卫星数: %d\n', length(hits), length(misses));
    
    if length(hits) < TRAJ_PARA.min_hit_sats
        fprintf('   警告: 有效波动卫星不足，无法拟合轨迹。\n');
        continue; 
    end
    
    % --- 4.3 轨迹拟合 (PCA) ---
    % 提取 Hit 卫星的 (x, y) 坐标
    P = [ [hits.x]', [hits.y]' ]; 
    
    % 中心化
    mean_P = mean(P);
    P_centered = P - mean_P;
    
    % PCA: 计算协方差矩阵的特征向量
    [coeff, ~, latent] = pca(P_centered);
    
    % 主成分方向 (第1主成分)
    dir_vec = coeff(:, 1)'; % [dx, dy]
    
    % --- 4.4 时序定向 (关键步骤) ---
    % 将所有 Hit 点投影到主成分轴上，看时间与位置的相关性
    % Projection = P_centered dot dir_vec
    projections = P_centered * dir_vec';
    times = [hits.peak_time_offset]';
    
    % 计算相关系数: 如果 位置越靠前的点 时间越晚，说明方向反了
    % 我们希望找到一个方向向量，使得沿该方向移动时，时间是增加的
    corr_time_space = corr(projections, times);
    
    final_dir = dir_vec;
    if isnan(corr_time_space) || corr_time_space < 0
        final_dir = -dir_vec; % 反转方向
    end
    
    % 计算轨迹的起点和终点 (用于画箭头)
    % 在轴线上延伸，覆盖投影范围
    min_proj = min(P_centered * final_dir');
    max_proj = max(P_centered * final_dir');
    
    start_pt = mean_P + (min_proj * final_dir);
    end_pt   = mean_P + (max_proj * final_dir);
    
    % --- 4.5 计算最终的角度 (0-360) ---
    % final_dir 是 [dx, dy] (East, North)
    % Azimuth = atan2(dx, dy)
    traj_az = atan2d(final_dir(1), final_dir(2));
    if traj_az < 0, traj_az = traj_az + 360; end
    
    fprintf('   >>> 推测手势方向: %.1f度 (基于 %d 颗星的时序相关性: %.2f)\n', ...
        traj_az, length(hits), corr_time_space);

    % --- 4.6 可视化绘制 ---
    draw_color = colors(mod(i-1, size(colors,1)) + 1, :);
    
    % 1. 画 Miss 卫星 (灰色叉叉)
    if ~isempty(misses)
        polarplot(pax, deg2rad([misses.az]), 90-[misses.el], 'x', ...
            'Color', [0.7 0.7 0.7], 'MarkerSize', 6, 'DisplayName', 'Miss');
    end
    
    % 2. 画 Hit 卫星 (圆点，颜色深浅代表时间早晚?)
    % 这里简单画大圆点
    polarplot(pax, deg2rad([hits.az]), 90-[hits.el], 'o', ...
        'MarkerFaceColor', draw_color, 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 8, 'DisplayName', 'Hit');
    
    % 3. 画轨迹箭头
    % 需要把 (x,y) 转回 (theta, r)
    [start_th, start_r] = cart2pol(start_pt(1), start_pt(2));
    [end_th, end_r]     = cart2pol(end_pt(1), end_pt(2));
    
    % cart2pol 返回的是数学极坐标 (逆时针，0在右)，需要转回地理极坐标
    % Theta_geo = pi/2 - Theta_math
    start_az_rad = pi/2 - start_th;
    end_az_rad   = pi/2 - end_th;
    
    % 绘制主轴线
    polarplot(pax, [start_az_rad, end_az_rad], [start_r, end_r], '-', ...
        'LineWidth', 3, 'Color', draw_color);
    
    % 绘制箭头头部 (简单画个点或者用 text 标记)
    polarplot(pax, end_az_rad, end_r, '^', ...
        'MarkerFaceColor', draw_color, 'MarkerSize', 10);
    
    % 4. 标注文字
    text(pax, end_az_rad, end_r, sprintf('  #%d: %.0f^o', i, traj_az), ...
        'Color', draw_color, 'FontSize', 12, 'FontWeight', 'bold', ...
        'BackgroundColor', 'w', 'EdgeColor', draw_color);
        
    % 5. (可选) 连接 Hit 点到拟合线的虚线，显示"拟合感"
    % (代码略，避免太乱)
    
end

hold(pax, 'off');
title(pax, sprintf('手势轨迹推演 (共%d个手势)', length(segments)));
fprintf('\n✅ 轨迹拟合完成！现在能够区分正反方向和具体轴线了。\n');










% % =========================================================================
% % step2_direction_estimation.m (终极修正版 v3)
% % 功能: 攻克难点2 - 在天空图背景上，叠加手势真实方向
% % 修正: 修复 ecef2enu 调用参数错误 (XYZ vs LLA)，确保点位精准重合
% % =========================================================================
% 
% %% 1. 检查环境
% if ~exist('segments', 'var') || isempty(segments)
%     error('错误: 未找到 segments 变量。请先运行 step1_segmentation_GVI.m');
% end
% 
% % 确保必要函数在路径中
% addpath(genpath('sky_plot')); 
% addpath(genpath('calculate_clock_bias_and_positon'));
% addpath(genpath('nav_parse'));
% 
% %% 2. 调用您的函数绘制背景天空图
% fprintf('--> 正在调用您的天空图函数作为背景...\n');
% 
% fig_handle = findobj('Type', 'figure', 'Name', '所有卫星天空图 (G/C/E/J)');
% if isempty(fig_handle)
%     calculate_and_plot_all_skyplot(obs_data, nav_data);
%     fig_handle = gcf;
% else
%     figure(fig_handle);
% end
% 
% pax = findobj(fig_handle, 'Type', 'polaraxes');
% if isempty(pax)
%     error('无法在图中找到 PolarAxes，请检查绘图函数是否正常。');
% end
% 
% % 强制纠正半径方向: 天顶(0)在中心，地平线(90)在边缘
% set(pax, 'RDir', 'normal'); 
% hold(pax, 'on'); 
% 
% %% 3. 参数设置
% PARA.top_n_sats = 3;  % 取前3颗波动最强的卫星
% 
% %% 4. 核心循环：逐段计算并叠加绘制
% % 颜色表
% colors = [1 0 0; 0 1 0; 0 0 1; 1 0 1; 0 1 1; 0.8 0.5 0]; 
% 
% fprintf('\n--> 开始叠加手势分析结果...\n');
% 
% for i = 1:length(segments)
%     seg = segments(i);
%     
%     % --- 4.1 计算受害卫星能量 ---
%     idx_range = seg.start_idx : seg.end_idx;
%     sub_volatility = volatility_matrix(idx_range, :);
%     sat_energies = sum(sub_volatility, 1, 'omitnan');
%     
%     % --- 4.2 筛选 Top-N ---
%     [sorted_energy, sort_idx] = sort(sat_energies, 'descend');
%     valid_top_count = min(PARA.top_n_sats, sum(sorted_energy > 0));
%     if valid_top_count == 0, continue; end
%     top_indices = sort_idx(1:valid_top_count);
%     
%     % --- 4.3 获取卫星位置 (Az, El) ---
%     [~, epoch_idx] = min(abs([obs_data.time] - seg.peak_time));
%     [rec_pos, rec_clk, sat_states] = calculate_receiver_position(obs_data, nav_data, epoch_idx);
%     
%     if isempty(rec_pos), continue; end
%     
%     % 【修正步骤 A】: 将接收机 ECEF 转换为 Geodetic (Lat, Lon, Alt)
%     % ecef2enu 需要 lat, lon 作为参考点，而不是 x, y, z
%     [rec_lat, rec_lon, rec_alt] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));
%     
%     az_list = []; el_list = []; weights = []; sat_names = {};
%     
%     fprintf('手势 #%d (时刻 %s):\n', i, datestr(seg.peak_time, 'HH:MM:SS'));
%     
%     for k = 1:length(top_indices)
%         s_idx = top_indices(k);
%         s_id = valid_sats{s_idx};
%         s_energy = sat_energies(s_idx);
%         
%         if isfield(sat_states, s_id)
%             sat_pos = sat_states.(s_id).position;
%             
%             % 【修正步骤 B】: 正确计算差值向量并调用 ecef2enu
%             vec = sat_pos - rec_pos; % 卫星指向接收机的向量 (反向?) 
%             % 不，skyplot代码用的是 sat - rec，得到的是从接收机指向卫星的向量(ENU)
%             
%             [e, n, u] = ecef2enu(vec(1), vec(2), vec(3), rec_lat, rec_lon, rec_alt);
%             
%             % 计算方位角 (GPS定义: atan2(e, n))
%             az_deg = atan2d(e, n);
%             if az_deg < 0, az_deg = az_deg + 360; end
%             
%             % 计算仰角
%             el_deg = asind(u / norm([e, n, u]));
%             
%             if el_deg >= 0 % 只处理地平线以上的
%                 az_list(end+1) = az_deg; %#ok<AGROW>
%                 el_list(end+1) = el_deg; %#ok<AGROW>
%                 weights(end+1) = s_energy; %#ok<AGROW>
%                 sat_names{end+1} = s_id; %#ok<AGROW>
%                 
%                 fprintf('   - 影响卫星: %s (Az:%.0f, El:%.0f)\n', s_id, az_deg, el_deg);
%             end
%         end
%     end
%     
%     if isempty(az_list), continue; end
%     
%     % --- 4.4 计算加权质心 (真实方向) ---
%     mean_el = sum(el_list .* weights) / sum(weights);
%     
%     % 方位角矢量加权 (基于 GPS 坐标系: N=cos, E=sin)
%     az_rads = deg2rad(az_list);
%     n_comp = sum(weights .* cos(az_rads));
%     e_comp = sum(weights .* sin(az_rads));
%     
%     mean_az = rad2deg(atan2(e_comp, n_comp)); % atan2(y, x) -> atan2(E, N)
%     if mean_az < 0, mean_az = mean_az + 360; end
%     
%     % --- 4.5 叠加绘图 ---
%     color_idx = mod(i-1, size(colors,1)) + 1;
%     draw_color = colors(color_idx, :);
%     
%     % 1. 绘制受害卫星 (小圆点) - 使用弧度
%     th_sats = deg2rad(az_list);
%     r_sats  = 90 - el_list;
%     polarplot(pax, th_sats, r_sats, '.', 'MarkerSize', 25, 'Color', draw_color);
%     
%     % 2. 绘制手势方向 (大五角星)
%     th_mean = deg2rad(mean_az);
%     r_mean  = 90 - mean_el;
%     polarplot(pax, th_mean, r_mean, 'p', 'MarkerSize', 20, 'LineWidth', 2, ...
%         'Color', draw_color, 'MarkerFaceColor', draw_color);
%     
%     % 3. 连线
%     for k = 1:length(th_sats)
%         polarplot(pax, [th_mean, th_sats(k)], [r_mean, r_sats(k)], '-', ...
%             'LineWidth', 1.0, 'Color', [draw_color 0.5]);
%     end
%     
%     % 4. 标注文字
%     text(pax, th_mean, r_mean, sprintf('  #%d: %.0f^o', i, mean_az), ...
%         'Color', draw_color, 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'w');
% end
% 
% hold(pax, 'off');
% fprintf('\n✅ 已修复 ecef2enu 调用参数错误，点位现在应该精准重合！\n');
