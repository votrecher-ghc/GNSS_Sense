% % ============== calculate_receiver_position.m (集成经纬度转换并内置打印) ==============
% function [receiver_pos, receiver_clk_err_sec, available_sat_states, lat_deg, lon_deg, alt_m] = calculate_receiver_position(obs_data, nav_data, epoch_idx)
%     % --- 函数输出参数增加了 lat_deg, lon_deg, alt_m ---
% 
%     c = 299792458.0; max_iterations = 10; tolerance = 1e-4;
%     receiver_pos_approx = [0; 0; 0]; receiver_clk_m_approx = 0;
%     
%     available_sat_states = struct(); % 只存储成功计算出状态的卫星
%     current_epoch = obs_data(epoch_idx);
%     t_obs = current_epoch.time;
%     observed_sats_all = fieldnames(current_epoch.data);
%     fprintf('--> 开始为历元 %d (时间: %s) 进行定位解算...\n', epoch_idx, datestr(t_obs));
%     % --- 第一次遍历：计算所有可用卫星的状态 ---
%     for i = 1:length(observed_sats_all)
%         sat_id = observed_sats_all{i};
%         if ~(startsWith(sat_id, 'G') || startsWith(sat_id, 'C')), continue; end
%         
%         sat_obs = current_epoch.data.(sat_id);
%         pr_code = '';
%         if startsWith(sat_id, 'G') && isfield(sat_obs.pseudorange, 'C1C'), pr_code = 'C1C'; end
%         if startsWith(sat_id, 'C') && isfield(sat_obs.pseudorange, 'C2I'), pr_code = 'C2I'; end
%         if isempty(pr_code), continue; end
%         
%         pseudorange = sat_obs.pseudorange.(pr_code);
%         if isnan(pseudorange), continue; end
%         try
%             [sat_pos, sat_vel, sat_clk_err] = calculate_satellite_state(t_obs, pseudorange, sat_id, nav_data);
%             available_sat_states.(sat_id).position = sat_pos;
%             available_sat_states.(sat_id).velocity = sat_vel;
%             available_sat_states.(sat_id).clock_error = sat_clk_err;
%         catch
%             % 如果某个卫星计算失败（如缺星历），则不加入解算，并给出提示
%             fprintf('⚠️  跳过卫星 %s 的状态计算。\n', sat_id);
%             continue;
%         end
%     end
%     
%     % --- 第二次遍历：使用状态可用的卫星进行定位 ---
%     sats_for_positioning = fieldnames(available_sat_states);
%     if length(sats_for_positioning) < 4
%         warning('‼️ 无法完成定位：历元 %d 中可用卫星数量不足4颗（实际 %d 颗）。', epoch_idx, length(sats_for_positioning));
%         receiver_pos = NaN(3,1); receiver_clk_err_sec = NaN;
%         lat_deg = NaN; lon_deg = NaN; alt_m = NaN;
%         return;
%     end
%     
%     for iter = 1:max_iterations
%         H = []; L = []; 
%         for i = 1:length(sats_for_positioning)
%             sat_id = sats_for_positioning{i};
%             sat_state = available_sat_states.(sat_id);
%             sat_pos = sat_state.position;
%             sat_clk_err = sat_state.clock_error;
%             if startsWith(sat_id, 'G'), pr_code = 'C1C'; else, pr_code = 'C2I'; end
%             pseudorange = current_epoch.data.(sat_id).pseudorange.(pr_code);
%             geom_range = norm(sat_pos - receiver_pos_approx);
%             l_i = pseudorange - (geom_range + receiver_clk_m_approx - c * sat_clk_err);
%             los_vec = (sat_pos - receiver_pos_approx) / geom_range;
%             h_i = [-los_vec', 1];
%             H(i, :) = h_i;
%             L(i, 1) = l_i;
%         end
%         dx = H \ L;
%         correction_pos = dx(1:3);
%         receiver_pos_approx = receiver_pos_approx + correction_pos;
%         receiver_clk_m_approx = receiver_clk_m_approx + dx(4);
%         fprintf('    迭代 %d: 位置修正量 = %.4f m\n', iter, norm(correction_pos));
%         if norm(correction_pos) < tolerance
%             fprintf('    迭代收敛!\n');
%             break;
%         end
%     end
%     receiver_pos = receiver_pos_approx;
%     receiver_clk_err_sec = receiver_clk_m_approx / c;
%     fprintf('--> 定位解算完成。\n');
%     
%     % --- 将计算出的ECEF坐标转换为大地坐标 (经纬度、高程) ---
%     fprintf('--> 开始将ECEF坐标转换为大地坐标...\n');
%     a = 6378137.0; f = 1 / 298.257223563; e_sq = f * (2 - f);
%     X = receiver_pos(1); Y = receiver_pos(2); Z = receiver_pos(3);
%     lon_rad = atan2(Y, X);
%     p = sqrt(X^2 + Y^2);
%     lat_rad_old = atan2(Z, p * (1 - e_sq));
%     tolerance_lat = 1e-12; delta_lat = tolerance_lat + 1;
%     while delta_lat > tolerance_lat
%         N = a / sqrt(1 - e_sq * sin(lat_rad_old)^2);
%         lat_rad_new = atan2(Z + N * e_sq * sin(lat_rad_old), p);
%         delta_lat = abs(lat_rad_new - lat_rad_old);
%         lat_rad_old = lat_rad_new;
%     end
%     lat_rad = lat_rad_new;
%     N = a / sqrt(1 - e_sq * sin(lat_rad)^2);
%     alt_m = (p / cos(lat_rad)) - N;
%     lat_deg = rad2deg(lat_rad);
%     lon_deg = rad2deg(lon_rad);
%     fprintf('--> 坐标转换完成。\n\n');
%     
%     % --- 新增：在函数内直接打印最终结果 ---
%     fprintf('****************** 历元 %d 定位结果 ******************\n', epoch_idx);
%     fprintf('  > ECEF (X, Y, Z)      : [%.3f, %.3f, %.3f] 米\n', receiver_pos(1), receiver_pos(2), receiver_pos(3));
%     fprintf('  > 大地坐标 (Lat, Lon, Alt): %.6f°, %.6f°, %.3f m\n', lat_deg, lon_deg, alt_m);
%     fprintf('  > 接收机钟差          : %.9f 秒\n', receiver_clk_err_sec);
%     fprintf('*********************************************************\n\n');
% end







% 
% % ============== calculate_receiver_position.m (简洁输出版) ==============
% function [receiver_pos, receiver_clk_err_sec, available_sat_states, lat_deg, lon_deg, alt_m] = calculate_receiver_position(obs_data, nav_data, epoch_idx)
%     % --- 函数输出参数增加了 lat_deg, lon_deg, alt_m ---
% 
%     c = 299792458.0; max_iterations = 10; tolerance = 1e-4;
%     receiver_pos_approx = [0; 0; 0]; receiver_clk_m_approx = 0;
%     
%     available_sat_states = struct(); % 只存储成功计算出状态的卫星
%     current_epoch = obs_data(epoch_idx);
%     t_obs = current_epoch.time;
%     observed_sats_all = fieldnames(current_epoch.data);
%     % fprintf('--> 开始为历元 %d (时间: %s) 进行定位解算...\n', epoch_idx, datestr(t_obs));
%     
%     % --- 第一次遍历：计算所有可用卫星的状态 ---
%     sats_for_positioning = {};
%     for i = 1:length(observed_sats_all)
%         sat_id = observed_sats_all{i};
%         if ~(startsWith(sat_id, 'G') || startsWith(sat_id, 'C')), continue; end
%         
%         sat_obs = current_epoch.data.(sat_id);
%         pr_code = '';
%         if startsWith(sat_id, 'G') && isfield(sat_obs.pseudorange, 'C1C'), pr_code = 'C1C'; end
%         if startsWith(sat_id, 'C') && isfield(sat_obs.pseudorange, 'C2I'), pr_code = 'C2I'; end
%         if isempty(pr_code), continue; end
%         
%         pseudorange = sat_obs.pseudorange.(pr_code);
%         if isnan(pseudorange), continue; end
%         try
%             % 使用最终修正版的 state 计算函数
%             [sat_pos, sat_vel, sat_clk_err] = calculate_satellite_state(t_obs, pseudorange, sat_id, nav_data);
%             available_sat_states.(sat_id).position = sat_pos;
%             available_sat_states.(sat_id).velocity = sat_vel;
%             available_sat_states.(sat_id).clock_error = sat_clk_err;
%             available_sat_states.(sat_id).pseudorange_raw = pseudorange; % 存储原始伪距
%             sats_for_positioning{end+1} = sat_id;
%         catch
%             % fprintf('⚠️  跳过卫星 %s 的状态计算。\n', sat_id);
%             continue;
%         end
%     end
%     
%     % --- 第二次遍历：使用状态可用的卫星进行定位 ---
%     sats_for_positioning = fieldnames(available_sat_states);
%     if length(sats_for_positioning) < 4
%         % warning('‼️ 无法完成定位：历元 %d 中可用卫星数量不足4颗（实际 %d 颗）。', epoch_idx, length(sats_for_positioning));
%         receiver_pos = NaN(3,1); receiver_clk_err_sec = NaN;
%         lat_deg = NaN; lon_deg = NaN; alt_m = NaN;
%         return;
%     end
%     
%     for iter = 1:max_iterations
%         H = []; L = []; 
%         for i = 1:length(sats_for_positioning)
%             sat_id = sats_for_positioning{i};
%             sat_state = available_sat_states.(sat_id);
%             sat_pos = sat_state.position;
%             sat_clk_err = sat_state.clock_error;
% 
%             pseudorange = sat_state.pseudorange_raw; % 使用存储的伪距
%             
%             % 地球自转改正
%             omega_e = 7.2921151467e-5;
%             travel_time = pseudorange / c;
%             R_e = [cos(omega_e * travel_time), sin(omega_e * travel_time), 0;
%                   -sin(omega_e * travel_time), cos(omega_e * travel_time), 0;
%                    0, 0, 1];
%             sat_pos_rotated = R_e * sat_pos;
%             
% %             geom_range = norm(sat_pos_rotated - receiver_pos_approx);
% %             
% %             % 伪距残差计算
% %             l_i = pseudorange - (geom_range + receiver_clk_m_approx - c * sat_clk_err);
% %             
% %             los_vec = (sat_pos_rotated - receiver_pos_approx) / geom_range;
% % ... 循环内部 ...
% 
% % ====================【新增：地球自转改正】====================
% omega_e = 7.2921151467e-5; % 地球自转角速度 (rad/s)
% travel_time = pseudorange / c; % 信号传播时间估算
% 
% % 构造地球自转矩阵
% R_e = [cos(omega_e * travel_time), sin(omega_e * travel_time), 0;
%       -sin(omega_e * travel_time), cos(omega_e * travel_time), 0;
%        0, 0, 1];
% 
% % 将卫星在“发射时刻”的位置，旋转到“接收时刻”的惯性系位置
% sat_pos_rotated = R_e * sat_pos;
% % =============================================================
% 
% % [修改] 使用旋转后的卫星位置来计算几何距离和方向向量
% geom_range = norm(sat_pos_rotated - receiver_pos_approx);
% l_i = pseudorange - (geom_range + receiver_clk_m_approx - c * sat_clk_err);
% los_vec = (sat_pos_rotated - receiver_pos_approx) / geom_range;
% % ...
% 
% 
% 
%             h_i = [-los_vec', 1];
%             H(i, :) = h_i;
%             L(i, 1) = l_i;
%         end
%         dx = H \ L;
%         receiver_pos_approx = receiver_pos_approx + dx(1:3);
%         receiver_clk_m_approx = receiver_clk_m_approx + dx(4);
%         
%         % fprintf('    迭代 %d: 位置修正量 = %.4f m\n', iter, norm(dx(1:3)));
%         if norm(dx(1:3)) < tolerance
%             % fprintf('    迭代收敛!\n');
%             break;
%         end
%     end
%     receiver_pos = receiver_pos_approx;
%     receiver_clk_err_sec = receiver_clk_m_approx / c;
%     % fprintf('--> 定位解算完成。\n');
%     
%     % --- 坐标转换 ---
%     a = 6378137.0; f = 1 / 298.257223563; e_sq = f * (2 - f);
%     X = receiver_pos(1); Y = receiver_pos(2); Z = receiver_pos(3);
%     lon_rad = atan2(Y, X);
%     p = sqrt(X^2 + Y^2);
%     lat_rad_old = atan2(Z, p * (1 - e_sq));
%     tolerance_lat = 1e-12; delta_lat = tolerance_lat + 1;
%     while delta_lat > tolerance_lat
%         N = a / sqrt(1 - e_sq * sin(lat_rad_old)^2);
%         lat_rad_new = atan2(Z + N * e_sq * sin(lat_rad_old), p);
%         delta_lat = abs(lat_rad_new - lat_rad_old);
%         lat_rad_old = lat_rad_new;
%     end
%     lat_rad = lat_rad_new;
%     N = a / sqrt(1 - e_sq * sin(lat_rad)^2);
%     alt_m = (p / cos(lat_rad)) - N;
%     lat_deg = rad2deg(lat_rad);
%     lon_deg = rad2deg(lon_rad);
%     % fprintf('--> 坐标转换完成。\n\n');
%     
%     % fprintf('****************** 历元 %d 定位结果 ******************\n', epoch_idx);
%     % fprintf('  > ECEF (X, Y, Z)      : [%.3f, %.3f, %.3f] 米\n', receiver_pos(1), receiver_pos(2), receiver_pos(3));
%     % fprintf('  > 大地坐标 (Lat, Lon, Alt): %.6f°, %.6f°, %.3f m\n', lat_deg, lon_deg, alt_m);
%     % fprintf('  > 接收机钟差          : %.9f 秒\n', receiver_clk_err_sec);
%     % fprintf('*********************************************************\n\n');
% end









% ============== calculate_receiver_position.m (新版: 支持 G/C/E/J 定位) ==============
function [receiver_pos, receiver_clk_err_sec, available_sat_states, lat_deg, lon_deg, alt_m] = calculate_receiver_position(obs_data, nav_data, epoch_idx)
    % 功能:
    %   - 使用 G, C, E, J 四系统的观测数据进行最小二乘定位。
    %   - 依赖: calculate_satellite_state.m (必须是支持 G/C/E/J 的新版本)
    
    c = 299792458.0; max_iterations = 10; tolerance = 1e-4;
    receiver_pos_approx = [0; 0; 0]; receiver_clk_m_approx = 0;
    
    available_sat_states = struct(); 
    current_epoch = obs_data(epoch_idx);
    t_obs = current_epoch.time;
    observed_sats_all = fieldnames(current_epoch.data);
    
    % --- 第一次遍历：计算所有可用卫星的状态 (G/C/E/J) ---
    sats_for_positioning = {};
    for i = 1:length(observed_sats_all)
        sat_id = observed_sats_all{i};
        sys = upper(sat_id(1));

        % 检查是否是G, C, E, J 系统
        if ~(sys == 'G' || sys == 'C' || sys == 'E' || sys == 'J')
            continue;
        end
        
        sat_obs = current_epoch.data.(sat_id);
        
        % 根据不同系统选择合适的伪距码
        % G (GPS L1)     -> C1C
        % E (Galileo E1) -> C1C
        % J (QZSS L1)    -> C1C
        % C (BDS B1I)    -> C2I
        pr_code = '';
        if (sys == 'G' || sys == 'E' || sys == 'J')
            if isfield(sat_obs.pseudorange, 'C1C')
                pr_code = 'C1C';
            % (如果您的 E/J 文件使用 C1X，可以在这里添加 elseif)
            end
        elseif sys == 'C'
            if isfield(sat_obs.pseudorange, 'C2I')
                pr_code = 'C2I';
            % (如果您的 BDS 文件使用 C1X，可以在这里添加 elseif)
            end
        end
        
        if isempty(pr_code), continue; end

        pseudorange = sat_obs.pseudorange.(pr_code);
        if isnan(pseudorange), continue; end
        
        try
            % !! 关键: 调用我们刚刚替换的、支持多系统的 'calculate_satellite_state.m'
            [sat_pos, sat_vel, sat_clk_err] = ...
                calculate_satellite_state(t_obs, pseudorange, sat_id, nav_data);
            
            available_sat_states.(sat_id).position = sat_pos;
            available_sat_states.(sat_id).velocity = sat_vel;
            available_sat_states.(sat_id).clock_error = sat_clk_err;
            available_sat_states.(sat_id).pseudorange_raw = pseudorange;
            sats_for_positioning{end+1} = sat_id;
        catch ME
            % fprintf('⚠️  在定位中跳过卫星 %s 的状态计算: %s\n', sat_id, ME.message);
            continue;
        end
    end
    
    % --- 第二次遍历：使用状态可用的卫星进行定位 ---
    sats_for_positioning = fieldnames(available_sat_states);
    if length(sats_for_positioning) < 4
        receiver_pos = NaN(3,1); receiver_clk_err_sec = NaN;
        lat_deg = NaN; lon_deg = NaN; alt_m = NaN;
        return;
    end
    
    for iter = 1:max_iterations
        H = []; L = []; 
        for i = 1:length(sats_for_positioning)
            sat_id = sats_for_positioning{i};
            sat_state = available_sat_states.(sat_id);
            sat_pos = sat_state.position;
            sat_clk_err = sat_state.clock_error;
            pseudorange = sat_state.pseudorange_raw;
            
            omega_e = 7.2921151467e-5;
            travel_time = pseudorange / c;
            R_e = [cos(omega_e * travel_time), sin(omega_e * travel_time), 0;
                  -sin(omega_e * travel_time), cos(omega_e * travel_time), 0;
                   0, 0, 1];
            sat_pos_rotated = R_e * sat_pos;

            geom_range = norm(sat_pos_rotated - receiver_pos_approx);
            l_i = pseudorange - (geom_range + receiver_clk_m_approx - c * sat_clk_err);
            los_vec = (sat_pos_rotated - receiver_pos_approx) / geom_range;

            h_i = [-los_vec', 1];
            H(i, :) = h_i;
            L(i, 1) = l_i;
        end
        
        try
            dx = H \ L;
        catch
            % 矩阵奇异，解算失败
            receiver_pos = NaN(3,1); receiver_clk_err_sec = NaN;
            lat_deg = NaN; lon_deg = NaN; alt_m = NaN;
            return;
        end
        
        receiver_pos_approx = receiver_pos_approx + dx(1:3);
        receiver_clk_m_approx = receiver_clk_m_approx + dx(4);
        
        if norm(dx(1:3)) < tolerance
            break;
        end
    end
    receiver_pos = receiver_pos_approx;
    receiver_clk_err_sec = receiver_clk_m_approx / c;
    
    % --- 坐标转换 (ECEF -> Geodetic) ---
    a = 6378137.0; f = 1 / 298.257223563; e_sq = f * (2 - f);
    X = receiver_pos(1); Y = receiver_pos(2); Z = receiver_pos(3);
    lon_rad = atan2(Y, X);
    p = sqrt(X^2 + Y^2);
    lat_rad_old = atan2(Z, p * (1 - e_sq));
    tolerance_lat = 1e-12; delta_lat = tolerance_lat + 1;
    while delta_lat > tolerance_lat
        N = a / sqrt(1 - e_sq * sin(lat_rad_old)^2);
        lat_rad_new = atan2(Z + N * e_sq * sin(lat_rad_old), p);
        delta_lat = abs(lat_rad_new - lat_rad_old);
        lat_rad_old = lat_rad_new;
    end
    lat_rad = lat_rad_new;
    N = a / sqrt(1 - e_sq * sin(lat_rad)^2);
    alt_m = (p / cos(lat_rad)) - N;
    lat_deg = rad2deg(lat_rad);
    lon_deg = rad2deg(lon_rad);
end