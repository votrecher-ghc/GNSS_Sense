function [valid_count, valid_sats_list] = count_high_elevation_satellites(obs_data, nav_data, zenith_threshold_deg, specific_epoch_idx)
% COUNT_HIGH_ELEVATION_SATELLITES 统计指定历元下【天顶角】符合条件的卫星
% 输入: zenith_threshold_deg 为最大天顶角阈值 (例如 60度，对应仰角30度)

    if nargin < 4 || isempty(specific_epoch_idx)
        specific_epoch_idx = floor(length(obs_data) / 2);
    end

    % 1. 获取接收机位置
    [rec_pos_ecef, ~, ~, rec_lat, rec_lon, ~] = calculate_receiver_position(obs_data, nav_data, specific_epoch_idx);
    
    current_epoch = obs_data(specific_epoch_idx);
    t_obs = current_epoch.time;
    sat_list = fieldnames(current_epoch.data);
    
    valid_count = 0;
    valid_sats_list = {};
    
    fprintf('\n%-8s | %-10s | %-10s | %-8s\n', 'SatID', 'Zenith(Deg)', 'Azimuth', 'Result');
    fprintf('----------------------------------------------\n');

    for i = 1:length(sat_list)
        sat_id = sat_list{i};
        sat_data = current_epoch.data.(sat_id);
        
        sys = upper(sat_id(1));
        pr_code = '';
        if (sys == 'G' || sys == 'E' || sys == 'J') && isfield(sat_data.pseudorange, 'C1C'), pr_code = 'C1C';
        elseif sys == 'C' && isfield(sat_data.pseudorange, 'C2I'), pr_code = 'C2I'; end
        
        if isempty(pr_code), continue; end
        pseudorange = sat_data.pseudorange.(pr_code);
        
        try
            [sat_pos_ecef, ~, ~] = calculate_satellite_state(t_obs, pseudorange, sat_id, nav_data);
        catch, continue; end
        
        diff_ecef = sat_pos_ecef - rec_pos_ecef;
        [E, N, U] = ecef2enu(diff_ecef(1), diff_ecef(2), diff_ecef(3), rec_lat, rec_lon, 0); 
        
        % ==========【修改点 1: 核心计算】==========
        range = norm([E, N, U]);
        % 原代码: ele_rad = asin(U / range);
        % 修改为: acos 计算天顶角 (0度在头顶)
        zen_rad = acos(U / range);  
        
        zen_deg = rad2deg(zen_rad);
        azi_rad = atan2(E, N);
        azi_deg = rad2deg(azi_rad);
        if azi_deg < 0, azi_deg = azi_deg + 360; end
        
        % ==========【修改点 2: 阈值判定】==========
        % 原逻辑: ele_deg >= threshold (仰角越高越好)
        % 新逻辑: zen_deg <= threshold (天顶角越小越好)
        if zen_deg <= zenith_threshold_deg
            valid_count = valid_count + 1;
            valid_sats_list{end+1} = sat_id;
            status = '✅ PASS';
        else
            status = '   OUT';
        end
        
        fprintf('%-8s | %6.2f deg | %6.2f deg | %s\n', sat_id, zen_deg, azi_deg, status);
    end
    fprintf('统计结果: 在历元 #%d，天顶角 < %.1f° 的卫星共有 %d 颗。\n', specific_epoch_idx, zenith_threshold_deg, valid_count);
end