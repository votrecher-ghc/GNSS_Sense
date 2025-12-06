% % ============== calculate_and_plot_all_skyplot.m (新增钟差绘图功能) ==============
% function [all_sats_data] = calculate_and_plot_all_skyplot(obs_data, nav_data)
% % 该版本在优化版基础上，新增了接收机钟差变化的绘图和统计。
% 
% fprintf('--> 开始执行整合分析流程 (优化版+钟差绘图)...\n');
% 
% % --- 1. [优化] 遍历所有历元，进行定位解算并存储所有结果 ---
% num_epochs = length(obs_data);
% % 结果矩阵列: [lat, lon, alt, x, y, z, clock_bias_sec]
% positioning_results = NaN(num_epochs, 7); 
% success_count = 0;
% 
% fprintf('\n=========================================================\n');
% fprintf('           步骤 1: 对所有 %d 个历元进行定位解算            \n', num_epochs);
% fprintf('=========================================================\n');
% h_wait = waitbar(0, '正在进行定位解算...');
% 
% for epoch_idx = 1:length(obs_data)
%     waitbar(epoch_idx / num_epochs, h_wait, sprintf('定位解算中... (%d / %d)', epoch_idx, num_epochs));
%     try
%         % [修改] 接收钟差输出 receiver_clk_err
%         [ecef_pos, receiver_clk_err, ~, lat, lon, alt] = ...
%             calculate_receiver_position(obs_data, nav_data, epoch_idx);
%         
%         if ~any(isnan(ecef_pos))
%             % [修改] 将钟差结果也存储起来
%             positioning_results(epoch_idx, :) = [lat, lon, alt, ecef_pos', receiver_clk_err];
%         end
%     catch ME
%         fprintf('历元 %d 定位失败: %s\n', epoch_idx, ME.message);
%     end
% end
% close(h_wait);
% 
% % 清理掉失败的历元 (NaN行)
% positioning_results = positioning_results(~isnan(positioning_results(:,1)), :);
% success_count = size(positioning_results, 1);
% 
% fprintf('\n>>> 定位解算完成！成功率: %d / %d 个历元 <<<\n\n', success_count, num_epochs);
% 
% % --- 2. [优化] 打印平均结果并绘制所有结果图 ---
% if success_count > 1
%     % 计算并打印平均值
%     mean_ecef_pos = mean(positioning_results(:, 4:6), 1);
%     [mean_lat, mean_lon, mean_alt] = ecef2geodetic(mean_ecef_pos(1), mean_ecef_pos(2), mean_ecef_pos(3));
%     mean_clk_err = mean(positioning_results(:, 7)); % [新增] 计算平均钟差
% 
%     fprintf('****************** 平均定位结果 ******************\n');
%     fprintf('  > 平均纬度: %.6f°\n', mean_lat);
%     fprintf('  > 平均经度: %.6f°\n', mean_lon);
%     fprintf('  > 平均高程: %.3f m\n', mean_alt);
%     fprintf('  > 平均钟差: %.9f 秒\n', mean_clk_err); % [新增] 打印平均钟差
%     fprintf('**************************************************\n\n');
% 
%     % 绘制轨迹图 (使用大地坐标)
%     figure('Name', '接收机定位轨迹图');
%     geoplot(positioning_results(:,1), positioning_results(:,2), 'b.-', 'MarkerSize', 10);
%     geobasemap satellite;
%     title(sprintf('接收机定位轨迹图 (%d 个点)', success_count));
%     grid on;
% 
%     % 绘制高程图 (使用大地坐标)
%     figure('Name', '接收机高程变化图');
%     plot(positioning_results(:,3), 'm-o', 'LineWidth', 1.5, 'MarkerSize', 4);
%     title(sprintf('接收机高程变化 (%d 个点)', success_count));
%     xlabel('成功定位的历元序号');
%     ylabel('高程 (米)');
%     grid on;
%     
%     % ====================【新增：绘制钟差变化图】====================
%     figure('Name', '接收机钟差变化图');
%     plot(positioning_results(:,7), 'g-o', 'LineWidth', 1.5, 'MarkerSize', 4);
%     title(sprintf('接收机钟差变化 (%d 个点)', success_count));
%     xlabel('成功定位的历元序号');
%     ylabel('钟差 (秒)');
%     grid on;
%     % =================================================================
%     
%     % [优化] 将平均ECEF和平均大地坐标作为天空图计算的参考
%     receiver_pos_ecef = mean_ecef_pos';
%     lat = mean_lat;
%     lon = mean_lon;
% 
% else
%     error('‼️ 没有足够多的成功定位点，无法继续进行天空图计算。');
% end
% 
% % --- 3. [优化] 使用平均位置，为所有卫星计算天空图数据 ---
% fprintf('\n=========================================================\n');
% fprintf('           步骤 2: 计算所有卫星的天空图轨迹            \n');
% fprintf('=========================================================\n');
% % (这部分代码保持不变，省略以保持简洁)
% % ...
% all_sat_ids = [];
% for i = 1:length(obs_data)
%     all_sat_ids = [all_sat_ids; fieldnames(obs_data(i).data)];
% end
% unique_sat_ids = unique(all_sat_ids);
% unique_sat_ids = unique_sat_ids(startsWith(unique_sat_ids, 'G') | startsWith(unique_sat_ids, 'C'));
% fprintf('    在数据中共找到 %d 颗不同的GPS/北斗卫星。\n', length(unique_sat_ids));
% fprintf('    使用平均参考位置 (ECEF): [%.2f, %.2f, %.2f]\n', receiver_pos_ecef);
% all_sats_data = struct();
% plotted_sats_count = 0;
% for sat_idx = 1:length(unique_sat_ids)
%     sid = unique_sat_ids{sat_idx};
%     az = []; el = [];
%     for epoch_i = 1:length(obs_data)
%         epoch = obs_data(epoch_i);
%         if ~isfield(epoch.data, sid), continue; end
%         
%         sat_obs = epoch.data.(sid);
%         pr_code = '';
%         if sid(1) == 'G' && isfield(sat_obs.pseudorange, 'C1C'), pr_code = 'C1C'; end
%         if sid(1) == 'C' && isfield(sat_obs.pseudorange, 'C2I'), pr_code = 'C2I'; end
%         if isempty(pr_code) || ~isfield(sat_obs.pseudorange, pr_code) || isnan(sat_obs.pseudorange.(pr_code)), continue; end
%         
%         try
%             [sat_pos, ~, ~] = calculate_satellite_state(epoch.time, sat_obs.pseudorange.(pr_code), sid, nav_data);
%             vec = sat_pos - receiver_pos_ecef;
%             [e, n, u] = ecef2enu(vec(1), vec(2), vec(3), lat, lon, 0);
%             
%             az_i = atan2d(e, n); if az_i < 0, az_i = az_i + 360; end
%             el_i = asind(u / norm([e, n, u]));
%             
%             if el_i >= 0
%                 az(end+1) = az_i;
%                 el(end+1) = el_i;
%             end
%         catch
%             continue;
%         end
%     end
%     if ~isempty(az)
%         plotted_sats_count = plotted_sats_count + 1;
%         all_sats_data(plotted_sats_count).ID = sid;
%         all_sats_data(plotted_sats_count).azimuth = az;
%         all_sats_data(plotted_sats_count).elevation = el;
%     end
% end
% 
% % --- 4. 绘制天空图 ---
% if plotted_sats_count == 0
%     error('‼️ 没有卫星可用于绘制天空图。');
% end
% fprintf('--> 数据计算完成，开始绘制最终天空图...\n');
% figure('Name', '所有卫星天空图 (带起止点标记)');
% pax = polaraxes;
% set(pax, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise', 'RAxisLocation', 90, 'RDir', 'reverse');
% rlim(pax, [0 90]);
% rticks(pax, [0 30 60 90]);
% rticklabels(pax, {'90', '60', '30', '0'});
% pax.RAxis.Label.String = '仰角 (°)';
% title('所有可见GPS/北斗卫星天空轨迹图');
% hold(pax, 'on');
% plot_handles = gobjects(plotted_sats_count, 1);
% legend_entries = cell(1, plotted_sats_count);
% colors = lines(plotted_sats_count);
% for i = 1:plotted_sats_count
%     theta = deg2rad(all_sats_data(i).azimuth);
%     r = 90 - all_sats_data(i).elevation;
%     plot_handles(i) = polarplot(pax, theta, r, '-', 'LineWidth', 2.0, 'Color', colors(i,:));
%     polarplot(pax, theta(1), r(1), 'o', 'MarkerEdgeColor', colors(i,:), 'MarkerSize', 8, 'LineWidth', 1.5);
%     polarplot(pax, theta(end), r(end), 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:), 'MarkerSize', 8);
%     text(pax, theta(end), r(end), ['  ' all_sats_data(i).ID], ...
%         'FontSize', 10, 'FontWeight', 'bold', 'Color', colors(i,:));
%     legend_entries{i} = all_sats_data(i).ID;
% end
% legend(plot_handles, legend_entries, 'Location', 'eastoutside');
% hold(pax, 'off');
% fprintf('✅ 所有分析和绘图任务完成！\n\n');
% 
% end




% 
% 
% % ============== calculate_and_plot_all_skyplot.m (简洁输出版) ==============
% function [all_sats_data] = calculate_and_plot_all_skyplot(obs_data, nav_data)
%     % 该版本只在最后输出关键的汇总信息和图形。
%     
%     fprintf('--> 开始执行整合分析流程...\n');
%     
%     % --- 1. 对所有历元进行定位解算 ---
%     num_epochs = length(obs_data);
%     positioning_results = NaN(num_epochs, 7); 
%     
%     h_wait = waitbar(0, '正在对所有历元进行定位解算，请稍候...');
%     for epoch_idx = 1:length(obs_data)
%         waitbar(epoch_idx / num_epochs, h_wait);
%         try
%             [ecef_pos, receiver_clk_err, ~, lat, lon, alt] = ...
%                 calculate_receiver_position(obs_data, nav_data, epoch_idx);
%             if ~any(isnan(ecef_pos))
%                 positioning_results(epoch_idx, :) = [lat, lon, alt, ecef_pos', receiver_clk_err];
%             end
%         catch ME
%             % 静默失败
%         end
%     end
%     close(h_wait);
%     
%     positioning_results = positioning_results(~isnan(positioning_results(:,1)), :);
%     success_count = size(positioning_results, 1);
%     
%     % --- 2. 打印关键的平均结果并绘图 ---
%     if success_count > 1
%         mean_ecef_pos = mean(positioning_results(:, 4:6), 1);
%         [mean_lat, mean_lon, mean_alt] = ecef2geodetic(mean_ecef_pos(1), mean_ecef_pos(2), mean_ecef_pos(3));
%         mean_clk_err = mean(positioning_results(:, 7));
% 
%         fprintf('\n================== [最终分析结果] ==================\n');
%         fprintf('成功解算率: %d / %d 个历元\n', success_count, num_epochs);
%         fprintf('------------------ 平均定位结果 ------------------\n');
%         fprintf('  > 平均纬度: %.6f°\n', mean_lat);
%         fprintf('  > 平均经度: %.6f°\n', mean_lon);
%         fprintf('  > 平均高程: %.3f m\n', mean_alt);
%         fprintf('  > 平均钟差: %.9f 秒\n', mean_clk_err);
%         fprintf('====================================================\n\n');
% 
% %         figure('Name', '接收机定位轨迹图');
% %         geoplot(positioning_results(:,1), positioning_results(:,2), 'b.-', 'MarkerSize', 10);
% %         geobasemap satellite;
% %         title(sprintf('接收机定位轨迹图 (%d 个点)', success_count));
% %         grid on;
% 
% %         figure('Name', '接收机高程变化图');
% %         plot(positioning_results(:,3), 'm-o', 'LineWidth', 1.5, 'MarkerSize', 4);
% %         title(sprintf('接收机高程变化 (%d 个点)', success_count));
% %         xlabel('成功定位的历元序号');
% %         ylabel('高程 (米)');
% %         grid on;
% 
%         figure('Name', '接收机钟差变化图');
%         plot(positioning_results(:,7), '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
%         title(sprintf('接收机钟差变化 (%d 个点)', success_count));
%         xlabel('成功定位的历元序号');
%         ylabel('钟差 (秒)');
%         grid on;
%         
%         receiver_pos_ecef = mean_ecef_pos';
%         lat = mean_lat;
%         lon = mean_lon;
%     else
%         error('‼️ 没有足够多的成功定位点，无法继续进行天空图计算。');
%     end
% 
%     % --- 3. 计算并绘制天空图 ---
%     fprintf('--> 正在计算并绘制卫星天空图...\n');
%     % (此部分代码与您提供的版本完全相同，省略以保持简洁)
%     all_sat_ids = [];
%     for i = 1:length(obs_data)
%         all_sat_ids = [all_sat_ids; fieldnames(obs_data(i).data)];
%     end
%     unique_sat_ids = unique(all_sat_ids);
%     unique_sat_ids = unique_sat_ids(startsWith(unique_sat_ids, 'G') | startsWith(unique_sat_ids, 'C'));
%     all_sats_data = struct();
%     plotted_sats_count = 0;
%     for sat_idx = 1:length(unique_sat_ids)
%         sid = unique_sat_ids{sat_idx};
%         az = []; el = [];
%         for epoch_i = 1:length(obs_data)
%             epoch = obs_data(epoch_i);
%             if ~isfield(epoch.data, sid), continue; end
%             sat_obs = epoch.data.(sid);
%             pr_code = '';
%             if sid(1) == 'G' && isfield(sat_obs.pseudorange, 'C1C'), pr_code = 'C1C'; end
%             if sid(1) == 'C' && isfield(sat_obs.pseudorange, 'C2I'), pr_code = 'C2I'; end
%             if isempty(pr_code) || ~isfield(sat_obs.pseudorange, pr_code) || isnan(sat_obs.pseudorange.(pr_code)), continue; end
%             try
%                 [sat_pos, ~, ~] = calculate_satellite_state(epoch.time, sat_obs.pseudorange.(pr_code), sid, nav_data);
%                 vec = sat_pos - receiver_pos_ecef;
%                 [e, n, u] = ecef2enu(vec(1), vec(2), vec(3), lat, lon, 0);
%                 az_i = atan2d(e, n); if az_i < 0, az_i = az_i + 360; end
%                 el_i = asind(u / norm([e, n, u]));
%                 if el_i >= 0, az(end+1) = az_i; el(end+1) = el_i; end
%             catch
%                 continue;
%             end
%         end
%         if ~isempty(az)
%             plotted_sats_count = plotted_sats_count + 1;
%             all_sats_data(plotted_sats_count).ID = sid;
%             all_sats_data(plotted_sats_count).azimuth = az;
%             all_sats_data(plotted_sats_count).elevation = el;
%         end
%     end
%     if plotted_sats_count > 0
%         figure('Name', '所有卫星天空图 (带起止点标记)');
%         pax = polaraxes;
%         set(pax, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise', 'RAxisLocation', 90, 'RDir', 'reverse');
%         rlim(pax, [0 90]);
%         rticks(pax, [0 30 60 90]);
%         rticklabels(pax, {'90', '60', '30', '0'});
%         pax.RAxis.Label.String = '仰角 (°)';
%         title('所有可见GPS/北斗卫星天空轨迹图');
%         hold(pax, 'on');
%         plot_handles = gobjects(plotted_sats_count, 1);
%         legend_entries = cell(1, plotted_sats_count);
%         colors = lines(plotted_sats_count);
%         for i = 1:plotted_sats_count
%             theta = deg2rad(all_sats_data(i).azimuth);
% 
% 
% 
% %             r = 90 - all_sats_data(i).elevation;
%             r = all_sats_data(i).elevation;
% 
%             plot_handles(i) = polarplot(pax, theta, r, '-', 'LineWidth', 2.0, 'Color', colors(i,:));
%             polarplot(pax, theta(1), r(1), 'o', 'MarkerEdgeColor', colors(i,:), 'MarkerSize', 8, 'LineWidth', 1.5);
%             polarplot(pax, theta(end), r(end), 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:), 'MarkerSize', 8);
%             text(pax, theta(end), r(end), ['  ' all_sats_data(i).ID], 'FontSize', 10, 'FontWeight', 'bold', 'Color', colors(i,:));
%             legend_entries{i} = all_sats_data(i).ID;
%         end
%         legend(plot_handles, legend_entries, 'Location', 'eastoutside');
%         hold(pax, 'off');
%     end
%     fprintf('✅ 所有分析和绘图任务完成！\n\n');
% end










% ============== calculate_and_plot_all_skyplot.m (新版: 支持 G/C/E/J 绘图) ==============
function [all_sats_data] = calculate_and_plot_all_skyplot(obs_data, nav_data)
    % 功能:
    %   - 计算并绘制 G, C, E, J 四系统的卫星天空图。
    %   - 依赖: 
    %     1. calculate_receiver_position.m (必须是支持 G/C/E/J 的新版本)
    %     2. calculate_satellite_state.m (必须是支持 G/C/E/J 的新版本)
    
    fprintf('--> 开始执行整合分析流程 (支持 G/C/E/J)...\n');
    
    % --- 1. 对所有历元进行定位解算 ---
    num_epochs = length(obs_data);
    positioning_results = NaN(num_epochs, 7); 
    
    h_wait = waitbar(0, '正在对所有历元进行定位解算 (G/C/E/J)，请稍候...');
    for epoch_idx = 1:length(obs_data)
        waitbar(epoch_idx / num_epochs, h_wait);
        try
            % 调用我们刚刚替换的、支持 G/C/E/J 的定位函数
            [ecef_pos, receiver_clk_err, ~, lat, lon, alt] = ...
                calculate_receiver_position(obs_data, nav_data, epoch_idx);
            
            if ~any(isnan(ecef_pos))
                positioning_results(epoch_idx, :) = [lat, lon, alt, ecef_pos', receiver_clk_err];
            end
        catch ME
            % 静默失败
        end
    end
    close(h_wait);
    
    positioning_results = positioning_results(~isnan(positioning_results(:,1)), :);
    success_count = size(positioning_results, 1);
    
    % --- 2. 打印关键的平均结果并绘图 ---
    if success_count > 1
        mean_ecef_pos = mean(positioning_results(:, 4:6), 1);
        [mean_lat, mean_lon, mean_alt] = ecef2geodetic(mean_ecef_pos(1), mean_ecef_pos(2), mean_ecef_pos(3));
        mean_clk_err = mean(positioning_results(:, 7));

        fprintf('\n================== [最终分析结果] ==================\n');
        fprintf('成功解算率: %d / %d 个历元\n', success_count, num_epochs);
        fprintf('------------------ 平均定位结果 ------------------\n');
        fprintf('  > 平均纬度: %.6f°\n', mean_lat);
        fprintf('  > 平均经度: %.6f°\n', mean_lon);
        fprintf('  > 平均高程: %.3f m\n', mean_alt);
        fprintf('  > 平均钟差: %.9f 秒\n', mean_clk_err);
        fprintf('====================================================\n\n');

        figure('Name', '接收机钟差变化图 (G/C/E/J)');
        plot(positioning_results(:,7), '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
        title(sprintf('接收机钟差变化 (%d 个点)', success_count));
        xlabel('成功定位的历元序号');
        ylabel('钟差 (秒)');
        grid on;
        
        receiver_pos_ecef = mean_ecef_pos';
        lat = mean_lat;
        lon = mean_lon;
    else
        fprintf('‼️ 定位解算失败率过高 (成功历元 < 2)，无法计算天空图。\n');
        fprintf('   请检查: \n');
        fprintf('   1. 观测文件 (%s) 是否包含 G/C/E/J 卫星的 C1C 或 C2I 观测值。\n', obs_data(1).time); % 假设 obs_data 非空
        fprintf('   2. 导航文件 (%s) 是否包含对应的星历数据。\n', '...');
        error('‼️ 没有足够多的成功定位点，无法继续进行天空图计算。');
    end

    % --- 3. 计算并绘制天空图 ---
    fprintf('--> 正在计算并绘制卫星天空图 (G/C/E/J)...\n');
    
    all_sat_ids = [];
    for i = 1:length(obs_data)
        all_sat_ids = [all_sat_ids; fieldnames(obs_data(i).data)];
    end
    unique_sat_ids = unique(all_sat_ids);
    
    % 筛选 G, C, E, J 系统的卫星
    unique_sat_ids = unique_sat_ids(startsWith(unique_sat_ids, 'G') | ...
                                    startsWith(unique_sat_ids, 'C') | ...
                                    startsWith(unique_sat_ids, 'E') | ...
                                    startsWith(unique_sat_ids, 'J'));
                                    
    fprintf('    在数据中共找到 %d 颗不同的G/C/E/J卫星。\n', length(unique_sat_ids));
    
    all_sats_data = struct();
    plotted_sats_count = 0;
    
    for sat_idx = 1:length(unique_sat_ids)
        sid = unique_sat_ids{sat_idx};
        az = []; el = [];
        
        for epoch_i = 1:length(obs_data)
            epoch = obs_data(epoch_i);
            if ~isfield(epoch.data, sid), continue; end
            sat_obs = epoch.data.(sid);
            
            % 根据不同系统选择合适的伪距码 (逻辑同 calculate_receiver_position)
            pr_code = '';
            sys = sid(1);
            if (sys == 'G' || sys == 'E' || sys == 'J')
                if isfield(sat_obs.pseudorange, 'C1C')
                    pr_code = 'C1C';
                end
            elseif sys == 'C'
                if isfield(sat_obs.pseudorange, 'C2I')
                    pr_code = 'C2I';
                end
            end
            
            if isempty(pr_code) || ~isfield(sat_obs.pseudorange, pr_code) || isnan(sat_obs.pseudorange.(pr_code))
                continue; 
            end
            
            try
                % !! 关键: 调用支持多系统的 'calculate_satellite_state.m'
                [sat_pos, ~, ~] = calculate_satellite_state(epoch.time, sat_obs.pseudorange.(pr_code), sid, nav_data);
                
                vec = sat_pos - receiver_pos_ecef;
                [e, n, u] = ecef2enu(vec(1), vec(2), vec(3), lat, lon, 0);
                
                az_i = atan2d(e, n); if az_i < 0, az_i = az_i + 360; end
                el_i = asind(u / norm([e, n, u]));
                
                if el_i >= 0, az(end+1) = az_i; el(end+1) = el_i; end
            catch ME
                % fprintf('⚠️  在天空图计算中跳过卫星 %s : %s\n', sid, ME.message);
                continue;
            end
        end
        
        if ~isempty(az)
            plotted_sats_count = plotted_sats_count + 1;
            all_sats_data(plotted_sats_count).ID = sid;
            all_sats_data(plotted_sats_count).azimuth = az;
            all_sats_data(plotted_sats_count).elevation = el;
        end
    end
    
    if plotted_sats_count > 0
        figure('Name', '所有卫星天空图 (G/C/E/J)');
        pax = polaraxes;
        set(pax, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise', 'RAxisLocation', 90, 'RDir', 'reverse');
        rlim(pax, [0 90]);
        rticks(pax, [0 30 60 90]);
        rticklabels(pax, {'90', '60', '30', '0'});
        pax.RAxis.Label.String = '仰角 (°)';
        title('所有可见 G/C/E/J 卫星天空轨迹图');
        hold(pax, 'on');
        
        plot_handles = gobjects(plotted_sats_count, 1);
        legend_entries = cell(1, plotted_sats_count);
        
        for i = 1:plotted_sats_count
            % 确定卫星系统的颜色
            sys_col = [0.5 0.5 0.5]; % 默认灰色
            if startsWith(all_sats_data(i).ID, 'G'), sys_col = [0 0.4470 0.7410]; % 蓝
            elseif startsWith(all_sats_data(i).ID, 'C'), sys_col = [0.8500 0.3250 0.0980]; % 橙
            elseif startsWith(all_sats_data(i).ID, 'E'), sys_col = [0.4660 0.6740 0.1880]; % 绿
            elseif startsWith(all_sats_data(i).ID, 'J'), sys_col = [0.4940 0.1840 0.5560]; % 紫
            end
            
            theta = deg2rad(all_sats_data(i).azimuth);
%             r = 90 - all_sats_data(i).elevation;
               r = all_sats_data(i).elevation;
            
            plot_handles(i) = polarplot(pax, theta, r, '-', 'LineWidth', 2.0, 'Color', sys_col);
            polarplot(pax, theta(1), r(1), 'o', 'MarkerEdgeColor', sys_col, 'MarkerSize', 8, 'LineWidth', 1.5);
            polarplot(pax, theta(end), r(end), 'o', 'MarkerFaceColor', sys_col, 'MarkerEdgeColor', sys_col, 'MarkerSize', 8);
            text(pax, theta(end), r(end), ['  ' all_sats_data(i).ID], 'FontSize', 10, 'FontWeight', 'bold', 'Color', sys_col);
            legend_entries{i} = all_sats_data(i).ID;
        end
        legend(plot_handles, legend_entries, 'Location', 'eastoutside');
        hold(pax, 'off');
    else
        warning('没有卫星可用于绘制天空图 (G/C/E/J)。');
    end
    fprintf('✅ 所有分析和绘图任务完成！\n\n');
end
