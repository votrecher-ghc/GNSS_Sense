% ============== test.m (最终调用脚本) ==============
clear; clc; 
close all;
obs_filepath = 'mixfingure_11_25_4.obs'; 
nav_filepath = 'arounds_11_25_1.nav'; 
% --- 2. 解析文件 ---
fprintf('--> 正在解析观测文件: %s\n', obs_filepath);
obs_data = parse_rinex_obs(obs_filepath);
fprintf('--> 正在解析导航文件: %s\n', nav_filepath);
nav_data = parse_rinex_nav_multi_gnss(nav_filepath);
fprintf('\n✅ 文件解析全部完成\n\n');
calculate_and_plot_all_skyplot(obs_data, nav_data);
 %%
% 
%  analyze_gpsense_ultimate(obs_data, nav_data, 'C08');
 %%
% =========================================================================
% =======================  循环处理所有卫星 --- 时间=========================
% =========================================================================
% 步骤 4.1: 获取在观测数据中出现过的所有卫星的唯一列表
fprintf('--> 正在获取所有观测到的卫星列表...\n');
all_sats = list_observed_satellites(obs_data);

% ====================【核心修改：允许 E 和 J 卫星】====================
% 步骤 4.2: 筛选出我们支持分析的G/C/E/J卫星
% analyze_gpsense_ultimate 函数目前支持'G','C','E','J'开头的卫星
sats_to_process = all_sats(startsWith(all_sats, 'G') | ...
                             startsWith(all_sats, 'C') | ...
                             startsWith(all_sats, 'E') | ...
                             startsWith(all_sats, 'J'));
num_sats_to_process = length(sats_to_process);
fprintf('--> 将对以下 %d 颗 G/C/E/J 卫星进行GPSense分析:\n', num_sats_to_process);
% =========================【修改结束】=========================

disp(sats_to_process'); % 打印将要处理的卫星列表
% 步骤 4.3: 建立 for 循环，为列表中的每一颗卫星执行分析
for i = 1:num_sats_to_process
    % 从列表中获取当前要处理的卫星ID
    current_satellite_id = sats_to_process{i};
    
    fprintf('\n\n======================================================\n');
    fprintf('======     正在分析卫星: %s (%d / %d)     ======\n', ...
             current_satellite_id, i, num_sats_to_process);
    fprintf('======================================================\n');
    
    % 调用您的终极分析函数
    % 函数会为这颗卫星生成一套完整的分析图表
    try
        analyze_gpsense_ultimate(obs_data, nav_data, current_satellite_id);
    catch ME
        % 增加一个错误捕获，如果某颗卫星分析失败，程序不会中断
        % 而是会打印错误信息，然后继续分析下一颗卫星
        fprintf('!!!!!! 对卫星 %s 的分析失败 !!!!!!\n', current_satellite_id);
        fprintf('错误信息: %s\n', ME.message);
    end
    
    pause(1); % 暂停1秒
end
fprintf('\n\n✅✅✅ 所有指定卫星的分析已全部完成！ ✅✅✅\n');
 
 %%
% =========================================================================
% =======================  循环处理所有卫星 --- 采样点=========================
% =========================================================================
% % 步骤 4.1: 获取在观测数据中出现过的所有卫星的唯一列表
% fprintf('--> 正在获取所有观测到的卫星列表...\n');
% all_sats = list_observed_satellites(obs_data);
% 
% % 步骤 4.2: 筛选出我们支持分析的GPS和北斗卫星
% % analyze_gpsense_ultimate 函数目前只支持'G'和'C'开头的卫星
% sats_to_process = all_sats(startsWith(all_sats, 'G') | startsWith(all_sats, 'C'));
% num_sats_to_process = length(sats_to_process);
% fprintf('--> 将对以下 %d 颗GPS/北斗卫星进行GPSense分析:\n', num_sats_to_process);
% disp(sats_to_process'); % 打印将要处理的卫星列表
% 
% % 步骤 4.3: 建立 for 循环，为列表中的每一颗卫星执行分析
% for i = 1:num_sats_to_process
%     % 从列表中获取当前要处理的卫星ID
%     current_satellite_id = sats_to_process{i};
%     
%     fprintf('\n\n======================================================\n');
%     fprintf('======     正在分析卫星: %s (%d / %d)     ======\n', ...
%              current_satellite_id, i, num_sats_to_process);
%     fprintf('======================================================\n');
%     
%     % 调用您的终极分析函数
%     % 函数会为这颗卫星生成一套完整的分析图表
%     try
%         analyze_gpsense_by_sample(obs_data, nav_data, current_satellite_id);
%     catch ME
%         % 增加一个错误捕获，如果某颗卫星分析失败，程序不会中断
%         % 而是会打印错误信息，然后继续分析下一颗卫星
%         fprintf('!!!!!! 对卫星 %s 的分析失败 !!!!!!\n', current_satellite_id);
%         fprintf('错误信息: %s\n', ME.message);
%     end
%     
%     pause(1); % 暂停1秒
% end
% 
% fprintf('\n\n✅✅✅ 所有指定卫星的分析已全部完成！ ✅✅✅\n');