% ============== test.m (最终调用脚本) ==============
clear; clc; 
close all;
obs_filepath = 'mixfingure_11_3_1.obs'; 
nav_filepath = 'arounds_11_3_1.nav'; 
% --- 2. 解析文件 ---
fprintf('--> 正在解析观测文件: %s\n', obs_filepath);
obs_data = parse_rinex_obs(obs_filepath);
fprintf('--> 正在解析导航文件: %s\n', nav_filepath);
nav_data = parse_rinex_nav_multi_gnss(nav_filepath);
% fprintf('\n✅ 文件解析全部完成\n\n');
calculate_and_plot_all_skyplot(obs_data, nav_data);

%%
step1_segmentation_GVI;
step2_direction_estimation;