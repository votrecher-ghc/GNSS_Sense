% ============== parse_rinex_nav_multi_gnss.m (保留G=1, C=2) ==============
function ephemeris = parse_rinex_nav_multi_gnss(filepath)
% PARSE_RINEX_NAV_MULTI_GNSS - 解析RINEX 3.xx格式的导航文件。
%
% 功能:
%   - 解析 G, C, R, E, J 系统。
%   - 严格保持与旧版 `parse_rinex_nav_gps_bds.m` 的兼容性：
%     - G (GPS)     -> 存储在第 1 列 (sys_idx = 1)
%     - C (BeiDou)  -> 存储在第 2 列 (sys_idx = 2)
%   - 新增系统索引:
%     - R (GLONASS) -> 存储在第 3 列 (sys_idx = 3)
%     - E (Galileo) -> 存储在第 4 列 (sys_idx = 4)
%     - J (QZSS)    -> 存储在第 5 列 (sys_idx = 5)

% --- 初始化存储 ---
% 为最多63颗卫星和 5个系统 初始化主元胞数组
ephemeris = cell(63, 5);

% --- 打开文件 ---
fid = fopen(filepath, 'r');
if fid == -1
    error('无法打开文件: %s', filepath);
end

% --- 1. 读取文件头部分 ---
try
    while true
        line = fgetl(fid);
        if contains(line, 'END OF HEADER')
            break;
        end
        if ~ischar(line) % 检查文件末尾
             error('RINEX 导航文件头不完整或未找到 "END OF HEADER"。');
        end
    end
catch ME
    fclose(fid);
    error('读取RINEX导航文件头部时出错: %s', ME.message);
end

% --- 2. 读取数据部分 ---
while ~feof(fid)
    line = fgetl(fid);
    if isempty(line) || ~ischar(line) || length(line) < 3
        continue;
    end
    
    % 统一科学计数法符号
    line = strrep(line, 'D', 'E');
    
    % 获取系统标识符
    sys_char = upper(line(1));
    
    % 初始化星历结构体
    eph = struct();
    sys_idx = []; % 初始化系统索引
    
    try
        % 解析所有系统的共同部分：PRN 和 历元时间 (Toc)
        eph.PRN = str2double(line(2:3));
        if isnan(eph.PRN)
           % fprintf('警告: 无法解析PRN: %s。跳过此条记录。\n', line(1:3));
           continue; % 跳过此条损坏的记录
        end
        eph.Toc.Year   = str2double(line(5:8));
        eph.Toc.Month  = str2double(line(10:11));
        eph.Toc.Day    = str2double(line(13:14));
        eph.Toc.Hour   = str2double(line(16:17));
        eph.Toc.Minute = str2double(line(19:20));
        eph.Toc.Second = str2double(line(22:23));

        % --- 根据不同系统，使用不同的解析逻辑 ---
        switch sys_char
            
            % --- Case 1: GPS (G) ---
            case 'G'
                eph.System = 'GPS'; 
                sys_idx = 1; % 保持第 1 列
                
                % 第1行: af0, af1, af2
                eph.af0 = str2double(line(24:42));
                eph.af1 = str2double(line(43:61));
                eph.af2 = str2double(line(62:80));
                
                % 读取接下来的 7 行
                for i = 1:7
                    line = strrep(fgetl(fid), 'D', 'E');
                    switch i
                        case 1; eph.IODE = str2double(line(5:23)); eph.Crs = str2double(line(24:42)); eph.Delta_n = str2double(line(43:61)); eph.M0 = str2double(line(62:80));
                        case 2; eph.Cuc = str2double(line(5:23)); eph.e = str2double(line(24:42)); eph.Cus = str2double(line(43:61)); eph.sqrtA = str2double(line(62:80));
                        case 3; eph.Toe = str2double(line(5:23)); eph.Cic = str2double(line(24:42)); eph.OMEGA0 = str2double(line(43:61)); eph.Cis = str2double(line(62:80));
                        case 4; eph.i0 = str2double(line(5:23)); eph.Crc = str2double(line(24:42)); eph.omega = str2double(line(43:61)); eph.OMEGA_DOT = str2double(line(62:80));
                        case 5; eph.IDOT = str2double(line(5:23)); eph.Codes_on_L2 = str2double(line(24:42)); eph.GPS_Week = str2double(line(43:61)); eph.L2_P_data_flag = str2double(line(62:80));
                        case 6; eph.SV_accuracy = str2double(line(5:23)); eph.SV_health = str2double(line(24:42)); eph.TGD = str2double(line(43:61)); eph.IODC = str2double(line(62:80));
                        case 7; eph.Transmission_time = str2double(line(5:23)); eph.Fit_interval = str2double(line(24:42));
                    end
                end

            % --- Case 2: BeiDou (C) ---
            case 'C'
                eph.System = 'BeiDou';
                sys_idx = 2; % 保持第 2 列

                % 第1行: A0, A1, A2 (北斗的钟差参数)
                eph.A0 = str2double(line(24:42));
                eph.A1 = str2double(line(43:61));
                eph.A2 = str2double(line(62:80));
                
                % 读取接下来的 7 行
                for i = 1:7
                    line = strrep(fgetl(fid), 'D', 'E');
                    switch i
                        case 1; eph.IODE = str2double(line(5:23)); eph.Crs = str2double(line(24:42)); eph.Delta_n = str2double(line(43:61)); eph.M0 = str2double(line(62:80));
                        case 2; eph.Cuc = str2double(line(5:23)); eph.e = str2double(line(24:42)); eph.Cus = str2double(line(43:61)); eph.sqrtA = str2double(line(62:80));
                        case 3; eph.Toe = str2double(line(5:23)); eph.Cic = str2double(line(24:42)); eph.OMEGA0 = str2double(line(43:61)); eph.Cis = str2double(line(62:80));
                        case 4; eph.i0 = str2double(line(5:23)); eph.Crc = str2double(line(24:42)); eph.omega = str2double(line(43:61)); eph.OMEGA_DOT = str2double(line(62:80));
                        case 5; eph.IDOT = str2double(line(5:23)); eph.Spare1 = str2double(line(24:42)); eph.BDS_Week = str2double(line(43:61)); eph.Spare2 = str2double(line(62:80));
                        case 6; eph.SV_accuracy = str2double(line(5:23)); eph.SatH1 = str2double(line(24:42)); eph.TGD1 = str2double(line(43:61)); eph.TGD2 = str2double(line(62:80)); % B1I/B3I
                        case 7; eph.Transmission_time = str2double(line(5:23)); eph.AODC = str2double(line(24:42));
                    end
                end
                
            % --- Case 3: GLONASS (R) ---
            case 'R'
                eph.System = 'GLONASS';
                sys_idx = 3; % 新增第 3 列
                
                % 第1行: 卫星钟差(tau_n), 相对频率(gamma_n), 消息时间(tk)
                eph.tau_n = str2double(line(24:42));
                eph.gamma_n = str2double(line(43:61));
                eph.tk = str2double(line(62:80));
                
                % 读取接下来的 3 行
                for i = 1:3
                    line = strrep(fgetl(fid), 'D', 'E');
                    switch i
                        case 1 % 第2行: X, X_dot, X_ddot, health
                            eph.pos_x = str2double(line(5:23));
                            eph.vel_x = str2double(line(24:42));
                            eph.acc_x = str2double(line(43:61));
                            eph.health = str2double(line(62:80));
                        case 2 % 第3行: Y, Y_dot, Y_ddot, freq_num
                            eph.pos_y = str2double(line(5:23));
                            eph.vel_y = str2double(line(24:42));
                            eph.acc_y = str2double(line(43:61));
                            eph.freq_num = str2double(line(62:80));
                        case 3 % 第4行: Z, Z_dot, Z_ddot, age
                            eph.pos_z = str2double(line(5:23));
                            eph.vel_z = str2double(line(24:42));
                            eph.acc_z = str2double(line(43:61));
                            eph.age = str2double(line(62:80));
                    end
                end
                
            % --- Case 4: Galileo (E) ---
            case 'E'
                eph.System = 'Galileo'; 
                sys_idx = 4; % 新增第 4 列
                
                % 第1行: af0, af1, af2
                eph.af0 = str2double(line(24:42));
                eph.af1 = str2double(line(43:61));
                eph.af2 = str2double(line(62:80));
                
                % 读取接下来的 7 行
                for i = 1:7
                    line = strrep(fgetl(fid), 'D', 'E');
                    switch i
                        case 1; eph.IODE = str2double(line(5:23)); eph.Crs = str2double(line(24:42)); eph.Delta_n = str2double(line(43:61)); eph.M0 = str2double(line(62:80));
                        case 2; eph.Cuc = str2double(line(5:23)); eph.e = str2double(line(24:42)); eph.Cus = str2double(line(43:61)); eph.sqrtA = str2double(line(62:80));
                        case 3; eph.Toe = str2double(line(5:23)); eph.Cic = str2double(line(24:42)); eph.OMEGA0 = str2double(line(43:61)); eph.Cis = str2double(line(62:80));
                        case 4; eph.i0 = str2double(line(5:23)); eph.Crc = str2double(line(24:42)); eph.omega = str2double(line(43:61)); eph.OMEGA_DOT = str2double(line(62:80));
                        case 5; eph.IDOT = str2double(line(5:23)); eph.Data_sources = str2double(line(24:42)); eph.GAL_Week = str2double(line(43:61)); % E/FNAV week
                        case 6; eph.SV_accuracy = str2double(line(5:23)); eph.SV_health = str2double(line(24:42)); eph.TGD_E1E5a = str2double(line(43:61)); eph.TGD_E1E5b = str2double(line(62:80));
                        case 7; eph.Transmission_time = str2double(line(5:23)); eph.Spare = str2double(line(24:42)); % 最后一个字段在E中是Spare
                    end
                end
                
            % --- Case 5: QZSS (J) ---
            case 'J'
                eph.System = 'QZSS'; 
                sys_idx = 5; % 新增第 5 列
                
                % 第1行: af0, af1, af2
                eph.af0 = str2double(line(24:42));
                eph.af1 = str2double(line(43:61));
                eph.af2 = str2double(line(62:80));
                
                % 读取接下来的 7 行 (格式同GPS LNAV)
                for i = 1:7
                    line = strrep(fgetl(fid), 'D', 'E');
                    switch i
                        case 1; eph.IODE = str2double(line(5:23)); eph.Crs = str2double(line(24:42)); eph.Delta_n = str2double(line(43:61)); eph.M0 = str2double(line(62:80));
                        case 2; eph.Cuc = str2double(line(5:23)); eph.e = str2double(line(24:42)); eph.Cus = str2double(line(43:61)); eph.sqrtA = str2double(line(62:80));
                        case 3; eph.Toe = str2double(line(5:23)); eph.Cic = str2double(line(24:42)); eph.OMEGA0 = str2double(line(43:61)); eph.Cis = str2double(line(62:80));
                        case 4; eph.i0 = str2double(line(5:23)); eph.Crc = str2double(line(24:42)); eph.omega = str2double(line(43:61)); eph.OMEGA_DOT = str2double(line(62:80));
                        case 5; eph.IDOT = str2double(line(5:23)); eph.Codes_on_L2 = str2double(line(24:42)); eph.GPS_Week = str2double(line(43:61)); % QZSS周(同GPS)
                        case 6; eph.SV_accuracy = str2double(line(5:23)); eph.SV_health = str2double(line(24:42)); eph.TGD = str2double(line(43:61)); eph.IODC = str2double(line(62:80));
                        case 7; eph.Transmission_time = str2double(line(5:23)); eph.Fit_interval = str2double(line(24:42));
                    end
                end
        end % 结束 switch sys_char

        % --- 存储解析完成的星历 ---
        if ~isempty(sys_idx) && eph.PRN > 0 && eph.PRN <= size(ephemeris, 1)
            if isempty(ephemeris{eph.PRN, sys_idx})
                ephemeris{eph.PRN, sys_idx} = eph;
            else
                ephemeris{eph.PRN, sys_idx}(end+1) = eph;
            end
        end
        
    catch ME
        % 如果某条记录解析失败，打印错误并尝试继续解析下一条
        fprintf('警告: 解析卫星 %s (系统 %s) 时发生错误: %s。跳过此条记录。\n', line(1:3), sys_char, ME.message);
        % 如果是因为 fgetl 失败 (例如读到了文件末尾)，则退出
        if strcmp(ME.identifier, 'MATLAB:FileIO:InvalidFid') || feof(fid)
            break;
        end
    end
end % 结束 while ~feof(fid)

% --- 关闭文件 ---
fclose(fid);
fprintf('✅ 多系统导航文件解析完成 (G=1, C=2, R=3, E=4, J=5)。\n');
end