function rinex_data = parse_rinex_obs_advanced(filepath)
% 解析 RINEX3 观测文件（ 常规 3.xx）
% • 支持多行观测码定义
% • 每观测值 16 字符，从第 4 字符开始
% • 自动拼接 "主数值 + 尾部数字" 或直接取主数值
% • 正确保留时间的“秒.微秒”精度

    debug = false;  % 若需调试可改为 true
    fid = fopen(filepath, 'r');
    if fid == -1
        error('❌ 无法打开文件 %s', filepath);
    end

    % ---------- 1. 解析头部 ----------
    obs_types = struct();
    while true
        ln = fgetl(fid);
        if contains(ln,'END OF HEADER'); break; end
        if contains(ln,'SYS / # / OBS TYPES')
            sys = strtrim(ln(1));
            if isempty(sys) || ~isletter(sys), continue; end
            if isfield(obs_types, sys), continue; end 
            n_total = str2double(ln(4:6));
            codes = strsplit(strtrim(ln(8:60)));
            codes = codes(~cellfun('isempty',codes)); 
            while length(codes) < n_total
                next_line = fgetl(fid);
                more_codes = strsplit(strtrim(next_line(8:60)));
                more_codes = more_codes(~cellfun('isempty',more_codes)); 
                codes = [codes, more_codes];
            end
            obs_types.(sys) = codes;
        end
    end

    % ---------- 2. 正文 ----------
    rinex_data = struct([]);
    epoch_idx = 0;
    sats_to_read_in_epoch = 0;

    while ~feof(fid)
        ln = fgetl(fid);
        if ~ischar(ln) || isempty(strtrim(ln)), continue; end

        if ln(1) == '>'  % 新历元
            epoch_idx = epoch_idx + 1;
            t = sscanf(ln(3:35), '%f');  % 解析时间：年 月 日 时 分 秒.毫秒
            rinex_data(epoch_idx).time = datetime(t(1), t(2), t(3), ...
                                                  t(4), t(5), t(6), ...
                                                  'Format','yyyy-MM-dd HH:mm:ss.SSSSSS');
            rinex_data(epoch_idx).data = struct();
            sats_to_read_in_epoch = str2double(ln(33:35)); 
            continue;
        end

        if sats_to_read_in_epoch > 0
            sat_ln = ln;
            sys = strtrim(sat_ln(1));
            prn = str2double(sat_ln(2:3));
            sat_id = sprintf('%s%02d', sys, prn);

            if ~isfield(obs_types, sys), continue; end
            codes = obs_types.(sys);  
            n = numel(codes);
            values = NaN(1, n);

            % ---- 补齐行长度 ----
            min_len = 3 + 16 * n;
            if length(sat_ln) < min_len
                sat_ln = [sat_ln, repmat(' ', 1, min_len - length(sat_ln))];
            end

            % ---- 解析观测值 ----
            for j = 1:n
                s = 4 + (j - 1) * 16;
                e = s + 15;
                raw = sat_ln(s:e);

%                 if all(isspace(raw))
%                     val = NaN;
%                     num_part = ''; flag_part = '';
%                 else
%                     num_part  = strtrim(raw(1:14));
%                     flag_part = strtrim(raw(15:16));
%                     val = str2double([num_part, flag_part]);
%                     if isnan(val)
%                         val = str2double([num_part, '.', flag_part]);
%                     end
%                     if isnan(val)
%                         val = str2double(num_part);
%                     end
%                 end



                if all(isspace(raw))
                    val = NaN;
                else
                    % ✅ 仅提取前14位主数值，忽略 LLI/SSI 标志
                    num_part = strtrim(raw(1:14));
                    val = str2double(num_part);
                end

                values(j) = val;

                if debug
                    fprintf('Epoch %4d %-3s 码:%-4s raw=[%s] num=[%s] flag=[%s] → %g\n',...
                        epoch_idx, sat_id, codes{j}, raw, num_part, flag_part, val);
                end
            end

            % ---- 分类存储 ----
            for j = 1:n
                code = codes{j};
                val = values(j);
                if startsWith(code,'C')
                    rinex_data(epoch_idx).data.(sat_id).pseudorange.(code) = val;
                elseif startsWith(code,'L')
                    rinex_data(epoch_idx).data.(sat_id).carrier_phase.(code) = val;
                elseif startsWith(code,'D')
                    rinex_data(epoch_idx).data.(sat_id).doppler.(code) = val;
                elseif startsWith(code,'S')
                    rinex_data(epoch_idx).data.(sat_id).snr.(code) = val;
                end
            end

            sats_to_read_in_epoch = sats_to_read_in_epoch - 1;
        end
    end

    fclose(fid);
    fprintf('✅ 成功解析 %d 个历元，包含系统: %s\n', ...
        epoch_idx, strjoin(fieldnames(obs_types), ', '));
end



