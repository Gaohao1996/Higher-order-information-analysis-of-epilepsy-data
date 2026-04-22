function [drivers_red, drivers_syn, g_red, g_syn, ...
    g_red_surr, g_syn_surr, ...
    sig_idx_red, sig_idx_syn, ...
    p_values_red, p_values_syn, ...
    U_final, R_final, S_final, ...
    sig_channels_red, sig_channels_syn, sig_channels_all, ...
    delayZ_red, delayZ_syn] = ...
    TE_syn_red_greedy_surr_full(x, i, j, p, d, phase, nsurr, alpha, delay_grid, max_Z_keep)
% 两阶段低开销版 + “历史 delayZ 累积保留”
% 先按原逻辑选 Z（不扫 delayZ），再仅对该 Z 扫 delay_grid。
% 已选Z的 delayZ 历史在冗余/协同两条路径分别累积保存，并用于置换与后续计算。

    if nargin < 10 || isempty(delay_grid), delay_grid = 1; end
    if nargin < 11 || isempty(max_Z_keep), max_Z_keep = inf; end

    % ---------- 基线 ----------
    base_set = [i j];
    TE0 = CTE_iteration(x(:, base_set), p, d, 1, 2, phase);  % 本地 1、2
    g0  = mean(TE0);

    % ---------- 预分配 ----------
    [~, N] = size(x);
    maxK = N;
    g_red  = nan(1, maxK);  g_red(1)  = g0;
    g_syn  = nan(1, maxK);  g_syn(1)  = g0;

    p_values_red = nan(1, maxK);  p_values_red(1) = 1;
    p_values_syn = nan(1, maxK);  p_values_syn(1) = 1;

    g_red_surr = nan(nsurr, maxK);
    g_syn_surr = nan(nsurr, maxK);

    sig_mask_red = false(1, maxK);
    sig_mask_syn = false(1, maxK);

    sig_channels_red = [];
    sig_channels_syn = [];

    % 累积保存：每条路径已选Z对应的 delayZ 历史
    delayZ_red = [];   % 与 sig_channels_red 对齐
    delayZ_syn = [];   % 与 sig_channels_syn 对齐

    drivers_red = base_set;
    drivers_syn = base_set;

    % ---------- 冗余：先定Z（不扫delayZ），再仅对该Z扫delayZ ----------
    cand_pool = setdiff(1:N, [i j]);
    step = 1;
    while ~isempty(cand_pool) && numel(sig_channels_red) < max_Z_keep
        % 阶段A：不扫 delayZ 的评分（取最小）
        deltas = zeros(1, numel(cand_pool));
        for h = 1:numel(cand_pool)
            dt    = [drivers_red, cand_pool(h)];
            TE_h  = CTE_iteration(x(:, dt), p, d, 1, 2, phase);  % 本地 1、2
            deltas(h) = mean(TE_h);
        end
        [~, pick_idx] = min(deltas);
        picked   = cand_pool(pick_idx);
        trial_dr = [drivers_red, picked];

        % 阶段B：只对 picked 扫描 delay_grid，使用“已选 delayZ 历史”
        [g_best, best_dZ] = local_best_delayZ_for_one( ...
            x, trial_dr, p, d, delay_grid, 'min', delayZ_red);

        % 置换（使用“历史 delayZ + 新 Z 的 best_dZ”）
        col = step + 1;
        if nsurr > 0
            g_surr_vec = local_surrogate_full_perZ_hist( ...
                x, p, d, trial_dr, nsurr, numel(trial_dr), [delayZ_red, best_dZ]);
            g_red_surr(:, col) = g_surr_vec;
        else
            g_surr_vec = [];
        end

        if g_best >= g_red(step), break; end

        step = step + 1;
        g_red(step) = g_best;

        if nsurr > 0
            p_values_red(step) = (sum(g_surr_vec <= g_best) + 1) / (nsurr + 1);
        else
            p_values_red(step) = NaN;
        end

        K = numel(cand_pool);
        alpha_B = alpha / max(1, K);

        if isnan(p_values_red(step)) || p_values_red(step) < alpha_B
            drivers_red = trial_dr;
            sig_mask_red(step) = true;
            sig_channels_red(end+1) = picked; %#ok<AGROW>
            delayZ_red(end+1)       = best_dZ; %#ok<AGROW>   % ★ 记录历史
            cand_pool(pick_idx) = [];
        else
            break;
        end
    end

    % ---------- 协同：先定Z（不扫delayZ），再仅对该Z扫delayZ ----------
    cand_pool = setdiff(1:N, [i j]);
    step = 1;
    while ~isempty(cand_pool) && numel(sig_channels_syn) < max_Z_keep
        % 阶段A：不扫 delayZ 的评分（取最大）
        deltas = zeros(1, numel(cand_pool));
        for h = 1:numel(cand_pool)
            dt    = [drivers_syn, cand_pool(h)];
            TE_h  = CTE_iteration(x(:, dt), p, d, 1, 2, phase);  % 本地 1、2
            deltas(h) = mean(TE_h);
        end
        [~, pick_idx] = max(deltas);
        picked   = cand_pool(pick_idx);
        trial_dr = [drivers_syn, picked];

        % 阶段B：只对 picked 扫描 delay_grid，使用“已选 delayZ 历史”
        [g_best, best_dZ] = local_best_delayZ_for_one( ...
            x, trial_dr, p, d, delay_grid, 'max', delayZ_syn);

        % 置换（使用“历史 delayZ + 新 Z 的 best_dZ”）
        col = step + 1;
        if nsurr > 0
            g_surr_vec = local_surrogate_full_perZ_hist( ...
                x, p, d, trial_dr, nsurr, numel(trial_dr), [delayZ_syn, best_dZ]);
            g_syn_surr(:, col) = g_surr_vec;
        else
            g_surr_vec = [];
        end

        if g_best <= g_syn(step), break; end

        step = step + 1;
        g_syn(step) = g_best;

        if nsurr > 0
            p_values_syn(step) = (sum(g_surr_vec >= g_best) + 1) / (nsurr + 1);
        else
            p_values_syn(step) = NaN;
        end

        K = numel(cand_pool);
        alpha_B = alpha / max(1, K);

        if isnan(p_values_syn(step)) || p_values_syn(step) < alpha_B
            drivers_syn = trial_dr;
            sig_mask_syn(step) = true;
            sig_channels_syn(end+1) = picked; %#ok<AGROW>
            delayZ_syn(end+1)       = best_dZ; %#ok<AGROW>   % ★ 记录历史
            cand_pool(pick_idx) = [];
        else
            break;
        end
    end

    % ---------- 输出整理 ----------
    sig_idx_red = find(sig_mask_red);
    sig_idx_syn = find(sig_mask_syn);

    g_red_valid = g_red(~isnan(g_red));
    g_syn_valid = g_syn(~isnan(g_syn));
    Tp     = g0;
    Tm     = min(g_red_valid);
    TM     = max(g_syn_valid);
    R_final = max(0, Tp - Tm);
    S_final = max(0, TM - Tp);
    U_final = Tm;

    sig_channels_all = unique([sig_channels_red(:); sig_channels_syn(:)])';

    last_col = max([2, 1 + numel(sig_idx_red) + numel(sig_idx_syn)]);
    g_red      = g_red(1:last_col);
    g_syn      = g_syn(1:last_col);
    p_values_red = p_values_red(1:last_col);
    p_values_syn = p_values_syn(1:last_col);
    if nsurr > 0
        g_red_surr = g_red_surr(:, 1:last_col);
        g_syn_surr = g_syn_surr(:, 1:last_col);
    end
end

function [g_opt, best_dZ] = local_best_delayZ_for_one(x, drivers_now, p, d, delay_grid, mode, delayZ_exist)
    X = x(:, drivers_now(1));   % 本地 i
    Y = x(:, drivers_now(2));   % 本地 j

    if numel(drivers_now) > 2
        Z_all = x(:, drivers_now(3:end));
        Dz0   = size(Z_all,2);
    else
        Z_all = [];
        Dz0   = 0;
    end

    % drivers_now(end) 是当前“新加入”的 Z
    Z_new = x(:, drivers_now(end));

    if nargin < 7 || isempty(delayZ_exist)
        delayZ_exist = [];           % 历史为空
    else
        delayZ_exist = delayZ_exist(:).';
    end

    if Dz0 > 1
        Z_exist = Z_all(:, 1:(end-1));    % 历史Z
    else
        Z_exist = [];
    end

    vals = zeros(1, numel(delay_grid));
    for ii = 1:numel(delay_grid)
        dZ = delay_grid(ii);

        if isempty(Z_exist)
            Zcat = Z_new;
            dvec = dZ;
        else
            Zcat = [Z_exist, Z_new];
            nZ   = size(Zcat,2);                    % 历史Z个数 + 1（新Z）
            dz   = [delayZ_exist, dZ];             % 期望与 nZ 匹配

            % —— 长度防御：不足则用1填充，过长则截断 ——
            if numel(dz) < nZ
                dz = [dz, ones(1, nZ - numel(dz))];
            elseif numel(dz) > nZ
                dz = dz(1:nZ);
            end
            dvec = dz;
        end

        vals(ii) = conditional_TE_delay(X, Y, Zcat, p, d, dvec);
    end

    if strcmpi(mode,'max')
        [g_opt, idx] = max(vals);
    else
        [g_opt, idx] = min(vals);
    end
    best_dZ = delay_grid(idx);
end


function g_surr_vec = local_surrogate_full_perZ_hist(x, p, d, drivers_now, nsurr, dr_idx, delayZ_full)
    % drivers_now: [i j Z1 ... Zk]；dr_idx 表示刚加入的是第 dr_idx 个元素
    base = x(:, drivers_now(1:dr_idx-1));
    g_surr_vec = zeros(nsurr,1);

    for k = 1:nsurr
        ys = surr_iaafft(x(:, drivers_now(dr_idx)));       % 打乱当前“新Z”
        if ~all(isfinite(ys)), ys = fillmissing(ys,'linear'); end
        xx = [base, ys];                                   % [i j Z1 ... Zk]

        if dr_idx <= 2
            Z = [];
            dvec = [];
        else
            Z = xx(:, 3:dr_idx);
            nZ = dr_idx - 2;                    % 当前 Z 的个数
            dz = delayZ_full(:).';              % 期望长度 = nZ

            % —— 长度防御：不足则用1填充，过长则截断 ——
            if numel(dz) < nZ
                dz = [dz, ones(1, nZ - numel(dz))];
            elseif numel(dz) > nZ
                dz = dz(1:nZ);
            end
            dvec = dz;
        end

        X = xx(:,1); Y = xx(:,2);
        g_surr_vec(k) = conditional_TE_delay(X, Y, Z, p, d, dvec);
    end
end

