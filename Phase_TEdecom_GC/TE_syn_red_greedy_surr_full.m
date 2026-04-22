function [drivers_red, drivers_syn, g_red, g_syn, ...
    g_red_surr, g_syn_surr, ...
    sig_idx_red, sig_idx_syn, ...
    p_values_red, p_values_syn, ...
    U_final, R_final, S_final, ...
    sig_channels_red, sig_channels_syn, sig_channels_all, ...
    delayZ_red, delayZ_syn] = ...
    TE_syn_red_greedy_surr_full(x, i, j, p, d, nsurr, alpha, delay_grid, max_Z_keep, blocks)
% TE_syn_red_greedy_surr_full (简化版，直接调用 conditional_TE_delay)
% - 所有“相位2列/幅度1列”的映射在函数外完成，这里只按 blocks 切片。
% - 冗余：单尾左尾；协同：单尾右尾。
% - delay_grid：按“变量”为单位，内部自动复制到该变量的所有列。
%
% 输入
%   x           : samples × dims（所有变量列拼接）
%   i, j        : 源/目标 “变量编号”（与 blocks 的顺序对齐）
%   p, d        : 模型阶数/延迟（X→Y）
%   nsurr       : 置换次数（0 表示不做置换）
%   alpha       : 显著性水平
%   delay_grid  : 仅对“当前加入的 Z 变量”的候选 delay（按变量）
%   max_Z_keep  : 最多保留多少个条件变量（上限）
%   blocks      : 1×M cell，blocks{k} = 该“变量”的列索引向量
%
% 输出（与原版一致）
%   drivers_red/syn : 变量编号列表（[i j Z1 ... Zk]）
%   g_red/g_syn     : 每步的 TE 值（k=1 是基线 TE）
%   g_red_surr/g_syn_surr : nsurr×K 的置换分布（第1列为空）
%   p_values_*      : 每步 p 值（第1列=1）
%   U/R/S           : 信息分解三项（基于 g_red/g_syn）
%   sig_channels_*  : 被接受进入模型的 Z 变量编号
%   delayZ_*        : 对应每个已接受 Z 的“每变量 delay”（与 sig_channels_* 对齐）

    if nargin < 9 || isempty(delay_grid), delay_grid = 1; end
    if nargin < 10 || isempty(max_Z_keep), max_Z_keep = inf; end
    if nargin < 11 || isempty(blocks)
        % 默认：每个变量占 1 列
        [~,D] = size(x);
        blocks = arrayfun(@(c) c, 1:D, 'uni', 0);
    end

    % ---------- 基线 ----------
    base_vars = [i j];
    g0 = mean(eval_TE_for_blocks(x, p, d, blocks, base_vars, [], []));

    % ---------- 预分配 ----------
    M = numel(blocks);                % 变量总数
    maxK = M;
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

    % 历史：每个已选 Z 的“每变量 delay”
    delayZ_red = [];
    delayZ_syn = [];

    drivers_red = base_vars;
    drivers_syn = base_vars;

    % ========== 冗余 ==========
    cand_pool = setdiff(1:M, base_vars);
    step = 1;
    while ~isempty(cand_pool) && numel(sig_channels_red) < max_Z_keep
        % A. 不扫 delay 的评分：更小更好
        deltas = zeros(1, numel(cand_pool));
        for h = 1:numel(cand_pool)
            dt_vars = [drivers_red, cand_pool(h)];
            TE_h = eval_TE_for_blocks(x, p, d, blocks, dt_vars(1:2), dt_vars(3:end), []);
            deltas(h) = mean(TE_h);
        end
        [~, pick_idx] = min(deltas);
        picked = cand_pool(pick_idx);

        % B. 仅对 picked 扫 delay_grid
        [g_best, best_dZ] = local_best_delayZ_for_one_var( ...
            x, p, d, blocks, drivers_red, picked, delay_grid, 'min', delayZ_red);

        % 置换：历史 delayZ + 新变量的 best_dZ
        col = step + 1;
        if nsurr > 0
            g_surr_vec = local_surrogate_full_perZ_hist_blocks( ...
                x, p, d, blocks, [drivers_red, picked], nsurr, ...
                [delayZ_red, best_dZ]);
            g_red_surr(:, col) = g_surr_vec;
        else
            g_surr_vec = [];
        end

        if g_best >= g_red(step), break; end

        step = step + 1;
        g_red(step) = g_best;

        if nsurr > 0
            % 左尾（更小更显著）
            p_values_red(step) = (sum(g_surr_vec <= g_best) + 1) / (nsurr + 1);
        else
            p_values_red(step) = NaN;
        end

        % Bonferroni
        K = numel(cand_pool);
        alpha_B = alpha / max(1, K);

        if isnan(p_values_red(step)) || p_values_red(step) < alpha_B
            drivers_red = [drivers_red, picked];
            sig_mask_red(step) = true;
            sig_channels_red(end+1) = picked; %#ok<AGROW>
            delayZ_red(end+1)       = best_dZ; %#ok<AGROW>
            cand_pool(pick_idx)     = [];
        else
            break;
        end
    end

    % ========== 协同 ==========
    cand_pool = setdiff(1:M, base_vars);
    step = 1;
    while ~isempty(cand_pool) && numel(sig_channels_syn) < max_Z_keep
        % A. 不扫 delay 的评分：更大更好
        deltas = zeros(1, numel(cand_pool));
        for h = 1:numel(cand_pool)
            dt_vars = [drivers_syn, cand_pool(h)];
            TE_h = eval_TE_for_blocks(x, p, d, blocks, dt_vars(1:2), dt_vars(3:end), []);
            deltas(h) = mean(TE_h);
        end
        [~, pick_idx] = max(deltas);
        picked = cand_pool(pick_idx);

        % B. 仅对 picked 扫 delay_grid
        [g_best, best_dZ] = local_best_delayZ_for_one_var( ...
            x, p, d, blocks, drivers_syn, picked, delay_grid, 'max', delayZ_syn);

        % 置换
        col = step + 1;
        if nsurr > 0
            g_surr_vec = local_surrogate_full_perZ_hist_blocks( ...
                x, p, d, blocks, [drivers_syn, picked], nsurr, ...
                [delayZ_syn, best_dZ]);
            g_syn_surr(:, col) = g_surr_vec;
        else
            g_surr_vec = [];
        end

        if g_best <= g_syn(step), break; end

        step = step + 1;
        g_syn(step) = g_best;

        if nsurr > 0
            % 右尾（更大更显著）
            p_values_syn(step) = (sum(g_surr_vec >= g_best) + 1) / (nsurr + 1);
        else
            p_values_syn(step) = NaN;
        end

        % Bonferroni
        K = numel(cand_pool);
        alpha_B = alpha / max(1, K);

        if isnan(p_values_syn(step)) || p_values_syn(step) < alpha_B
            drivers_syn = [drivers_syn, picked];
            sig_mask_syn(step) = true;
            sig_channels_syn(end+1) = picked; %#ok<AGROW>
            delayZ_syn(end+1)       = best_dZ; %#ok<AGROW>
            cand_pool(pick_idx)     = [];
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
    g_red        = g_red(1:last_col);
    g_syn        = g_syn(1:last_col);
    p_values_red = p_values_red(1:last_col);
    p_values_syn = p_values_syn(1:last_col);
    if nsurr > 0
        g_red_surr = g_red_surr(:, 1:last_col);
        g_syn_surr = g_syn_surr(:, 1:last_col);
    end
end

% ======================= 本文件内局部工具 =======================

% function dvec = expand_delay_for_blocks(blocks, Z_ids, delay_per_var)
%     if isempty(Z_ids), dvec = []; return; end
%     assert(numel(delay_per_var) == numel(Z_ids), 'delay/Z 数不匹配');
%     dvec = [];
%     for k = 1:numel(Z_ids)
%         ncol = numel(blocks{Z_ids(k)});
%         dvec = [dvec, repmat(delay_per_var(k), 1, ncol)];
%     end
% end
function dvec = expand_delay_for_blocks(blocks, Z_ids, delay_per_var)
% 将“每变量 delay”扩展到“每列 delay”
% —— 自动防御 delay/Z 数不匹配（greedy + surrogate 必需）

    if isempty(Z_ids)
        dvec = [];
        return;
    end

    nZ = numel(Z_ids);

    % ===== 长度防御 =====
    if isempty(delay_per_var)
        delay_per_var = ones(1, nZ);
    elseif numel(delay_per_var) < nZ
        delay_per_var = [delay_per_var, ones(1, nZ - numel(delay_per_var))];
    elseif numel(delay_per_var) > nZ
        delay_per_var = delay_per_var(1:nZ);
    end

    % ===== 展开到列 =====
    dvec = [];
    for k = 1:nZ
        ncol = numel(blocks{Z_ids(k)});
        dvec = [dvec, repmat(delay_per_var(k), 1, ncol)];
    end
end


function TEv = eval_TE_for_blocks(x, p, d, blocks, base_vars, Z_vars, dvec_per_col)
    X = x(:, cols_of(blocks, base_vars(1)));
    Y = x(:, cols_of(blocks, base_vars(2)));

    if isempty(Z_vars)
        Z = []; dvec = [];
    else
        Z = x(:, cols_of(blocks, Z_vars));
        if isempty(dvec_per_col)
            dvec = expand_delay_for_blocks(blocks, Z_vars, ones(1,numel(Z_vars)));
        else
            dvec = dvec_per_col;
        end
    end
    TEv = conditional_TE_delay(X, Y, Z, p, d, dvec);
end

function [g_opt, best_dZ] = local_best_delayZ_for_one_var( ...
        x, p, d, blocks, drivers_now_vars, new_var_id, delay_grid, mode, delayZ_hist)

    if nargin < 9 || isempty(delayZ_hist), delayZ_hist = []; end

    base_vars  = drivers_now_vars(1:2);
    Z_hist_ids = drivers_now_vars(3:end);   % 已选历史Z（变量编号）

    % —— 对历史 delay 向量做长度防护，使其与 Z_hist_ids 数量匹配 ——
    if numel(delayZ_hist) < numel(Z_hist_ids)
        delayZ_hist = [delayZ_hist, ones(1, numel(Z_hist_ids) - numel(delayZ_hist))];
    elseif numel(delayZ_hist) > numel(Z_hist_ids)
        delayZ_hist = delayZ_hist(1:numel(Z_hist_ids));
    end

    vals = zeros(1, numel(delay_grid));
    for ii = 1:numel(delay_grid)
        dZ = delay_grid(ii);                                % 本次仅对“新变量”扫描的 delay
        delay_per_var = [delayZ_hist, dZ];                  % [历史Z的delay, 新Z的delay]
        dvec = expand_delay_for_blocks(blocks, [Z_hist_ids, new_var_id], delay_per_var);

        % ★ 不再传 phase；eval_TE_for_blocks 只有 7 个参数
        TEv  = eval_TE_for_blocks(x, p, d, blocks, base_vars, [Z_hist_ids, new_var_id], dvec);
        vals(ii) = mean(TEv);
    end

    if strcmpi(mode,'max'), [g_opt, idx] = max(vals);
    else,                  [g_opt, idx] = min(vals);
    end
    best_dZ = delay_grid(idx);
end




% function [g_opt, best_dZ] = local_best_delayZ_for_one_var( ...
%         x, p, d, blocks, drivers_now_vars, new_var_id, delay_grid, mode, delayZ_hist)
% 
%     if nargin < 9 || isempty(delayZ_hist), delayZ_hist = []; end
%     base_vars = drivers_now_vars(1:2);
%     Z_hist_ids = drivers_now_vars(3:end);
% 
%     vals = zeros(1, numel(delay_grid));
%     for ii = 1:numel(delay_grid)
%         dZ = delay_grid(ii);
%         delay_per_var = [delayZ_hist, dZ];   % 每变量 delay
%         dvec = expand_delay_for_blocks(blocks, [Z_hist_ids, new_var_id], delay_per_var);
%         TEv  = eval_TE_for_blocks(x, p, d, blocks, base_vars, [Z_hist_ids, new_var_id], dvec);
%         vals(ii) = mean(TEv);
%     end
% 
%     if strcmpi(mode,'max'), [g_opt, idx] = max(vals);
%     else,                  [g_opt, idx] = min(vals);
%     end
%     best_dZ = delay_grid(idx);
% end

function g_surr_vec = local_surrogate_full_perZ_hist_blocks( ...
        x, p, d, blocks, drivers_now_vars, nsurr, delayZ_per_var)
% 仅随机化“新加入”的那个变量块（其所有列），其他不动

    Z_vars    = drivers_now_vars(3:end);
    dvec = expand_delay_for_blocks(blocks, Z_vars, delayZ_per_var);

    base_vars = drivers_now_vars(1:2);
    baseX = cols_of(blocks, base_vars(1));
    baseY = cols_of(blocks, base_vars(2));
    Zcols  = cols_of(blocks, Z_vars);

    new_var  = Z_vars(end);
    new_cols = cols_of(blocks, new_var);

    g_surr_vec = zeros(nsurr,1);
    for k = 1:nsurr
        xx = x;

        % IAAFT：块内每列分别 IAAFT（你也可换成统一随机循环移位）
        for cc = 1:numel(new_cols)
            ys = surr_iaafft(xx(:, new_cols(cc)));
            if ~all(isfinite(ys)), ys = fillmissing(ys,'linear'); end
            xx(:, new_cols(cc)) = ys;
        end
        % 统一循环移位版本（可选）：
        % lag = randi(size(x,1)); xx(:, new_cols) = circshift(xx(:, new_cols), lag);

        X = xx(:, baseX);
        Y = xx(:, baseY);
        if isempty(Zcols), Z = []; dvec_col = [];
        else,              Z = xx(:, Zcols); dvec_col = dvec;
        end
        g_surr_vec(k) = conditional_TE_delay(X, Y, Z, p, d, dvec_col);
    end
end
