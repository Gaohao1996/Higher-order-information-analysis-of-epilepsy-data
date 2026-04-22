function [drivers_red, drivers_syn, g_red, g_syn, ...
    g_red_surr, g_syn_surr, ...
    sig_idx_red, sig_idx_syn, ...
    p_values_red, p_values_syn, ...
    U_final, R_final, S_final, ...
    sig_channels_red, sig_channels_syn, sig_channels_all] = ...
    TE_syn_red_greedy_surr_full(x, i, j, p, d, phase, nsurr, alpha)
% 贪婪 + IAAFT置换 + 早停 + 完整surrogate保留（即使无提升也计算一次）
% 与你的 CTE_iteration/conditional_TE_delay 兼容：
%   - i/j 在本函数内始终用“局部索引 1,2”
%   - delayZ 由下游默认值 0 处理
%
% 输入
%   x:      samples x channels （原始/相位前数据矩阵）
%   i, j:   源/目标通道索引（基于 x 的列）
%   p, d:   嵌入阶、延迟（对齐到 Y 上）
%   phase:  1=相位数据（CTE_iteration 内部相位展开）
%   nsurr:  置换次数（0 则跳过置换）
%   alpha:  显著性阈值（Bonferroni 每步/剩余候选数）
%
% 依赖：CTE_iteration / surr_iaafft / conditional_TE_delay / gccmi_ccc

% ---------- 基线 ----------
base_set = [i j];
TE0 = CTE_iteration(x(:, base_set), p, d, 1, 2, phase);  % 本地 1、2
g0  = mean(TE0);

% ---------- 预分配 ----------
[~, N] = size(x);
maxK = N;                                     % 宽松上限（每次最多加 1 个通道）
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

drivers_red = base_set;
drivers_syn = base_set;

% ---------- 冗余：取最小 ----------
cand_pool = setdiff(1:N, [i j]);
step = 1;
while ~isempty(cand_pool)
    % 扫描候选，挑更小
    deltas = zeros(1, numel(cand_pool));
    for h = 1:numel(cand_pool)
        dt    = [drivers_red, cand_pool(h)];
        TE_h  = CTE_iteration(x(:, dt), p, d, 1, 2, phase);  % 本地 1、2
        deltas(h) = mean(TE_h);
    end
    [g_best, pick_idx] = min(deltas);
    picked   = cand_pool(pick_idx);
    trial_dr = [drivers_red, picked];

    % 置换（即使无提升也计算一次，写在第 step+1 列）
    col = step + 1;
    if nsurr > 0
        g_surr_vec = local_surrogate_full(x, p, d, trial_dr, nsurr, phase, numel(trial_dr));
        g_red_surr(:, col) = g_surr_vec;
    else
        g_surr_vec = [];
    end

    % 无提升就退出（但置换已记录）
    if g_best >= g_red(step), break; end

    % 更新与显著性
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
        cand_pool(pick_idx) = [];
    else
        break;
    end
end

% ---------- 协同：取最大 ----------
cand_pool = setdiff(1:N, [i j]);
step = 1;
while ~isempty(cand_pool)
    % 扫描候选，挑更大
    deltas = zeros(1, numel(cand_pool));
    for h = 1:numel(cand_pool)
        dt    = [drivers_syn, cand_pool(h)];
        TE_h  = CTE_iteration(x(:, dt), p, d, 1, 2, phase);  % 本地 1、2
        deltas(h) = mean(TE_h);
    end
    [g_best, pick_idx] = max(deltas);
    picked   = cand_pool(pick_idx);
    trial_dr = [drivers_syn, picked];

    % 置换（即使无提升也计算一次，写在第 step+1 列）
    col = step + 1;
    if nsurr > 0
        g_surr_vec = local_surrogate_full(x, p, d, trial_dr, nsurr, phase, numel(trial_dr));
        g_syn_surr(:, col) = g_surr_vec;
    else
        g_surr_vec = [];
    end

    % 无提升就退出
    if g_best <= g_syn(step), break; end

    % 更新与显著性
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
        cand_pool(pick_idx) = [];
    else
        break;
    end
end

% ---------- 输出整理 ----------
sig_idx_red = find(sig_mask_red);
sig_idx_syn = find(sig_mask_syn);

% U/R/S（极值近似）
g_red_valid = g_red(~isnan(g_red));
g_syn_valid = g_syn(~isnan(g_syn));
Tp     = g0;
Tm     = min(g_red_valid);
TM     = max(g_syn_valid);
R_final = max(0, Tp - Tm);
S_final = max(0, TM - Tp);
U_final = Tm;

sig_channels_all = unique([sig_channels_red(:); sig_channels_syn(:)])';

% 裁剪输出到已用长度
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


function g_surr_vec = local_surrogate_full(x, p, d, drivers_now, nsurr, phase, dr_idx)
% 注意：这里不再传 i/j，全程使用本地索引 1、2
base = x(:, drivers_now(1:dr_idx-1));
g_surr_vec = zeros(nsurr,1);
for k = 1:nsurr
    ys = surr_iaafft(x(:, drivers_now(dr_idx)));
    % IAAFT 健壮性（很少见，但防御一下）
    if ~all(isfinite(ys)), ys = fillmissing(ys,'linear'); end
    xx = [base, ys];
    gc = CTE_iteration(xx, p, d, 1, 2, phase);      % 本地 [i j] = [1 2]
    g_surr_vec(k) = mean(gc);
end
end
