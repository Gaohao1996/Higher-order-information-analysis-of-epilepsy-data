close all;
clear;
clc;

%% Initialization
%% 1.Data import
file_index = 2;
dataPath ='D:\Matlab_project\dataset\Epilepsy_dataset\seizure8';
filename1 = sprintf('sz%d_ict_clean.mat', file_index );
Data = load(fullfile(dataPath, filename1));
% pre_eeg= load('sz1_pre_clean.mat').pre_eeg_clean;
eeg_data  = Data.ict_eeg;
% eeg_alltime = [pre_eeg; ict_eeg]; % [sample (M) X channels (N)]

% Data structure
fs = 400;  %sample rate
[M, N] =size(eeg_data);
Channels = 1:N;
Samples = 1:M;

%% initialize matrix for TE decompostion results accross channels and trials
% results = zero(); %(channels X trials X windows)

%% Parameter setting (time duration, filter, variable index, time window, parmeter optimization, surrogate tests)
%time duation ( epoch in [1/fs, M/fs])
start_sample = 1;
end_sample = M;
t_start = start_sample/fs;
t_end = end_sample/fs;
t_duration = [t_start,t_end];

%% fix seed
% rng(1,'twister');

%% Zero values check in EEG data (necessary for Gaussian copula calculation)
%Due to finite ADC resolution, the broadband iEEG exhibited discrete amplitude levels
% (median quantization step ≈ 0.90 in acquisition units), 
% which led to a high proportion of identical samples after delay embedding. 
% We applied minimal de-quantization by adding Gaussian noise with σ = 5
% of the estimated quantization step prior to copula normalization, to ensure valid Gaussian-copula CMI estimation.”
[eeg_data, rep] = preprocess_ieeg_for_te(eeg_data, 0.05);

% （可选）打印一下效果
fprintf('Before: median uniq=%.3f, frac_bad=%.2f, q_step≈%.6f\n', ...
    median(rep.before.uniq_frac), rep.before.frac_bad, rep.before.q_step_global);
fprintf('After : median uniq=%.3f, frac_bad=%.2f\n', ...
    median(rep.after.uniq_frac), rep.after.frac_bad);


%% oscillation types
%  - 1: Delta: 0.5-4 Hz
%  - 2: Theta: 4-8 Hz
%  - 3: Alpha: 8-13 Hz
%  - 4: Beta: 13 - 30 Hz
%  - 5: Gamma: 30 - 80 Hz
%  - 6: HFO: 80 - 200 Hz

Oscillation_index = 2;

%% Settting source,target and conditioned variables
% select ROI for TE
source_channel = [75 76];
condition_channel = [65 66 67 68 69 70];  %Grid_channels;%setdiff(Grid_channels,EZ_channels);
target_channel = [1,2,3,4,9,10,11,12,17,18,19,26];
%original indexes for source, target and conditioned variables
source_target_channels = [source_channel,target_channel];
Channels_idx = [source_channel, target_channel, condition_channel];

%new indexes for source, target and conditioned variables
sc_num = length(source_channel);
tc_num = length(target_channel);
cc_num = length(condition_channel);
souce_index = 1:sc_num;
target_index = sc_num+1:sc_num+tc_num;
condition_index = sc_num+tc_num+1:length(Channels_idx);

%% downsampling (if need)

%b. mutual information calculation (optional)
% Default selection for MI (change the default option in  function named 'gccmi_ccc'  if necessary)
%biascorrect = true; % whether bias correction should be applied to the esimtated MI
%demeaned = true;  % already been copula-nomarlized so that no need to change
%cov = true; % when the covariance matrix is illconditioned use the 'false' button to reduce it (shrinkage matrix)

% TE decomposition and surrogate tests (test amount and p value setting)
nsurr=199;
alpha = 0.05;

%% Phase mode (1: phase TE ; 0: normal TE)
% phase =1;

% Parameters optimization for TE(order and delay)
% a. order and delay
model_orders = 2:6; % not too long for test
delays = 8:60;  %based on the required cycles for oscillations

% 2. Filter the oscillation (optional, filter parameter need to change in data_filter if sample size changes)
Filtered_EEG = data_filter(eeg_data, fs, t_duration, Channels, Oscillation_index);

%% Extract phase information from the filtered signal
%phase data from certain band signal
% Phase_EEG = angle(hilbert(Filtered_EEG)); %Theta Band Signal
% x_all = Phase_EEG(:,Channels_idx);
% x_all = Filtered_EEG(:,Channels_idx);

%% variable mapping
% 1) amp-amp
X = []; blocks = {};
for ch = Channels_idx
    a = Filtered_EEG(:,ch);
%     a = abs(hilbert(Filtered_EEG(:,ch))); % 或直接用 Filtered(:,ch)
    blocks{end+1} = size(X,2)+1;
    X = [X, a];
end
x_all = X;

% 2) phase-phase
% phi = angle(hilbert(Filtered_EEG));
% %构造 x 与 blocks（每个变量2列：sinφ, cosφ）
% X = []; blocks = {};
% for ch = Channels_idx
%     s = sin(phi(:, ch));
%     c = cos(phi(:, ch));
%     blocks{end+1} = (size(X,2) + (1:2));
%     X = [X, s, c]; %#ok<AGROW>
% end
% x_all = X;

%3)phase-amp
% H   = hilbert(Filtered_EEG);
% phi = angle(H);
% amp = abs(H);
% 
% X = []; blocks = {};
% for idx = 1:numel(Channels_idx)
%     ch = Channels_idx(idx);
%     if ismember(idx, souce_index)         % 源相位→2列
%         blocks{end+1} = (size(X,2)+(1:2));
%         X = [X, sin(phi(:,ch)), cos(phi(:,ch))];
%     elseif ismember(idx, target_index)    % 目标幅度→1列
%         blocks{end+1} = size(X,2)+1;
%         X = [X, amp(:,ch)];
%     else                                   % 条件（举例：也用相位）
%         blocks{end+1} = (size(X,2)+(1:2));
%         X = [X, sin(phi(:,ch)), cos(phi(:,ch))];
%     end
% end
% x_all = X;



%% ----------TE decompostion ----------
nsurr_small = 80;         % 全时段“摸底”小置换
cand_T = [4 5 6 8 10];    % 候选窗口秒数
overlap_rate = 0.5;       % 滑动重叠比例
rep_mode = 'K1';          % 代表统计量：'K1' 或 'max'（K=1 或 max_k）

% [best_p0, best_d0, ~, ~, ~] = optimize_te_aic_delay_on_Y(x_all(:,souce_index), x_all(:,target_index), [], model_orders, delays, phase);

[best_p0, best_d0, ~, ~, ~] = optimize_te_aic_delay_on_Y( ...
    x_all(:, cols_of(blocks, souce_index)), ...
    x_all(:, cols_of(blocks, target_index)), [], 2:6, 8:60);

f_theta = 6;
delay_grid = unique(max(round(fs * (0.25:0.25:1.5) / f_theta), 1));

% 全时段（正式 nsurr）
% 3) 你可以把 max_Z_keep 设成条件通道数量（或你想要的上限）
max_Z_keep = numel(condition_index);   % 或者 3、2 等

% 4) 调用
[drivers_redF, drivers_synF, g_redF, g_synF, ...
 g_red_surrF, g_syn_surrF, sig_idx_redF, sig_idx_synF, ...
 p_values_redF, p_values_synF, U, R, S, ...
 sig_channels_redF, sig_channels_synF, ~, ...
 delayZ_redF, delayZ_synF] = TE_syn_red_greedy_surr_full(x_all, souce_index, target_index, ...
        best_p0, best_d0, nsurr, alpha, delay_grid, max_Z_keep, blocks);
 
%     [drivers_red, drivers_syn, g_red,...
%   g_syn, g_red_surr, g_syn_surr, ...sig_idx_red, sig_idx_syn, ...
%  p_values_red, p_values_syn, U, R, S, ...
%   sig_channels_red, sig_channels_syn, ~, delayZ_red, delayZ_syn] = ...
%     TE_syn_red_greedy_surr_full(x_all, souce_index, target_index, ...
%         best_p0, best_d0, nsurr, alpha, delay_grid, max_Z_keep, blocks);
% [drivers_redF, drivers_synF, g_redF, g_synF, ...
%  g_red_surrF, g_syn_surrF, sig_idx_redF, sig_idx_synF, ...
%  p_values_redF, p_values_synF, U, R, S, ...
%  sig_channels_redF, sig_channels_synF, ~, ...
%  delayZ_redF, delayZ_synF] = ...
%     TE_syn_red_greedy_surr_full( ...
%         x_all, souce_index, target_index, best_p0, best_d0, ...
%         phase, nsurr, alpha, delay_grid, max_Z_keep);

% === 这里定义 full_sig_set（冗余/协同显著通道并集）===
full_sig_set = unique([sig_channels_redF(:); sig_channels_synF(:)])';
% 如果你希望用“真实通道号”而不是列索引用于一致性比较：
% full_sig_set = unique( Channels_idx(full_sig_set) )';

% [drivers_red0, drivers_syn0, ~, ~, ~, ~, sig_idx_red0, sig_idx_syn0, ~, ~, ~, ~, ~, ...
%     sig_channels_red0, sig_channels_syn0, full_sig_set0] = ...
%     TE_syn_red_greedy_surr_full(x_all, souce_index, target_index, best_p0, best_d0, phase, nsurr_small, alpha);


% %% ===== 通用：取本次实际使用到的列数 =====
last_red = find(~isnan(g_redF), 1, 'last');      % g_redF/g_synF 是 1×K
last_syn = find(~isnan(g_synF), 1, 'last');
if isempty(last_red), last_red = numel(g_redF); end
if isempty(last_syn), last_syn = numel(g_synF); end

% plot_te_stepwise(g_synF, g_syn_surrF, 'Synergistic information flow');
% plot_te_stepwise(g_redF, g_red_surrF, 'Redundant information flow');


%% ----- 冗余（无连线的置换点） -----
figure; hold on;
x_red = 1:last_red;

% Raw
plot(x_red, g_redF(1:last_red), '-*k', 'LineWidth', 2);

% Surrogate 点（从第2列）
scols_red = intersect(2:last_red, 1:size(g_red_surrF,2));
if ~isempty(scols_red)
    Ypts  = g_red_surrF(:, scols_red);                          % nsurr × K
    Xpts  = repmat(scols_red, size(Ypts,1), 1);                 % nsurr × K
    hSurr = scatter(Xpts(:), Ypts(:), 8, 'r', 'filled', ...     % 不连线
                    'MarkerFaceAlpha', 0.15, 'MarkerEdgeAlpha', 0.15);
end

% —— 冗余显著点标注 ——
for k = 1:numel(sig_idx_redF)
    t  = sig_idx_redF(k);
    ch = Channels_idx(sig_channels_redF(k));
    text(t, g_redF(t), sprintf('%d', ch), ...
         'VerticalAlignment','top', 'HorizontalAlignment','center', ...
         'FontSize',10, 'Color','g', 'FontWeight','bold');
end
plot(x_red(sig_idx_redF), g_redF(sig_idx_redF), ...
     'p', ...                         % 五角星 / 也可用 'o' 或 '*'
     'MarkerSize', 10, ...
     'MarkerEdgeColor', 'g', ...
     'MarkerFaceColor', 'g', ...
     'LineStyle', 'none');
xlim([0, max(3, last_red)]);
xlabel('Number of Conditioned Variables');
ylabel(sprintf('TE values with order %d and lag %d', best_p0, best_d0));
title('Redundant information flow test at time window');
legend('Raw signal','Surrogate Signal','Location','best'); grid on;
adaptive_yaxis(g_redF(1:last_red), g_red_surrF(:, 1:last_red));

%% ----- 协同（无连线的置换点） -----
figure; hold on;
x_syn = 1:last_syn;

% Raw
plot(x_syn, g_synF(1:last_syn), '-*k', 'LineWidth', 2);

% Surrogate 点（从第2列）
scols_syn = intersect(2:last_syn, 1:size(g_syn_surrF,2));
if ~isempty(scols_syn)
    Ypts  = g_syn_surrF(:, scols_syn);
    Xpts  = repmat(scols_syn, size(Ypts,1), 1);
    scatter(Xpts(:), Ypts(:), 8, 'r', 'filled', 'MarkerFaceAlpha', 0.15, 'MarkerEdgeAlpha', 0.15);
end

% —— 协同显著点标注 ——
for k = 1:numel(sig_idx_synF)
    t  = sig_idx_synF(k);
    ch = Channels_idx(sig_channels_synF(k));
    text(t, g_synF(t), sprintf('%d', ch), ...
         'VerticalAlignment','top', 'HorizontalAlignment','center', ...
         'FontSize',10, 'Color','g', 'FontWeight','bold');
end
plot(x_syn(sig_idx_synF), g_synF(sig_idx_synF), ...
     'p', ...                         % 五角星 / 也可用 'o' 或 '*'
     'MarkerSize', 10, ...
     'MarkerEdgeColor', 'g', ...
     'MarkerFaceColor', 'g', ...
     'LineStyle', 'none');
xlim([0, max(3, last_syn)]);
xlabel('Number of Conditioned Variables');
ylabel(sprintf('TE values with order %d and lag %d', best_p0, best_d0));
title('Synergistic information flow test at time window');
legend('Raw signal','Surrogate Signal','Location','best'); grid on;
adaptive_yaxis(g_synF(1:last_syn), g_syn_surrF(:, 1:last_syn));



%% ====== 2) window length estimation固定 p*, d*） ======
% summary_tbl = [];   % [WinSec, NumWins, MeanConsistency, FracWinsWithSig, MeanZ_red, MeanZ_syn]

% for win_sec = cand_T
%     window_length = ceil(win_sec * fs);
%     overlap = round(overlap_rate * window_length);
%     
%     % 切窗：channels x samples x wins -> samples x channels x wins
%     x_windows = time_window_shifting(x_all, window_length, overlap);
%     x_windows = permute(x_windows, [2 1 3]);
%     WN = size(x_windows, 3);
%     
%     % 每窗记录：显著集合 + z*
%     max_keep = numel(condition_index);                          % 每窗最多存多少显著通道编号（可调）
%     win_sig = zeros(WN, max_keep);
%     z_red = nan(WN,1); z_syn = nan(WN,1);
%     has_sig = false(WN,1);
%     
%     for w = 1:WN
%         xw = x_windows(:,:,w);
%         % 注意参数顺序：phase, nsurr, alpha
%         [drivers_red, drivers_syn, g_red, g_syn, ...
%             g_red_surr, g_syn_surr, sig_idx_red, sig_idx_syn, ~, ~, ...
%             ~, ~, ~, sig_channels_red, sig_channels_syn] = ...
%             TE_syn_red_greedy_surr_full(xw, souce_index, target_index, ...
%             best_p0, best_d0, ...
%             phase, nsurr, alpha);
%         
%         % === 显著集合（直接用返回的显著通道，合并去重） ===
%         sig_set_w = unique([sig_channels_red(:); sig_channels_syn(:)])';
%         sig_set_w = sig_set_w(sig_set_w > 0);       % 保险：过滤占位0
%         has_sig(w) = ~isempty(sig_set_w);
%         
%         if ~isempty(sig_set_w)
%             ncap = min(numel(sig_set_w), max_keep);
%             win_sig(w,1:ncap) = sig_set_w(1:ncap);
%             % 多余位置保持为 0 作为占位；compute_consistency 里请忽略0
%         end
%         
%         % ===== 安全工具 =====
%         % clamp p 到 (1/(nsurr+1), 1-1/(nsurr+1))，避免 norminv 出 NaN/Inf
%         clamp_p = @(p) max(min(p, 1 - 1/(nsurr+1)), 1/(nsurr+1));
%         
%         % 优先用 condition_index(t-1) 取列；若该列不存在/全NaN，再用 t+1
%         pick_col = @(S, t) ...
%             ( exist('condition_index','var') && numel(condition_index) >= (t-1) && ...
%             any(~isnan(S(:, condition_index(t-1)))) ) .* condition_index(t-1) + ...
%             ~( exist('condition_index','var') && numel(condition_index) >= (t-1) && ...
%             any(~isnan(S(:, condition_index(t-1)))) ) .* (t+1);
%         
%         aggfun = @(v) max(v,[],'omitnan');   % 或 median
%         
%         % ====== 冗余（只看被接受步；若没有，就兜底看已计算步） ======
%         last_red = find(~isnan(g_red), 1, 'last');
%         if ~isempty(sig_idx_red)
%             zlist = nan(1, numel(sig_idx_red));
%             for ii = 1:numel(sig_idx_red)
%                 t = sig_idx_red(ii);                % t=2,3,...
%                 c = pick_col(g_red_surr, t);
%                 if c <= size(g_red_surr,2) && ~isnan(g_red(t))
%                     s = g_red_surr(:, c); s = s(~isnan(s));
%                     if ~isempty(s)
%                         p = (sum(s <= g_red(t)) + 1) / (numel(s) + 1);
%                         zlist(ii) = -norminv(clamp_p(p));   % 越大越显著
%                     end
%                 end
%             end
%             z_red(w) = aggfun(zlist);
%         elseif ~isempty(last_red) && last_red >= 2
%             % 没有“被接受步”时的兜底（可删掉以实现“严格只看被接受步”）
%             K = last_red;
%             zlist = nan(1, K);
%             for t = 2:K
%                 c = pick_col(g_red_surr, t);
%                 if c <= size(g_red_surr,2) && ~isnan(g_red(t))
%                     s = g_red_surr(:, c); s = s(~isnan(s));
%                     if ~isempty(s)
%                         p = (sum(s <= g_red(t)) + 1) / (numel(s) + 1);
%                         zlist(t) = -norminv(clamp_p(p));
%                     end
%                 end
%             end
%             z_red(w) = aggfun(zlist);
%         else
%             z_red(w) = NaN;
%         end
%         
%         % ====== 协同 ======
%         last_syn = find(~isnan(g_syn), 1, 'last');
%         if ~isempty(sig_idx_syn)
%             zlist = nan(1, numel(sig_idx_syn));
%             for ii = 1:numel(sig_idx_syn)
%                 t = sig_idx_syn(ii);
%                 c = pick_col(g_syn_surr, t);
%                 if c <= size(g_syn_surr,2) && ~isnan(g_syn(t))
%                     s = g_syn_surr(:, c); s = s(~isnan(s));
%                     if ~isempty(s)
%                         p = (sum(s >= g_syn(t)) + 1) / (numel(s) + 1);
%                         zlist(ii) = -norminv(clamp_p(p));
%                     end
%                 end
%             end
%             z_syn(w) = aggfun(zlist);
%         elseif ~isempty(last_syn) && last_syn >= 2
%             K = last_syn;
%             zlist = nan(1, K);
%             for t = 2:K
%                 c = pick_col(g_syn_surr, t);
%                 if c <= size(g_syn_surr,2) && ~isnan(g_syn(t))
%                     s = g_syn_surr(:, c); s = s(~isnan(s));
%                     if ~isempty(s)
%                         p = (sum(s >= g_syn(t)) + 1) / (numel(s) + 1);
%                         zlist(t) = -norminv(clamp_p(p));
%                     end
%                 end
%             end
%             z_syn(w) = aggfun(zlist);
%         else
%             z_syn(w) = NaN;
%         end
%         
%         
%     end
%     
%     % 与全时段集合比的一致性（逐窗 recall）
%     %  确保 compute_consistency 内部忽略 win_sig==0 的占位元素
%     stats = compute_consistency(win_sig, full_sig_set, 2/3);
%     mean_consist = stats.mean_consistency;
%     frac_win_sig = mean(has_sig);
%     
%     summary_tbl = [summary_tbl; [win_sec, WN, mean_consist, frac_win_sig, ...
%         mean(z_red,'omitnan'), mean(z_syn,'omitnan')]];
% end
% 
% disp(array2table(summary_tbl, ...
%     'VariableNames', {'WinSec','NumWins','MeanConsistency','FracWinsWithSig','MeanZ_red','MeanZ_syn'}));





%% ====== 3) 可选：把最佳窗口（例如 MeanConsistency最高者）用于后续画图/分析 ======
% [~, idx_best] = max(summary_tbl(:,3)); best_win = summary_tbl(idx_best,1);
% 之后你就用 best_win 重新切窗，做结果展示即可。








