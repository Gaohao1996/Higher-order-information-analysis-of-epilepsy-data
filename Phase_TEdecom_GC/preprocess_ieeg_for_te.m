function [X_pre, report] = preprocess_ieeg_for_te(X_raw, jitter_frac)
% X_raw: samples x channels（建议先转 double）
if nargin<2, jitter_frac = 0.05; end
X_raw = double(X_raw);

% 第一次体检
S1 = tie_stats_matrix(X_raw);

% 自适应决定是否去量化
qg = S1.q_step_global;           % 例如≈0.8997803
[X_pre, L] = conditional_dequant_dither(X_raw, jitter_frac, qg);

% 第二次体检（验证效果）
S2 = tie_stats_matrix(X_pre);

% 汇总报告
report.before = S1;
report.after  = S2;
report.log    = L;
end


