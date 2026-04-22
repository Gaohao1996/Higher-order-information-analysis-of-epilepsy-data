function S = tie_stats_matrix(X)
% X: samples x channels (double)
X = double(X);
[N,D] = size(X);
S.N = N; S.D = D;
S.uniq_frac = zeros(1,D);
S.repeat_adj = zeros(1,D);
S.q_step = nan(1,D);
for d = 1:D
    xd = X(:,d);
    u  = unique(xd);
    S.uniq_frac(d)  = numel(u)/N;            % 唯一值比例
    S.repeat_adj(d) = mean([0; diff(xd)]==0);% 相邻样本完全相等比例
    if numel(u)>2
        du = diff(u); du = du(du>0);
        if ~isempty(du), S.q_step(d) = median(du); end
    end
end
S.frac_bad = mean(S.uniq_frac<0.9 | S.repeat_adj>0.1); % 不达标列占比
S.q_step_global = median(S.q_step(~isnan(S.q_step)));   % 全通道统一步长估计
end

