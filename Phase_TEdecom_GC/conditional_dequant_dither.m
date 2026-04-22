function [X2, loginfo] = conditional_dequant_dither(X, jitter_frac, q_step_global)
% X: samples x channels (double)
% jitter_frac: 抖动幅度占量化步长比例（默认 0.05 = 5%）
% q_step_global: 可留空，自动从每列估计；如你已测得≈0.8997803可在此传入
if nargin<2 || isempty(jitter_frac), jitter_frac = 0.05; end
X = double(X);
[N,D] = size(X);
X2 = X;
loginfo.dithered = false(1,D);
loginfo.sigma    = zeros(1,D);
loginfo.q_step   = nan(1,D);

for d = 1:D
    xd = X(:,d);
    ufrac  = numel(unique(xd))/N;
    repadj = mean([0; diff(xd)]==0);
    if (ufrac<0.9) || (repadj>0.1)
        % 估计该列的 q_step（优先用全局估计）
        if nargin>=3 && ~isempty(q_step_global)
            qstep = q_step_global;
        else
            u = unique(xd);
            du = diff(u); du = du(du>0);
            if ~isempty(du), qstep = median(du); else, qstep = NaN; end
        end
        if ~isfinite(qstep) || qstep<=0
            qstep = 1e-8*std(xd); % 兜底
        end
        sig = jitter_frac * qstep;
        X2(:,d) = xd + sig*randn(N,1);
        loginfo.dithered(d) = true;
        loginfo.sigma(d)    = sig;
        loginfo.q_step(d)   = qstep;
    end
end
end


