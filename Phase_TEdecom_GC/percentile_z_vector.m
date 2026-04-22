function z_all = percentile_z_vector(raw_vec, surr_mat)
% raw_vec: 1 x K    （每个K步的TE原始值）
% surr_mat: nsurr x K（对应置换分布）
% 输出：z_all: 1 x K 的百分位z（z* = norminv(p))
    if isrow(raw_vec), raw_vec = raw_vec(:).'; end
    [ns, K] = size(surr_mat);
    z_all = nan(1,K);
    for k = 1:K
        r = raw_vec(k);
        s = surr_mat(:,k);
        p = (1 + sum(s <= r)) / (ns + 1);   % 一侧（raw > surr）
        p = min(max(p, 1e-6), 1-1e-6);      % 数值裁剪，避免 inf
        z_all(k) = norminv(p);
    end
end

