function stats = compute_consistency(win_sig, full_sig, presence_thresh)
% win_sig: K x N，逐窗显著通道编号（0代表无）；full_sig: 1 x M
    if nargin<3, presence_thresh = 2/3; end
    K = size(win_sig,1);
    Sfull = unique(full_sig(full_sig>0)); nfull = numel(Sfull);
    consistency = nan(K,1); precision = nan(K,1); f1 = nan(K,1);
    for w = 1:K
        Sw = unique(win_sig(w, win_sig(w,:)>0));
        tp = numel(intersect(Sw, Sfull));
        if nfull>0, consistency(w) = tp / nfull; end
        if ~isempty(Sw), precision(w) = tp / numel(Sw); end
        if ~isnan(consistency(w)) && ~isnan(precision(w)) && (consistency(w)+precision(w)>0)
            f1(w) = 2*consistency(w)*precision(w)/(consistency(w)+precision(w));
        end
    end
    stable_rate = []; stable_frac = NaN;
    if nfull>0
        present_mat = false(K, nfull);
        for j = 1:nfull
            c = Sfull(j);
            present_mat(:,j) = any(win_sig == c, 2);
        end
        stable_rate = mean(present_mat, 1);
        stable_frac = mean(stable_rate >= presence_thresh);
    end
    stats.consistency_per_window = consistency;
    stats.precision_per_window   = precision;
    stats.f1_per_window          = f1;
    stats.mean_consistency       = mean(consistency,'omitnan');
    stats.mean_precision         = mean(precision,'omitnan');
    stats.mean_f1                = mean(f1,'omitnan');
    stats.full_set               = Sfull;
    stats.stable_rate_per_channel= stable_rate;
    stats.stable_frac_ge_thresh  = stable_frac;
end