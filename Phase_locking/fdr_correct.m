function [adj_pval, sig_mask, sig_counts] = fdr_correct(pval, alpha)
% Final robust FDR correction function (3D p-value matrix)
% Inputs:
%   plv_pval: channels x channels x time_windows
%   alpha: FDR threshold (e.g., 0.05)
% Outputs:
%   adj_pval: FDR-corrected p-values (3D)
%   sig_mask: logical mask (adj_pval < alpha)
%   sig_counts: number of significant connections per window

if nargin < 2
    alpha = 0.05;
end

[ch_n, ~, win_n] = size(pval);
adj_pval = nan(ch_n, ch_n, win_n);
sig_mask = false(ch_n, ch_n, win_n);
sig_counts = zeros(1, win_n);

for w = 1:win_n
    pmat = pval(:,:,w);
    mask = triu(true(ch_n, ch_n), 1);  % 上三角不含对角线
    pvec = pmat(mask);

    if all(isnan(pvec))
        warning(['Time window ', num2str(w), ' has all NaNs. Skipping.']);
        continue;
    end

    % --- FDR correction using your 4-output fdr_bh ---
    [~, ~, ~, adj_pvec] = fdr_bh(pvec, alpha, 'pdep', 'no');

    if length(adj_pvec) ~= nnz(mask)
        warning(['Time window ', num2str(w), ': size mismatch between pvec and adj_pvec.']);
        continue;
    end

    % --- Build full symmetric matrix (zeros init to avoid NaN pollution)
    % ---%important part
    adj = zeros(ch_n, ch_n);
    adj(mask) = adj_pvec;
    adj = adj + adj';  % Mirror to lower triangle
    adj(1:ch_n+1:end) = NaN;  % Optionally mask diagonal

    % --- Store results ---
    adj_pval(:,:,w) = adj;
    sig = adj < alpha;
    sig_mask(:,:,w) = sig;
    sig_counts(w) = nnz(triu(sig, 1));  % Number of significant edges

    fprintf('FDR correction: %d / %d tests significant at q = %.3f\n', ...
        sig_counts(w), ch_n*(ch_n-1)/2, alpha);
end
end

