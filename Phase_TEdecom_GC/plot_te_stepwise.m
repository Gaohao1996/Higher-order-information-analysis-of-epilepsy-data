function plot_te_stepwise(raw_vec, surr_mat, titleStr)
% raw_vec: 1×K 或 K×1（每步加入一个条件后的 TE）
% surr_mat: Nsurr×K（第1列通常是 NaN；从第2列才有值）

    raw_vec = raw_vec(:)';                     % 行向量
    K = numel(raw_vec);

    % 置换均值（从第2列开始），并让它对齐在对应的 x=k
    cols = find(any(isfinite(surr_mat),1));    % 有效的列
    mu_surr = nan(1,K);                        % 先全 NaN
    if ~isempty(cols)
        mu_surr(cols) = mean(surr_mat(:,cols), 1, 'omitnan');
    end

    % --- 开始画 ---
    figure; hold on; box on
    x = 1:K;

    % Raw（保留 NaN 为缺测）
    plot(x, raw_vec, '-*k', 'LineWidth', 1.8, 'DisplayName','Raw');

    % Surrogate（与 k 对齐；只画非 NaN）
    idx = isfinite(mu_surr);
    plot(x(idx), mu_surr(idx), 'o', 'MarkerSize', 6, ...
        'DisplayName','Surrogate mean');

    % 对齐 y 轴范围（自动线性/对数）
    vals = [raw_vec(:); mu_surr(:)];
    vals = vals(isfinite(vals));
    if isempty(vals)
        ylim([0 1]);  % 兜底
    else
        vmin = min(vals); vmax = max(vals);
        if vmax <= 0           % 全非正（极少见），给个兜底范围
            ylim([vmin*1.2, 1]);
        else
            if vmin <= 0, vmin = min(vals(vals>0)); end  % 对数轴不能 ≤0
            if vmax/vmin > 100   % 跨 2 个数量级以上就用对数轴
                set(gca, 'YScale', 'log');
                ylim([max(vmin*0.8, eps), vmax*1.2]);
            else
                ylim([min(0, vmin*0.9), vmax*1.2]);
            end
        end
    end
    xlim([0.8, K+0.2])
    xlabel('Number of conditioned variables (k)')
    ylabel('TE')
    title(titleStr)
    legend('Location','best')

    % 明确标注“无置换”/“无TE”的步
    kk_no_surr = setdiff(1:K, cols);
    if ~isempty(kk_no_surr)
        for kk = kk_no_surr
            text(kk, ylim(gca)*[0;1]*0 + min(ylim), '\downarrow no surr', ...
                'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
                'FontSize',8, 'Color',[.4 .4 .4])
        end
    end
    kk_no_raw = find(~isfinite(raw_vec));
    if ~isempty(kk_no_raw)
        plot(kk_no_raw, repmat(min(ylim),size(kk_no_raw)), 'x', ...
            'MarkerSize', 8, 'LineWidth', 1.2, 'DisplayName','no TE')
    end
end
