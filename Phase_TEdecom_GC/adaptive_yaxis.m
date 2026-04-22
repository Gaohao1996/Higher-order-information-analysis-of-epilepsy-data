function adaptive_yaxis(raw_vec, surr_mat)
% raw_vec: 1×K 或 K×1 的 TE（如 g_redF(1:last_red)）
% surr_mat: Nsurr×K 的置换矩阵（如 g_red_surrF(:, 1:last_red)）

    raw_vec = raw_vec(:);
    % 收集所有可用值
    V = [raw_vec; surr_mat(:)];
    V = V(isfinite(V));

    if isempty(V)
        ylim([0 1]); return;
    end

    % 若全非正/包含零，后续对数轴要处理
    vmin = min(V); vmax = max(V);

    % —— 阈值：量级跨度大就用 log 轴（你也可把 100 改成 50/200）
    useLog = (vmax > 0) && (vmax / max(min(V(V>0)), eps) > 100);

    if useLog
        % 确保对数轴没有 ≤0 的点：将非正值替换为极小正数用于显示（不改数据，只改轴）
        set(gca,'YScale','log');
        vpos_min = min(V(V>0));                         % 最小正值
        lo = max(vpos_min/1.5, 1e-12);                  % 下界稍放松
        hi = vmax*1.3;
        ylim([lo hi]);
        ytickformat('%.1e');
    else
        pad = 0.1*(vmax - vmin + eps);                  % 10% 边距
        lo = min(0, vmin - pad);                        % 允许含 0 以便直观看差异
        hi = vmax + pad;
        ylim([lo hi]);
        if hi>0 && hi<1e-2, ytickformat('%.1e'); end    % 数值很小时用科学计数
    end
end

