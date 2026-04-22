function [best_L, best_delay, best_TE, best_AIC, AIC_matrix] = optimize_te_aic_delay_on_Y(x, y, z, window_lengths, delays, phase)
% 优化基于 TE 或 CTE 的阶数 L 与延迟 delay，使用 AIC 作为选择标准
%
% 输入：
%   x, y        - source 和 target（相位）变量
%   z           - 条件变量（[] 表示不使用条件变量，退化为 TE）
%   window_lengths - 模型阶数 L 的候选集合（如 1:10）
%   delays         - 延迟 d 的候选集合（如 1:20）
%
% 输出：
%   best_L, best_delay - 最优阶数与延迟
%   best_TE            - 对应传递熵值（bit）
%   best_AIC           - 对应 AIC 值
%   AIC_matrix         - 所有 (L,d) 组合对应的 AIC 值

    % Copula 映射相位数据 → 2D 向量
    if phase ==1
        X = phase2vector(x);
        Y = phase2vector(y);
    else 
        X = x;
        Y = y;
    end
        
    if ~isempty(z)
        Z = phase2vector(z);
    else
        Z = [];
    end

    best_AIC = Inf;
    best_L = NaN;
    best_delay = NaN;
    best_TE = NaN;

    num_L = length(window_lengths);
    num_d = length(delays);
    AIC_matrix = NaN(num_L, num_d);

    for j = 1:num_d
        d = delays(j);

        for i = 1:num_L
            L = window_lengths(i);

            try
                TE_bits = conditional_TE_delay(X, Y, Z, L, d, []);  % bit
            catch
                warning('TE计算出错 (L=%d, delay=%d)', L, d);
                continue;
            end

            % 无效结果跳过
            if ~isfinite(TE_bits) || TE_bits < 0
                continue;
            end

            % 有效样本数
            N_eff = size(Y, 1) - L - d;
            if N_eff <= 0, continue; end

            % log-likelihood（以 nats 为单位）
            logL = N_eff * TE_bits * log(2);  % bit → nats

            % 参数数 = 嵌入向量维度（X, Y_past, Z）
            k = L * size(X,2);  % X_past
            k = k + L * size(Y,2);  % Y_past
            if ~isempty(Z)
                k = k + L * size(Z,2);  % Z_past
            end

            [aic, ~] = aicbic(logL, k, N_eff);

            % 更新最优
            if aic < best_AIC
                best_AIC = aic;
                best_L = L;
                best_delay = d;
                best_TE = TE_bits;
            end

            AIC_matrix(i, j) = aic;
        end
    end
end
