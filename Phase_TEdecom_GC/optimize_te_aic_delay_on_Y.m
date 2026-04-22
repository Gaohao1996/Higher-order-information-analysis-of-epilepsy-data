function [best_L, best_delay, best_TE, best_AIC, AIC_matrix] = ...
    optimize_te_aic_delay_on_Y(x, y, z, window_lengths, delays)
% AIC-based optimization of (L, delay) for TE / CTE
% Assumes x,y,z are already properly embedded (e.g. phase mapped to sin/cos)

    X = x;
    Y = y;
    Z = z;

    best_AIC = Inf;
    best_L = NaN;
    best_delay = NaN;
    best_TE = NaN;

    num_L = numel(window_lengths);
    num_d = numel(delays);
    AIC_matrix = NaN(num_L, num_d);

    for jd = 1:num_d
        d = delays(jd);

        for il = 1:num_L
            L = window_lengths(il);

            try
                TE_bits = conditional_TE_delay(X, Y, Z, L, d, []);
            catch
                continue;
            end

            if ~isfinite(TE_bits) || TE_bits < 0
                continue;
            end

            % 有效样本数
            N_eff = size(Y,1) - L - d;
            if N_eff <= 0, continue; end

            % log-likelihood (nats)
            logL = N_eff * TE_bits * log(2);

            % 参数数（自由度）
            k = L * size(X,2);   % X_past
            k = k + L * size(Y,2); % Y_past
            if ~isempty(Z)
                k = k + L * size(Z,2); % Z_past
            end

            [aic, ~] = aicbic(logL, k, N_eff);
            AIC_matrix(il, jd) = aic;

            if aic < best_AIC
                best_AIC = aic;
                best_L = L;
                best_delay = d;
                best_TE = TE_bits;
            end
        end
    end
end
