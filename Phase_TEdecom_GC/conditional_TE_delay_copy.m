function TE_bits = conditional_TE_delay(X, Y, Z, L, delay, delayZ)
    if nargin<5 || isempty(delay),  delay  = 0; end
    if nargin<6 || isempty(delayZ), delayZ = 0; end

    N = size(Y,1);

    % 同时满足 Y(t+delay) 与 Z 的滞后
    t0    = L + 1 + max(0, delayZ);
    t_idx = t0 : (N - delay);
    if numel(t_idx) < 1
        TE_bits = NaN;  warning('有效样本不足');  return;
    end

    Y_future = Y(t_idx + delay, :);
    X_past   = build_embedding(X, t_idx, L, 0);
    Y_past   = build_embedding(Y, t_idx, L, 0);

    if isempty(Z)
        condition = Y_past;
    else
        Z_past = build_embedding(Z, t_idx, L, delayZ);
        condition = [Y_past, Z_past];
    end

    % —— 稳健性护栏：仅在降秩时加极小抖动 ——
    M = [X_past condition];
    C = cov(M);
    if rank(C) < size(C,1)
        condition = condition + 1e-10*randn(size(condition));
    end

    TE_bits = gccmi_ccc(Y_future, X_past, condition);
    
%     fprintf('Neff=%d, dimX=%d, dimCond=%d, rank=%d\n', ...
%     numel(t_idx), size(X_past,2), size(condition,2), rank(cov([X_past condition])));

end
