function E = build_embedding(U, t_idx, L, delay)
% 构造过去嵌入：
%  - U: [N x D]
%  - t_idx: 有效时刻索引
%  - L: 阶数（取 delay, delay+1, ..., delay+L-1）
%  - delay: 可为标量（所有列同起始滞后）或 1xD 向量（每列各自起始滞后）
%
% 生成：
%   E = [ U_d(t - (delay_d + 0)), U_d(t - (delay_d + 1)), ..., U_d(t - (delay_d + L - 1)) ]_d  拼接

    if isvector(U), U = U(:); end
    [N, D] = size(U);

    % 统一成行向量
    if isscalar(delay)
        delay = repmat(delay, 1, D);
    else
        delay = delay(:).';
        if numel(delay) ~= D
            error('build_embedding: delay 向量长度应等于 U 的列数（%d）。', D);
        end
    end

    % 预分配
    E = zeros(numel(t_idx), D*L);

    col = 1;
    for d = 1:D
        d0 = delay(d);
        for ell = 0:(L-1)
            lag = d0 + ell;                     % 本列的具体滞后
            E(:, col) = U(t_idx - lag, d);
            col = col + 1;
        end
    end
end

