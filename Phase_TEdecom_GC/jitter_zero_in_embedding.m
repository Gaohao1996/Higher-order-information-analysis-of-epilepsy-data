function X_jittered = jitter_zero_in_embedding(X, eps_scale)
% JITTER_ZERO_IN_EMBEDDING - add jitter to 0 element in the matrix
%   X_jittered = jitter_zero_in_embedding(X, eps_scale)
% 输入：
%   X         - 输入矩阵，样本 x 特征（可为嵌入后的数据）
%   eps_scale - jitter 幅度（默认 1e-6），表示扰动标准差
% 输出：
%   X_jittered - 添加 jitter 后的矩阵，仅 0 值被扰动

    if nargin < 2
        eps_scale = 1e-12;
    end

    noise = eps_scale * randn(size(X));

    zero_mask = (X == 0);
    X_jittered = X;
    X_jittered(zero_mask) = X_jittered(zero_mask) + noise(zero_mask);

%     n_zeros = sum(zero_mask(:));
%     if n_zeros > 0
%         fprintf('[jitter_zero_in_embedding] 添加 jitter 到 %d 个零值点。\n', n_zeros);
%     else
%         fprintf('[jitter_zero_in_embedding] 没有检测到零值，无需处理。\n');
%     end
end

