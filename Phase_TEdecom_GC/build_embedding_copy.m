function emb = build_embedding(data, t_idx, L, delay)
% BUILD_EMBEDDING 构造过去嵌入向量
%
%   emb = build_embedding(data, t_idx, L, delay)
%
% 输入：
%   data   - [N x D] 原始时间序列（每列一个变量）
%   t_idx  - [Neff x 1] 时间索引（如 L+1:N-delay）
%   L      - 嵌入阶数（过去时刻数量）
%   delay  - 延迟，从 t-delay 开始嵌入 L 个点
%
% 输出：
%   emb    - [Neff x (L*D)] 嵌入向量，每行是 t-delay, ..., t-delay-L+1 的拼接

    [~, D] = size(data);
    Neff = length(t_idx);

    % 检查是否越界
    min_required_index = min(t_idx - delay - (L - 1));
    if min_required_index < 1
        error('build_embedding: 嵌入越界，最早需要的数据索引为 %d', min_required_index);
    end

    % 初始化输出
    emb = zeros(Neff, L * D);

    % 构建嵌入向量：每列对应一个延迟的变量
    for l = 1:L
        idx_lag = t_idx - delay - (l - 1);
        emb(:, (l - 1)*D + (1:D)) = data(idx_lag, :);
    end
end
