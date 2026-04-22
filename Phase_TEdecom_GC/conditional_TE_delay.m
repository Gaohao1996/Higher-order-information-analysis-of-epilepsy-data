function TE_bits = conditional_TE_delay(X, Y, Z, L, delay, delayZ)
% 条件传递熵 TE = I(Y_future ; X_past | Y_past, Z_past)
% 约定：build_embedding(U, t_idx, L, delayStart) 取 U(t - (delayStart + 0..L-1))
%       即 delayStart=0 表示从 1步前开始取 L 阶（与你原先一致）

    if nargin < 5 || isempty(delay),  delay  = 1;  end   % Y 的预测步长，默认 1
    if nargin < 6, delayZ = []; end

    N = size(Y,1);
    if size(X,1) ~= N, error('X 和 Y 的样本长度不一致'); end
    if ~isempty(Z) && size(Z,1) ~= N, error('Z 和 Y 的样本长度不一致'); end

    Dx = size(X,2); Dy = size(Y,2);
    Dz = size(Z,2);

    % 统一 delayZ：标量 -> 重复到 Dz；空 -> []
    if Dz == 0
        delayZ = [];
    else
        if isempty(delayZ)
            delayZ = zeros(1, Dz);             % 默认与 Y 同步的过去（从 1步前开始）
        elseif isscalar(delayZ)
            delayZ = repmat(delayZ, 1, Dz);
        else
            delayZ = delayZ(:).';
            if numel(delayZ) ~= Dz
                error('delayZ 长度应等于 Z 的列数（%d）。', Dz);
            end
        end
    end

    % 计算最大滞后，确保所有索引有效：
    % - 对 X_past / Y_past 我们用 delayStart=0 -> 需要 L-1 的历史
    % - 对 Z_past 用 delayZ(k) 作为起点 -> 需要 delayZ(k) + (L-1) 的历史
    if Dz == 0
        maxLagZ = 0;
    else
        maxLagZ = max(delayZ + (L - 1));
    end
    maxLag = max([delay, L-1, L-1, maxLagZ]);  % Y_future 的 delay 也要计入

    % 有效时刻
    t_idx = (1 + maxLag) : (N - delay);
    Neff = numel(t_idx);
    if Neff < 1
        TE_bits = NaN;
        warning('有效样本不足（L=%d, delay=%d, maxLag=%d, N=%d）', L, delay, maxLag, N);
        return;
    end

    % 构造变量：
    % 约定：delayStart=0 -> 取 U(t-1), U(t-2), ..., U(t-L)
    Y_future = Y(t_idx + delay, :);
    X_past   = build_embedding(X, t_idx, L, 0);
    Y_past   = build_embedding(Y, t_idx, L, 0);

    if Dz == 0
        condition = Y_past;
    else
        % 注意：build_embedding 现在支持向量 delayStart（per-Z）
        Z_past = build_embedding(Z, t_idx, L, delayZ);
        condition = [Y_past, Z_past];
    end

    % 估计 TE（与你现有流程一致）
    TE_bits = gccmi_ccc(Y_future, X_past, condition);
end
