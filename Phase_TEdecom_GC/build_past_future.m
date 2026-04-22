function [Y_future, Y_past, X_past] = build_past_future(Y, X, L, delay)
% BUILD_PAST_FUTURE  基于窗口阶数 L 和延迟 delay 构造 TE 的三个矩阵
% 约定：
%   对每个 t ∈ {L+delay+1, ..., N}:
%     Y_future(t) = Y(t)
%     Y_past  包含 [Y(t-1), Y(t-2), ..., Y(t-L)]
%     X_past  包含 [X(t-delay), X(t-delay-1), ..., X(t-delay-L+1)]
%
% 输入：
%   Y : [N x Dy]   目标多维序列
%   X : [N x Dx]   源多维序列
%   L : 标量，嵌入阶数
%   delay : 标量，滞后（通常 = 1）
%
% 输出：
%   Y_future : [Neff x Dy]
%   Y_past   : [Neff x (L*Dy)]
%   X_past   : [Neff x (L*Dx)]

    if size(X,1) ~= size(Y,1)
        error('X 与 Y 行数（样本数）必须相同');
    end
    N  = size(Y,1);
    Dy = size(Y,2);
    Dx = size(X,2);

    Neff = N - L - delay;
    if Neff < 1
        error('有效样本数不足：N-L-delay 必须 ≥ 1');
    end

    t_idx   = (L + delay + 1 : N)';        % 未来时刻索引
    Y_future = Y(t_idx, :);

    Y_past = zeros(Neff, L*Dy);
    X_past = zeros(Neff, L*Dx);

    % Y_past: t-1 ... t-L
    for l = 1:L
        Y_past(:, (l-1)*Dy + (1:Dy)) = Y(t_idx - l, :);
    end
    % X_past: t-delay ... t-delay-L+1
    for l = 1:L
        X_past(:, (l-1)*Dx + (1:Dx)) = X(t_idx - delay - l + 1, :);
    end
end

