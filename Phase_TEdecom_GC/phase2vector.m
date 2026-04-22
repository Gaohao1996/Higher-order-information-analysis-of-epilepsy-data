function Xgcmi = phase2vector(theta)
% PHASE2GCMI  将相位矩阵（弧度）转为 GCMI 友好的二维表示 (cos,sin)
% 输入：
%   theta : [T x C]，T 为时间点数或样本数，C 为通道数，元素为相位（弧度）
% 输出：
%   Xgcmi : [T x (2C)]，列按 [cos1 sin1 cos2 sin2 ...]
    if ~ismatrix(theta)
        error('theta 必须是 [T x C] 矩阵');
    end
    [T, C] = size(theta);
    Xgcmi = zeros(T, 2*C);
    for c = 1:C
        Xgcmi(:, 2*c-1) = cos(theta(:, c));
        Xgcmi(:, 2*c  ) = sin(theta(:, c));
    end
end

