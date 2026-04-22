function z_rep = pick_rep(z_all, mode)
% 选择代表统计量：'K1' 或 'max'
    if nargin<2 || isempty(mode), mode = 'K1'; end
    switch lower(mode)
        case 'max'
            z_rep = max(z_all);
        otherwise % 'k1'
            z_rep = z_all(1);
    end
end

