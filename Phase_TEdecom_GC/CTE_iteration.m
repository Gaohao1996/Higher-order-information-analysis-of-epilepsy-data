function TEv = CTE_iteration(x, p, d, idx_i, idx_j, ~, blocks, drivers_vars, dvec_cols)
% idx_i/idx_j 仍固定是 1 和 2（兼容你原来的写法）
% drivers_vars: [i j Z...]（变量编号）
% dvec_cols: Z 的“列级”延迟（若为空则无 Z 或由上层计算）

    base_vars = drivers_vars(1:2);
    Z_vars    = drivers_vars(3:end);

    Xcols = cols_of(blocks, base_vars(1));
    Ycols = cols_of(blocks, base_vars(2));

    X = x(:, Xcols);
    Y = x(:, Ycols);

    if isempty(Z_vars)
        Z = [];
        dvec = [];
    else
        Zcols = cols_of(blocks, Z_vars);
        Z = x(:, Zcols);
        if nargin < 9 || isempty(dvec_cols)
            % 如果没给 dvec，就默认每个变量 delay=1 再展开
            dvec = expand_delay_for_blocks(blocks, Z_vars, ones(1,numel(Z_vars)));
        else
            dvec = dvec_cols;
        end
    end

    TEv = conditional_TE_delay(X, Y, Z, p, d, dvec);
end




% function te = CTE_iteration(x,p,d,i,j,phase,delayZ)
%     if nargin<7 || isempty(delayZ), delayZ = 0; end
% 
%     if phase == 1
%         % 先展开，再在“列空间”里排除 i/j
%         X_all = phase2vector(x);
% 
%         i_m = linear_mapping_vector(i);   % i 的全部 cos/sin 列
%         j_m = linear_mapping_vector(j);   % j 的全部 cos/sin 列
% 
%         allcols = 1:size(X_all,2);
%         ind_m = setdiff(allcols, [i_m(:).' j_m(:).']);  % 关键修复：展开后再 setdiff
% 
%         X = X_all(:, i_m);
%         Y = X_all(:, j_m);
%         Z = X_all(:, ind_m);
%     else
%         X_all = x;
%         [~, n] = size(X_all);
%         ind = setdiff(1:n, [i j]);
% 
%         X = X_all(:, i);
%         Y = X_all(:, j);
%         Z = X_all(:, ind);
%     end
% 
%     te = conditional_TE_delay(X, Y, Z, p, d, delayZ);
% end
% 
% 
% 
% 

% 
