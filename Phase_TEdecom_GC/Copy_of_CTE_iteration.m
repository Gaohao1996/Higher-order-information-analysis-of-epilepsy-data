function te = CTE_iteration(x,p,d,i,j,phase,delayZ)
    if nargin<7 || isempty(delayZ), delayZ = 0; end

    if phase == 1
        % 先展开，再在“列空间”里排除 i/j
        X_all = phase2vector(x);

        i_m = linear_mapping_vector(i);   % i 的全部 cos/sin 列
        j_m = linear_mapping_vector(j);   % j 的全部 cos/sin 列

        allcols = 1:size(X_all,2);
        ind_m = setdiff(allcols, [i_m(:).' j_m(:).']);  % 关键修复：展开后再 setdiff

        X = X_all(:, i_m);
        Y = X_all(:, j_m);
        Z = X_all(:, ind_m);
    else
        X_all = x;
        [~, n] = size(X_all);
        ind = setdiff(1:n, [i j]);

        X = X_all(:, i);
        Y = X_all(:, j);
        Z = X_all(:, ind);
    end

    te = conditional_TE_delay(X, Y, Z, p, d, delayZ);
end




% function te=CTE_iteration(x,p,d,i,j,phase)
% % x data S (samples) x n (variables) dati 
% % i driver j target
% % p order of the model
% [~, n]=size(x);
% ind=setdiff(1:n,[i j]);
% %% X,Y and Z variables
% if phase == 1
%     X_all = phase2vector(x);
%     i = linear_mapping_vector(i);
%     j = linear_mapping_vector(j);
%     ind = linear_mapping_vector(ind);
% else
%      X_all =x;
% end
%     
% X = X_all(:,i);
% Y = X_all(:,j);
% Z = X_all(:,ind);
% te= conditional_TE_delay(X,Y,Z,p,d);


