function cols = cols_of(blocks, var_ids)
%COLS_OF  把“变量编号列表”映射成列号列表（支持相位2列、幅度1列）
% blocks{k} = 该变量在 x 中占用的列索引向量
    if isempty(var_ids), cols = []; return; end
    cols = cell2mat(blocks(var_ids));
end


