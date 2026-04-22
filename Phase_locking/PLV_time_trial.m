function PLV = PLV_time_trial(data,win_len, step, mode)
% data form: channel X samples X trials
ndat = data ./ abs(data);
[nc, ns, nt] = size(data);
win_interval = 1 : step : (ns - win_len + 1);
win_num = length(win_interval);
idx = 1;
PLV = zeros(nc,nc,win_num);

% mode of PLV calculation (1: ciPLV; 0:PLV(default))
if mode ==1
    for t = win_interval
        tmp = squeeze(ndat(:, t:t+win_len-1, :));      % ch × win_len × trial
        Z   = reshape(tmp, nc, []);                    % ch × M, M=win_len*trials
        C   = (Z * Z') / size(Z, 2);                   % 复数“相干矩阵”的PLV版本
%         ImC = imag(C);
%         % ciPLV: 逐元素归一化（避免除零）
%         denom = sqrt(max(1 - (real(C)).^2, 0));        % 或者 1 - (abs(C).^2 - ImC.^2)
%         PLV_win = ImC ./ (denom + 1e-12);
        ImC = imag(C);
        denom = sqrt(max(1 - (real(C)).^2, 0));
        PLV_win = abs(ImC) ./ (denom + 1e-12);  % 注意这里加了 abs
        PLV(:,:,idx) = PLV_win; %ciPLV
        idx = idx + 1;
    end
else  
    for s = win_interval
        tmp = squeeze(ndat(:, s:s+win_len-1, :));   % (channel × win_len × trial)
        % 合并时间维与trial维 → 相当于更多样本
        tmp2 = reshape(tmp, nc, []);
        PLV(:,:,idx) = abs(tmp2 * tmp2') / size(tmp2, 2);
        idx = idx + 1;
    end
end

end


