function [plv_zscore, plv_pval] = plv_significance_time_trial(data, plv_real, n_perm, win_len, step, mode)
% data: samples x channels x trials
% fs: sampling rate (e.g., 400)
% win_size: window length in seconds (e.g., 1.5)
% step_size: step size in seconds (e.g., 0.1)
% n_perm: number of permutations (e.g., 200)
[s_num, ch_num, trial_num] = size(data);
data = permute(data, [2 1 3]); % turn to the shape of channels x samples x trials (for randomizing phase)
win_interval = 1 : step : (s_num - win_len + 1);
win_num = length(win_interval);
plv_zscore = zeros(ch_num, ch_num, win_num);
plv_pval = zeros(ch_num, ch_num, win_num);

shuffled_data = zeros(ch_num, win_len);
shuffled_analytic_data = zeros(ch_num,win_len,trial_num);

for w = 1:win_num
    % 构造null分布
    plv_null = zeros(ch_num, ch_num, n_perm);
    for p = 1:n_perm
        for Nt = 1:trial_num
            %randomize phase within each time window in each trail
            shuffled_data = phase_randomize(data(:,win_interval(w):win_interval(w)+win_len-1,Nt))';  % samples x channels
%         shuffled_data = permute(shuffled_data, [2 1 3]); % samples (win) x channels x trials
            shuffled_analytic_data(:,:,Nt) = hilbert(shuffled_data)'; % channels x samples
        end
        plv_null(:,:,p) = PLV_time_trial(shuffled_analytic_data,win_len, win_len, mode); % calculate the null hyposis for PLV in a window   
    end

    % 4. 计算z-score和p值
    plv_mean = mean(plv_null, 3);
    plv_std = std(plv_null, [], 3);
    plv_zscore(:,:,w) = (plv_real(:,:,w) - plv_mean) ./ plv_std;

    % p值：右尾检验
    for i = 1:ch_num-1
        for j = i+1:ch_num
            real_val = plv_real(i,j,w);
            null_vals = squeeze(plv_null(i,j,:));
            plv_pval(i,j,w) = mean(real_val <= null_vals); % = sum(null_vals >= real_val) / n_perm;
            plv_pval(j,i,w) = plv_pval(i,j,w);
        end
    end
end

end




