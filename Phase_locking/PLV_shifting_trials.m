close all;
clear;
clc;
%data load
Data = load('Epilepsy_pre.mat').Preictal_data;

%% ===== 参数设置 =====
fs          = 400;      % 示例：采样率
onset_time  = 10;        %  t=10s
Nt  = numel(Data);
[Ns, Nc] = size(Data(1).signal);
Title= 1/fs:1/fs:(Ns/fs);
T_duration = [Title(1), Title(end)];
T_ref = Title-onset_time;

eeg_data = zeros(Ns,Nc,Nt);  %channels x samples x trials
%oscillation types 
%  - 1: Delta: 0.5-4 Hz
%  - 2: Theta: 4-8 Hz
%  - 3: Alpha: 8-13 Hz
%  - 4: Beta: 13 - 30 Hz
%  - 5: Gamma: 30 - 80 Hz
%  - 6: HFO: 80 - 200 Hz

Oscillation_index = 2;
for i = 1:Nt
   eeg_data(:,:,i) = data_filter(Data(i).signal, fs, T_duration, 1:Nc, Oscillation_index); 
end

analytic_eeg = zeros(Nc,Ns,Nt);
% analytical signal from eeg 
for i = 1:Nt
    analytic_eeg(:,:,i) = permute(hilbert(eeg_data(:,:,i)), [2 1 3]);
end

%Shifting windows for PLV
% δ: 5s
% α/θ（4–12 Hz）：win_len ≥ 2–4 s ；step 0.2–0.5 s；BW≈2–4 Hz;4s
% β（13–30 Hz）：win_len 0.8–2 s ；step 0.1–0.3 s；BW≈5–10 Hz;1s
% γ（30–80 Hz）：win_len 0.3–1 s；step 0.05–0.2 s；BW≈10–20 Hz
% HFO（80–250+ Hz）：win_len 50–150 ms；step 10–30 ms；BW≈40–100 Hz

%% win length (sec) setting for different oscillations (win_sec x BW x trials)
if  Oscillation_index == 1
    win_sec = 5;
    Band_name = 'Delta';
    select_windows = [1, 2, 3];
elseif Oscillation_index == 2
    win_sec = 3;
    Band_name = 'Theta';
    select_windows = [1, 2, 3,4,5];
elseif Oscillation_index == 3
    win_sec = 3;
    Band_name = 'Alpha';
    select_windows = [1, 2, 3,4,5];
else
    win_sec = 3;
    Band_name = 'Beta';
    select_windows = [1, 2, 3,4,5];
end

%% time window setting
win_len = win_sec*fs;         % samples (delta:5s; theta:4s; alpha:4s;Beta:1)
overlap_rate = 0.5;
step    = overlap_rate*win_len;         % hop

%% PLV mode\
mode = 0;
if mode == 1
    mode_name = 'ciPLV';
else
    mode_name = 'PLV';
end

plv = PLV_time_trial(analytic_eeg,win_len,step,mode); %channels x channels x window_number
% eeg_data_timewindow = time_window_shifting(eeg_data, win_len, step);
win_num = size(plv,3);

% % significance test of PLV
n_perm = 200; 
[plv_zscore, plv_pval] = plv_significance_time_trial(eeg_data, plv, n_perm,win_len,step,mode);
any_nan = sum(isnan(plv_pval(:)));

%%False Discovery Rate(FDR) test
[adj_pval, sig_mask, sig_counts] = fdr_correct(plv_pval, 0.05);

%%weighted adjacent matrix for creating functional epilepsy network
PLV_adjacent_weighted = plv.*sig_mask; 

% time window period
win_interval = 1 : step : (Ns - win_len + 1);

%%Generate the heatmap plot for PLV
for i = 1:win_num
    time_period = [win_interval(i)/fs, (win_interval(i)+win_len-1)/fs];
    figure;
    imagesc(PLV_adjacent_weighted(:, :, i));  % Heatmap plot
    set(gca, 'YDir', 'normal');              % 让Y轴从下到上递增
    axis square;
    colorbar;
    caxis([0 1]);  % Normalization
    xlabel('Channel');
    ylabel('Channel');
    title(sprintf('%s in %s band from %.2f(s) to %.2f(s) at window %d', mode_name, Band_name, time_period(1), time_period(2), i))
end

%% Select results of important period  to plot
win_num_selected = length(select_windows);
picture_index = 1;
% figure;
% for i = select_windows
    time_period = [win_interval(i)/fs, (win_interval(i)+win_len-1)/fs];
    subplot(1,win_num_selected, picture_index)
        imagesc(PLV_adjacent_weighted(:, :, i));  % Heatmap plot
    set(gca, 'YDir', 'normal');              % 让Y轴从下到上递增
    colorbar;
    caxis([0 1]);  % Normalization
    xlabel('Channel');
    ylabel('Channel');
    title(sprintf('%s in %s band from %.2f(s) to %.2f(s) at window %d',mode_name,Band_name, time_period(1), time_period(2), i))
    picture_index = picture_index+1;
% end

figure('Position',[100 100 1800 600]);   % 加大 figure 尺寸
pic = tiledlayout(1,win_num_selected,'Padding','compact','TileSpacing','compact');

for i = 1:win_num_selected
    time_period = [win_interval(select_windows(i))/fs, (win_interval(select_windows(i))+win_len-1)/fs];
    nexttile
    imagesc(PLV_adjacent_weighted(:, :, select_windows(i)));  % Heatmap plot
    set(gca, 'YDir', 'normal');              % 让Y轴从下到上递增
    caxis([0 1]);  % Normalization
    xlabel('Channel');
    ylabel('Channel');
    title(sprintf('From %.2f(s) to %.2f(s) at window %d', time_period(1), time_period(2), select_windows(i)))
    axis square
end
cb = colorbar;
cb.Layout.Tile = 'east';
annotation('textbox',[0 0.9 1 0.06], ...   % [x y w h]
    'String', sprintf('%s in %s band', mode_name, Band_name), ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'EdgeColor','none','FontSize',15,'FontWeight','normal');



