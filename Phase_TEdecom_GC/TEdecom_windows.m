close all;
clear;
clc;

%% Initialization
index = 1:7;
Trial_num = numel(index);
fs = 400;  %sample rate
%% 4.Time window shifting for dynamical analysis

% TE decomposition 
U_information =zeros (Trial_num,3);
R_information = zeros(Trial_num,3);
S_information = zeros(Trial_num,3);


%% parameter setting
%oscillation types
%  - 1: Delta: 0.5-4 Hz
%  - 2: Theta: 4-8 Hz
%  - 3: Alpha: 8-13 Hz
%  - 4: Beta: 13 - 30 Hz
%  - 5: Gamma: 30 - 80 Hz
%  - 6: HFO: 80 - 200 Hz
% Oscillation_index = 2;

%% Settting source,target and conditioned variables
% select ROI for TE
source_channel = [75 76];
condition_channel = [65 66 67 68 69 70];  %Grid_channels;%setdiff(Grid_channels,EZ_channels);
target_channel = [1,2,3,4,9,10,11,12,17,18,19,26];
%original indexes for source, target and conditioned variables
source_target_channels = [source_channel,target_channel];
Channels_idx = [source_channel, target_channel, condition_channel];

%new indexes for source, target and conditioned variables
sc_num = length(source_channel);
tc_num = length(target_channel);
cc_num = length(condition_channel);
souce_index = 1:sc_num;
target_index = sc_num+1:sc_num+tc_num;
condition_index = sc_num+tc_num+1:length(Channels_idx);


%% downsampling (prepare)
%b. mutual information calculation (optional)
% Default selection for MI (change the default option in  function named 'gccmi_ccc'  if necessary)
%biascorrect = true; % whether bias correction should be applied to the esimtated MI
%demeaned = true;  % already been copula-nomarlized so that no need to change
%cov = true; % when the covariance matrix is illconditioned use the 'false' button to reduce it (shrinkage matrix)

% TE decomposition and surrogate tests (test amount and p value setting)
nsurr=199;
alpha = 0.05;

% Phase mode (1: calculate phase TE ; 0: normal TE)
phase =0;

% Parameters optimization for TE(order and delay)
% a. order and delay
model_orders = 1:5; % not too long for test
delays = 50:100;  %based on the required cycles for oscillations

for Trial_index = index
    file_index = Trial_index;
    dataPath ='D:\Matlab_project\dataset\Epilepsy_dataset\seizure8';
    filename1 = sprintf('sz%d_ict_clean.mat', file_index );
    Data = load(fullfile(dataPath, filename1));
    eeg_data  = Data.ict_eeg;
end




%% 1.Data import
% for Trial_index = index
%     file_index = Trial_index;
%     dataPath ='D:\Matlab_project\dataset\Epilepsy_dataset\seizure8';
%     filename1 = sprintf('sz%d_ict_clean.mat', file_index );
%     Data = load(fullfile(dataPath, filename1));
%     % pre_eeg= load('sz1_pre_clean.mat').pre_eeg_clean;
%     eeg_data  = Data.ict_eeg;
%     % eeg_alltime = [pre_eeg; ict_eeg]; % [sample (M) X channels (N)]
%     
%     % Data structure
%     [M, N] =size(eeg_data);
%     Channels = 1:N;
%     Samples = 1:M;
%     
%     %% initialize matrix for TE decompostion results accross channels and trials
%     % results = zero(); %(channels X trials X windows)
%     
%     %% Parameter setting (time duration, filter, variable index, time window, parmeter optimization, surrogate tests)
%     %time duation ( epoch in [1/fs, M/fs])
%     start_sample = 1;
%     end_sample = M;
%     t_start = start_sample/fs;
%     t_end = end_sample/fs;
%     t_duration = [t_start,t_end];
%     
%     %% Zero values check in EEG data (necessary for Gaussian copula calculation)
%     % Check if the repeated 0 value exceed the threshold in data
%     % Adding jitter to the matrix for copula calculation if necessary (more zeros than expected)
%     eeg_check =  eeg_data;
%     repeat_ratio = zeros(1,N);
%     
%     for i = 1:N
%         zero_ratio = sum(eeg_check(:,i) == 0) ./ M;  % x: samples x channels
%         repeat_ratio(1,i) = zero_ratio;
%         if any(zero_ratio > 0.1)
%             warning('zero repeated ratio in this channels is too high(>0.1)');
%             % ict_eeg = jitter_zero_in_embedding(ict_eeg);
%         end
%     end
%     
%     for i = 1:M
%         time_zero_ratio = sum(eeg_check(i,:) == 0) ./N;
%         if any(time_zero_ratio > 0.1)
%             warning('zero repeated ratio at this time point is too high(>0.1)');
%             % ict_eeg = jitter_zero_in_embedding(ict_eeg);
%         end
%     end
%     
%     
%     %% 2. Filter the oscillation (optional, filter parameter need to change in data_filter if sample size changes)
% %     Filtered_EEG = data_filter(eeg_data, fs, t_duration, Channels, Oscillation_index);
%     Filtered_EEG = eeg_data;
%     %% Extract phase information from the filtered signal
%     %phase data from certain band signal
%     % Phase_EEG = angle(hilbert(Filtered_EEG)); %Theta Band Signal
%     % Phase_EEG = angle(hilbert(Filtered_EEG));
%     % x_all = Phase_EEG(:,Channels_idx);
%     x_all = Filtered_EEG(:,Channels_idx);
%     [best_p, best_d, ~, ~, ~] = optimize_te_aic_delay_on_Y(x_all(:,souce_index), x_all(:,target_index),[], model_orders, delays,phase);
%    
%     % ---------- Time shifting TE decompostion ----------
%     win_sec = 5;
%     window_length = ceil(win_sec*fs);%size(EEG_all_time,1);   %e.g. ceil(0.5*fs)
%     overlap_rate = 0.5;
%     overlap = overlap_rate*window_length; %e.g. 0.5*window_length
%     x_windows = time_window_shifting(x_all, window_length, overlap); %channels x samples x time windows
%     x_windows = permute(x_windows, [2 1 3]); % samples x channels x time windows
%     [~,~,window_number] = size(x_windows);
% 
%     
%     
%     %% 5.TE decomposition  
%     for window_index = 1: window_number
%         xw = x_windows(:,:,window_index);
%         
%         [drivers_red, drivers_syn, g_red, g_syn, ...
%             g_red_surr, g_syn_surr, ...
%             sig_idx_red, sig_idx_syn, ...
%             p_values_red, p_values_syn, ...
%             U_final, R_final, S_final, ...
%             sig_channels_red, sig_channels_syn, sig_channels_all] = ...
%             TE_syn_red_greedy_surr_full(xw, souce_index, target_index, ...
%             best_p, best_d, ...
%             phase, nsurr, alpha);
%         U_information(Trial_index,window_index) = U_final;
%         R_information(Trial_index,window_index) =R_final;
%         S_information(Trial_index,window_index) = S_final;
%     end
%     
% end
% 
% %% Box plot
% % 输入：U_information, R_information, S_information  (均为 [numTrials x numWins])
% [numTrials, numWins] = size(U_information);
% 
% % 展平
% U_vec = U_information(:);
% R_vec = R_information(:);
% S_vec = S_information(:);
% 
% % 每条数据所属的窗口编号（1..numWins）
% win_idx = repelem((1:numWins)', numTrials);
% 
% % 类型编号：1=U, 2=R, 3=S
% type_U = ones(numTrials*numWins,1);
% type_R = 2*ones(numTrials*numWins,1);
% type_S = 3*ones(numTrials*numWins,1);
% 
% % 合并
% values = [U_vec; R_vec; S_vec];
% win_all = [win_idx; win_idx; win_idx];
% type_all = [type_U; type_R; type_S];
% 
% % 为每个 (窗口, 类型) 生成唯一的组 ID：1..(numWins*3)
% group_id = (win_all-1)*3 + type_all;
% 
% % 在每个窗口位置左右错位：U at w-0.25, R at w, S at w+0.25
% offset = 0.25;
% pos = zeros(numWins*3,1);
% for w = 1:numWins
%     pos((w-1)*3+(1:3)) = [w - offset, w, w + offset];
% end
% 
% % 颜色
% colorU = [0, 0.447, 0.741];    % 蓝
% colorR = [0.85, 0.325, 0.098]; % 橙
% colorS = [0.466, 0.674, 0.188];% 绿
% colors_mat = repmat([colorU; colorR; colorS], numWins, 1);  % 为每个组着色
% 
% % 绘图
% figure('Color','w'); hold on;
% boxplot(values, group_id, ...
%     'positions', pos, ...
%     'colors', colors_mat, ...
%     'symbol', 'k+', ...       % 离群点符号
%     'whisker', 1.5, ...
%     'widths', 0.18);
% 
% % X 轴只在整数窗口上打刻度
% xlim([1-0.6, numWins+0.6]);
% xticks(1:numWins);
% xticklabels(string(1:numWins));
% xlabel('Window');
% ylabel('Information Value');
% title('U / R / S Information Distribution by Window');
% grid on;
% 
% % 图例（用虚拟线条生成）
% hU = plot(nan, nan, '-', 'Color', colorU, 'LineWidth', 6);
% hR = plot(nan, nan, '-', 'Color', colorR, 'LineWidth', 6);
% hS = plot(nan, nan, '-', 'Color', colorS, 'LineWidth', 6);
% legend([hU, hR, hS], {'U','R','S'}, 'Location','best');




