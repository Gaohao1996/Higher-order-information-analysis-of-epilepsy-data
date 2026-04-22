clear;
clc;

eeg_data = load('eeg_alltime.mat').eeg_data;
fs = 400;  
[M, N] =size(eeg_data);
Channels = 1:N;
Samples = 1:M;
start_sample = 1;
end_sample = M;
t_start = start_sample/fs;
t_end = end_sample/fs; 
t_duration = [t_start,t_end];

Oscillation_index = 1;
Filtered_EEG = data_filter(eeg_data, fs, t_duration, Channels, Oscillation_index); 

% time-frequency power analysis
%STFT 
window_length = ceil(3* fs);   % multiple cycles at least (Delta: 3; Theta:1.5; Alpha:0.7; Beta:0.4; Gamma: 0.2; HFO:0.07)
overlap = floor(0.5 * window_length);   % 避免非整数
nfft = 2^nextpow2(window_length);

% -------- 先跑一个通道拿帧数 --------
x0 = Filtered_EEG(:,1);
[~,F,T,P0] = spectrogram(x0, window_length, overlap, nfft, fs, 'yaxis');
num_frames = size(P0, 2);

% 预分配（按真实帧数）
total_power_band = zeros(N, num_frames);

% 频带（示例：Delta 0.5–4 Hz；根据你的滤波器来改）
fL = 0.5; 
fH = 4;

%PSD
for ch = 1:N
    signal_band = Filtered_EEG(:,ch)'; 
    [S,F,T,P] = spectrogram(signal_band, window_length, overlap, nfft, fs, 'yaxis');
     
    % 若某些通道因为长度等原因导致帧数不同，做一次对齐
    nf = size(P,2);
    if nf ~= num_frames
        % 以最小帧数为准对齐（也可选择截断或重新预分配）
        K = min(nf, num_frames);
        P = P(:, 1:K);
        if ch == 1
            total_power_band = zeros(N, K);
            num_frames = K;  % 更新全局帧数
        end
    end
    
    % 仅在带内积分（线性域）
    idx = (F >= fL) & (F <= fH);
    if ~any(idx)
        total_power_band(ch, :) = NaN;  % 带宽与频率栅格不重叠时给 NaN
        continue
    end
    df = mean(diff(F));
    band_power = sum(P(idx, :), 1) * df;   % 线性功率积分
    band_db    = 10*log10(band_power + eps);
    total_power_band(ch, :) = band_db;
end

% Threshold setting for potential SOZ channels
threshold = mean(total_power_band(:)) + 2* std(total_power_band(:));
time2exceed_threshold = ones(N,1);
for i = 1:N
    for k = 1:size(total_power_band,2)
            if total_power_band(i,k)>= threshold
                break;
            end
    end
    time2exceed_threshold(i,1) = k; 
end
[rows, cols] = find(total_power_band > threshold); 
potential_channel = unique(rows);
disp('Potential EZ channels:');
disp(potential_channel);

% plot PSD
figure;
imagesc(T, 1:N, total_power_band);
xlabel('Time (s)');
ylabel('Channel');
title('Delta Power Spectrum Change Across Channels');
colorbar;
colormap jet;
axis xy;
