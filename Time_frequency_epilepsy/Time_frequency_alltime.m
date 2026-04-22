clear;
clc;

Data = load('Epilepsy_alltime.mat').Epilepsy_data;

%% ===== 参数设置 =====
fs          = 400;      % 示例：采样率
onset_time  = 10;        %  t=10s
overlap_r   = 0.5;       % 重叠比例（0.5~0.75 常用）
use_pow2fft = true;      
Oscillation_index = 2;
bands = [0.5 4;
    4   8;
    8  13;
    13 30;
    30 80;
    80 200];

%% Bands range
fL = bands(Oscillation_index,1);
fH = bands(Oscillation_index,2);


%% multiple cycles at least (Delta: 3; Theta:1.5; Alpha:0.7; Beta:0.4; Gamma: 0.2; HFO:0.07)
if Oscillation_index == 1
    window_sec = 3;
elseif Oscillation_index == 2
    window_sec = 1.5;
elseif Oscillation_index == 3
    window_sec = 0.7;
elseif Oscillation_index == 4
    window_sec = 0.4;
elseif Oscillation_index == 5
    window_sec = 0.2;
else
    window_sec = 0.07;
end
    
% 目标频段（示例：ripple 80–250 Hz；请按需修改）
%oscillation types 
%  - 1: Delta: 0.5-4 Hz
%  - 2: Theta: 4-8 Hz
%  - 3: Alpha: 8-13 Hz
%  - 4: Beta: 13 - 30 Hz
%  - 5: Gamma: 30 - 80 Hz
%  - 6: HFO: 80 - 200 Hz

% 基线区间（相对 onset 的时间）
baseL = -10; 
baseH =  -2;

% 可选平滑（对最终通道×时间矩阵做一点时间滑动均值；0 表示不平滑）
smooth_win_sec = 0;  % 例如 0.1（100ms）；0 = 不平滑


%% ===== 预备：确定 STFT 频率/时间栅格 =====
K  = numel(Data);
[ Ns, Nc ] = size(Data(1).signal);
win_len = max(1, round(window_sec * fs));
overlp = floor(overlap_r * win_len);
% if use_pow2fft
%     nfft = 2^nextpow2(win_len);
% else
%     nfft = win_len;
% end
nfft    = use_pow2fft * 2^nextpow2(win_len) + (~use_pow2fft) * win_len;
% 用第1个 trial 的第1个通道拿参考栅格
[~, F_ref, T_ref, P0] = spectrogram(Data(1).signal(:,1), win_len, overlp, nfft, fs);
T_rel = T_ref - onset_time;                 % 相对发作起始的时间轴
Nt_ref = numel(T_rel);

% 频带索引与 df
idxBand = (F_ref >= fL) & (F_ref <= fH);
df      = mean(diff(F_ref));
if ~any(idxBand)
    error('频段 [%g, %g] Hz 与 STFT 频率栅格不重叠，请检查窗长/nfft。', fL, fH);
end

% 基线帧索引
base_idx = (T_rel >= baseL) & (T_rel <= baseH);
if ~any(base_idx)
    error('基线区间 [%g, %g] s 在时间轴上为空，请检查 onset_time / baseL / baseH。', baseL, baseH);
end


%% ===== 主计算：得到 (通道 × 时间) 的跨-trial 平均 dB 变化 =====
band_dB_mean = nan(Nc, Nt_ref);        % 输出矩阵：y=通道，x=时间（单位 dB 相对基线）

for ch = 1:Nc
    acc = zeros(1, Nt_ref);            % 累加跨 trial 的 dB 轨迹
    n_used = 0;

    for k = 1:K
        x = Data(k).signal(:, ch);

        % 本 trial 的 STFT
        [~, Fk, Tk, Pk] = spectrogram(x, win_len, overlp, nfft, fs);

        % 与参考时间栅格对齐：截到公共最小帧数
        Nt = min(size(Pk,2), Nt_ref);
        if Nt < Nt_ref
            % 截断参考时间轴与基线索引（仅在第一次遇到更短 trial 时）
            if n_used == 0 && ch == 1
                T_rel   = T_rel(1:Nt);
                base_idx = base_idx(1:Nt);
                Nt_ref   = Nt;
                acc      = zeros(1, Nt_ref);
                band_dB_mean = nan(Nc, Nt_ref);
            end
            Pk = Pk(:, 1:Nt);
        end

        % 带内“积分”功率（线性域）
        band_lin = sum(Pk(idxBand, :), 1) * df;   % 1 × Nt

        % 基线均值（线性域）
        mu_base = mean(band_lin(base_idx));

        % 相对基线 dB 变化
        band_rel_dB = 10*log10(band_lin / (mu_base + eps) + eps);

        % 累加
        acc(1:Nt) = acc(1:Nt) + band_rel_dB(1:Nt);
        n_used = n_used + 1;
    end

    band_dB_mean(ch, :) = acc / max(n_used, 1);
end

% 可选时间平滑（对每个通道的时间序列做 movmean）
if smooth_win_sec > 0
    w = max(1, round(smooth_win_sec * fs * ( (Nt_ref-1) / (T_rel(end)-T_rel(1)+eps) ))); % 近似换算
    for ch = 1:Nc
        band_dB_mean(ch, :) = movmean(band_dB_mean(ch, :), w, 2);
    end
end


%% ===== 绘图：通道 × 时间 的热图 =====
ez_ch = [1,2,3,4,9,10,11,12,17,18,19,26,65,66,67,68,69,70,75,76];     % 你的 EZ 通道编号（按原索引）
figure; imagesc(T_rel, 1:Nc, band_dB_mean); axis xy;
colormap(parula); colorbar; xlabel('Time from onset (s)'); ylabel('Channel');

% 画黑色虚线在 EZ 通道上
hold on;
for c = ez_ch
    plot([T_rel(1) T_rel(end)], [c c], 'r--', 'LineWidth', 1);
end
plot([0 0], ylim, 'k--','LineWidth', 2); % onset
hold off;



% figure('Color','w');
% imagesc(T_rel, 1:Nc, band_dB_mean); axis xy;
% xlabel('Time from onset (s)');
% ylabel('Channel');
% title(sprintf('Band power change (dB rel. baseline), %g–%g Hz', fL, fH));
% hcb = colorbar; ylabel(hcb, 'dB');
% % 可选：固定色轴，便于不同病人/通道对比
% % caxis([-3 3]);  % 例如固定在 ±3 dB
% colormap(parula);
% 
% % 可选：标注发作起点
% hold on; plot([0 0], ylim, 'k--', 'LineWidth', 1); hold off;



