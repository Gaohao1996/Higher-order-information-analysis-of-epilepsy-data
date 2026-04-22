clear;
clc;

Data = load('Epilepsy_alltime.mat').Epilepsy_data;

%% ===== 参数设置 =====
fs          = 400;      % 采样率
onset_time  = 10;       % t = 10 s 为发作起始
overlap_r   = 0.5;      % 重叠比例（0.5~0.75 常用）
use_pow2fft = true;     % if nfft= 2^nextpow2(window_length)

% 目标频段（这里以 Alpha 举例，可在外层写个循环多频段跑）
%oscillation types 
% - 1: Delta: 0.5-4 Hz 
% - 2: Theta: 4-8 Hz 
% - 3: Alpha: 8-13 Hz 
% - 4: Beta: 13 - 30 Hz 
% - 5: Gamma: 30 - 80 Hz 
% - 6: HFO: 80 - 200 Hz
bands = [0.5 4;
    4   8;
    8  13;
    13 30;
    30 80;
    80 200];
nbands = size(bands,1);
band_Z_mean_all = cell(nbands,1);
T_rel_all       = cell(nbands,1);

for b = 1:nbands
  fL = bands(b,1); 
  fH = bands(b,2);
  f0 = (fL + fH)/2;       % 中心频率

% 自适应窗：每窗 K 个周期
K_cycles   = 5;                         % 例如 5 个周期
window_sec = K_cycles / f0;            % = 5 / f0 秒
fprintf('Band [%g %g] Hz, window = %.3f s (~%d samples)\n',...
        fL,fH,window_sec,round(window_sec*fs));

% 基线区间（相对 onset 的时间）
baseL = -10;
baseH =  -2;

% 可选平滑（对最终通道×时间矩阵做一点时间滑动均值；0 表示不平滑）
smooth_win_sec = 0;

%% ===== 预备：确定 STFT 频率/时间栅格（本频段） =====
Ktr  = numel(Data);
[ Ns, Nc ] = size(Data(1).signal);

win_len = max(1, round(window_sec * fs));
noverlp = floor(overlap_r * win_len);
if use_pow2fft
    nfft = 2^nextpow2(win_len);
else
    nfft = win_len;
end

% 用第1个 trial 的第1个通道拿参考栅格
[~, F_ref, T_ref, ~] = spectrogram(Data(1).signal(:,1), win_len, noverlp, nfft, fs);
T_rel = T_ref - onset_time;          % 相对发作起始的时间轴
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

%% ===== 主计算：得到 (通道 × 时间) 的跨-trial 平均 Z 分数 =====
band_Z_mean = nan(Nc, Nt_ref);    % 输出矩阵：y=通道，x=时间（单位 Z-score）
for ch = 1:Nc
    accZ   = zeros(1, Nt_ref);    % 累加跨 trial 的 Z 轨迹
    n_used = 0;

    for k = 1:Ktr
        x = Data(k).signal(:, ch);

        % 本 trial 的 STFT
        [~, Fk, Tk, Pk] = spectrogram(x, win_len, noverlp, nfft, fs);

        % 与参考时间栅格对齐：截到公共最小帧数
        Nt = min(size(Pk,2), Nt_ref);
        if Nt < Nt_ref
            % 截断参考时间轴与基线索引（仅在第一次遇到更短 trial 时）
            if n_used == 0 && ch == 1
                T_rel    = T_rel(1:Nt);
                base_idx = base_idx(1:Nt);
                Nt_ref   = Nt;
                accZ     = zeros(1, Nt_ref);
                band_Z_mean = nan(Nc, Nt_ref);
            end
            Pk = Pk(:, 1:Nt);
        end

        % 带内“积分”功率（线性域）
        band_lin = sum(Pk(idxBand, :), 1) * df;   % 1 × Nt

        % --------- 关键改动 1：基线 Z-score，而不是 dB ---------
        base_vals = band_lin(base_idx);
        mu_base   = mean(base_vals);
        sigma_base= std(base_vals);

        if sigma_base < 1e-12
            % 避免除零；若基线极平，可直接跳过该 trial 或设为 0
            band_Z = zeros(1, Nt);
        else
            band_Z = (band_lin - mu_base) / sigma_base;   % 1 × Nt
        end
        % ------------------------------------------------------

        % 累加
        accZ(1:Nt) = accZ(1:Nt) + band_Z(1:Nt);
        n_used = n_used + 1;
    end

    band_Z_mean(ch, :) = accZ / max(n_used, 1);
end

    % 可选时间平滑（对每个通道的时间序列做 movmean）
    if smooth_win_sec > 0
        % 这里用 T_rel 长度粗略换算平滑点数
        w = max(1, round(smooth_win_sec * fs * ((Nt_ref-1) / (T_rel(end)-T_rel(1)+eps))));
        for ch = 1:Nc
            band_Z_mean(ch, :) = movmean(band_Z_mean(ch, :), w, 2);
        end
    end
    
    % ==== 存入“所有频段”的 cell 里 ====
    band_Z_mean_all{b} = band_Z_mean;   % 大小：Nc × Nt_b
    T_rel_all{b}       = T_rel;         % 对应的时间轴
end

%% interpolation and resample
% 定义一个公共时间轴，例如所有频段都统一到 0.05 s 步长
T_common = -10 : 0.05 : 10;     % 你也可以换成别的间隔

nbands = size(bands,1);
band_Z_mean_all_timealigned = cell(nbands,1);   % “4” 就理解成 final/对齐版

for b = 1:nbands
    Z  = band_Z_mean_all{b};    % Nc × Nt_b
    Tb = T_rel_all{b};          % 1 × Nt_b

    % 对每个通道插值到 T_common
    Z_aligned = interp1(Tb', Z', T_common, 'linear', 'extrap')';   % → Nc × Nt_common

    band_Z_mean_all_timealigned{b} = Z_aligned;
end

%% Visualization
ez_ch = [1,2,3,4,9,10,11,12,17,18,19,26,65,66,67,68,69,70,75,76];     
figure;
for b = 1:nbands
    subplot(2,3,b);
    imagesc(T_common, 1:Nc, band_Z_mean_all_timealigned{b});
    axis xy; caxis([-5 5]); colormap(parula);
    title(sprintf('%g–%g Hz (Z)', bands(b,1), bands(b,2)));
    hold on; 
    % 画黑色虚线在 EZ 通道上
    for c = ez_ch
        plot([T_rel(1) T_rel(end)], [c c], 'r--', 'LineWidth', 1);
    end
    plot([0 0], ylim, 'k--'); hold off;
end
colorbar('Position',[0.93 0.11 0.015 0.8]);



