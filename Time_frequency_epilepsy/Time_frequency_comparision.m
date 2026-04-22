clear;
clc;

Data = load('Epilepsy_alltime.mat').Epilepsy_data;

%% 参数
fs          = 400;
onset_time  = 10;              % t=10s 为发作起始
win_sec     = 3;             % 固定 3s 窗
overlap_r   = 0.5;             % 50% 重叠 → hop = 1.5s
use_pow2fft = true;

% 频段列表（行向量：fL fH），可按需增删
bands = [2 4; 4 8; 8 13; 13 30; 30 80; 80 200];

% 基线区间（相对 onset）
baseL = -10; baseH = -2;

% —— 选择归一化方式：'z' 或 'db'
norm_mode = 'z';   % 'z' = Z-score；'db' = dB 相对基线

%% 固定窗 STFT 栅格（对所有频段通用）
Ktr  = numel(Data);
[Ns, Nc] = size(Data(1).signal);

win_len = max(1, round(win_sec * fs));
noverlp = floor(overlap_r * win_len);
nfft    = use_pow2fft * 2^nextpow2(win_len) + (~use_pow2fft) * win_len;

% 参考时间/频率轴
[~, F_ref, T_ref, ~] = spectrogram(Data(1).signal(:,1), win_len, noverlp, nfft, fs);
T_rel = T_ref - onset_time;
Nt_ref = numel(T_rel);

% 基线时间索引（固定，对所有频段通用）
base_idx = (T_rel >= baseL) & (T_rel <= baseH);
if ~any(base_idx), error('基线区间为空，请检查 baseL/baseH/onset_time'); end

% 结果容器：每个频段一个 cell，元素大小 Nc × Nt_ref
nb = size(bands,1);
band_mean_all = cell(nb,1);   % 存 Z 或 dB

for b = 1:nb
    fL = bands(b,1); fH = bands(b,2);
    idxBand = (F_ref >= fL) & (F_ref <= fH);
    if ~any(idxBand)
        error('频段 [%g %g] Hz 与固定窗 STFT 频率栅格不重叠；检查 nfft/win_len', fL, fH);
    end
    df = mean(diff(F_ref));

    M = nan(Nc, Nt_ref);  % 本频段的 (通道×时间) 跨 trial 平均

    for ch = 1:Nc
        acc = zeros(1, Nt_ref);
        n_used = 0;

        for k = 1:Ktr
            x = Data(k).signal(:, ch);
            [~, ~, ~, Pk] = spectrogram(x, win_len, noverlp, nfft, fs);

            % 对齐帧数（不同 trial 长度差异时截到公共最小）
            Nt = min(size(Pk,2), Nt_ref);
            if Nt < Nt_ref
                if n_used == 0 && ch == 1
                    T_rel    = T_rel(1:Nt);
                    base_idx = base_idx(1:Nt);
                    Nt_ref   = Nt;
                    acc      = zeros(1, Nt_ref);
                    M        = nan(Nc, Nt_ref);
                end
                Pk = Pk(:, 1:Nt);
            end

            % 带内积分功率（线性域）
            band_lin = sum(Pk(idxBand,:), 1) * df;   % 1×Nt

            % 归一化：Z 或 dB（先在线性域计算基线统计）
            mu_base = mean(band_lin(base_idx));
            switch lower(norm_mode)
                case 'z'
                    sigma_base = std(band_lin(base_idx));
                    if sigma_base < 1e-12
                        y = zeros(1, Nt);
                    else
                        y = (band_lin - mu_base) / sigma_base;
                    end
                case 'db'
                    y = 10*log10( band_lin / (mu_base + eps) + eps );
                otherwise
                    error('norm_mode 仅支持 ''z'' 或 ''db''');
            end

            acc(1:Nt) = acc(1:Nt) + y(1:Nt);
            n_used = n_used + 1;
        end

        M(ch, :) = acc / max(n_used,1);
    end

    band_mean_all{b} = M;   % 存本频段的跨试次平均 (Z 或 dB)
end

% —— 可视化例子（统一色轴方便横向对比）
ez_ch = [1,2,3,4,9,10,11,12,17,18,19,26,65,66,67,68,69,70,75,76];   
figure('Color','w');
for b = 1:nb
    subplot(2,3,b);
    imagesc(T_rel, 1:Nc, band_mean_all{b}); axis xy;
    title(sprintf('PSD change in %g–%g Hz (%s) with %ds win', bands(b,1), bands(b,2), upper(norm_mode),win_sec));
    xlabel('Time from onset (s)'); 
    ylabel('Channel'); 
    hold on; 
    for c = ez_ch
        plot([T_rel(1) T_rel(end)], [c c], 'r--', 'LineWidth', 1);
    end
    plot([0 0], ylim, 'k--'); 
    hold off;
    if norm_mode == "z", caxis([-3 3]); else, % dB 按需设定，如 [-6 12]
        % caxis([-6 12]);
    end
end
colorbar('Position',[0.93 0.11 0.015 0.8]); colormap(parula);
