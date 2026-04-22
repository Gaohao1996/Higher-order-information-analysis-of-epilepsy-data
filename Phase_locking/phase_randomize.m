function data_rand = phase_randomize(data)
% data: channels x time
[ch_n, n] = size(data);
data_rand = zeros(size(data));

for ch = 1:ch_n
    X = fft(data(ch,:));
    mag = abs(X);
    rand_phase = rand(1, n) * 2 * pi - pi;
    if mod(n,2) == 0
        rand_phase(n/2+1) = 0; % Nyquist
        rand_phase(1) = 0; % DC
        rand_phase(n/2+2:end) = -fliplr(rand_phase(2:n/2));
    else
        rand_phase(1) = 0;
        rand_phase(n) = 0;
        rand_phase((n+1)/2+1:end) = -fliplr(rand_phase(2:(n-1)/2));
    end
    
    X_rand = mag .* exp(1i * rand_phase);
    data_rand(ch,:) = real(ifft(X_rand));
end
