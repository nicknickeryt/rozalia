clc; clear; 

dane = readRigolbin("RigolDS2.bin");

Fs = 625e6;          % częstotliwość próbkowania (625 MSa/s)
N = length(dane);    % liczba próbek

w = hamming(N);
Y = fft(dane .* w');

% poprawna normalizacja amplitudy
A = abs(Y) / sum(w);

% poprawka na jednostronne widmo
A(2:end-1) = 2*A(2:end-1);

% dBV (1 V RMS odniesienie)
Y_dBV = 20*log10(A / sqrt(2));


R = 50; % albo 51

V_rms = A / sqrt(2);
P = (V_rms.^2) / R;

Y_dBm = 10*log10(P / 1e-3);

% oś częstotliwości
f = (0:N-1)*(Fs/N);

% tylko dodatnie częstotliwości
half = 1:floor(N/2);
f = f(half);
Y = Y(half);

% normalizacja amplitudy
Y = abs(Y)/N;

% zakres FM
fmin = 87e6;
fmax = 108.5e6;

idx = (f >= fmin) & (f <= fmax);

f_plot = f(idx);
Y_plot = Y_dBm(idx);

% --- DETEKCJA PIKÓW ---
minHeight = max(Y_plot) - 50;   % próg (np. 40 dB poniżej maksimum)
minDist = 200e3;               % min odstęp między pikami (100 kHz)

% zamiana na indeksy (bo findpeaks działa na próbkach)
df = f_plot(2) - f_plot(1);
minDist_pts = round(minDist / df);

[pks, locs] = findpeaks(Y_plot, ...
    'MinPeakHeight', minHeight, ...
    'MinPeakDistance', minDist_pts);

f_peaks = f_plot(locs);

% --- WYKRES ---
figure;
plot(f_plot/1e6, Y_plot); hold on;
plot(f_peaks/1e6, pks, 'ro');

xlabel('Częstotliwość [MHz]');
ylabel('Amplituda [dBm]');
title('FFT w paśmie FM + piki');
grid on;

% --- OPISY NAD PIKAMI ---
for i = 1:length(pks)
    text(f_peaks(i)/1e6, pks(i), ...
        sprintf('%.2f MHz\n [%.2f MHz]\n%.1f dBm', f_peaks(i)/1e6 + 10.7, f_peaks(i)/1e6, pks(i)), ...
        'VerticalAlignment','bottom', ...
        'HorizontalAlignment','center', ...
        'FontSize',8);
end