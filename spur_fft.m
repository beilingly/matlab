clear; close all; clc;
% calc PSD with FFT 

data = textread('..\modelsim\v3p0\output_x_ofd.txt');

tcross = cumsum(data(:,1));
tcross = tcross - tcross(1);
pcoeff = polyfit(1:length(tcross), tcross, 1);
tcross_ideal = polyval(pcoeff, 1:length(tcross));
tcross_ideal = tcross_ideal';
period_osc = pcoeff(1); % average period
freq_osc = 1/period_osc;
tcross_error = tcross - tcross_ideal;
phase_error = 2*pi/period_osc * tcross_error;
phase_error = phase_error - mean(phase_error);

% fft
Nx = length(phase_error);
win = rectwin(Nx);
N = 2^nextpow2(length(phase_error));
x = phase_error.*win;
y_fft = fft(x, N);
y_power = (1/N^2) * abs(y_fft).^2;

% freq define
fres = freq_osc/N;
f = 0:fres:(N-1)*fres;
f = f';

% fig
figure(1);
plot(f/1e6, pow2db(y_power), 'Color', [1 1 1]*0.6, 'LineWidth', 2);
set(gcf, 'unit', 'normalized', 'position', [0.05,0.5,0.3,0.4]);
set(gca, 'YLim', [-100 0]);
ylabel('Power(dB)');
xlabel('Freq(MHz)');
grid minor;
hold on;

%% calc use pwelch
[freq,ps,freq_osc,jitter_total,spur_idx,spur_level] = spurpwelch(tcross, 'seg_num', 1, 'spur_fmax', 100e3);
% fig
figure(2);
plot(freq/1e6, pow2db(ps), 'Color', [1 1 1]*0.6, 'LineWidth', 2);
set(gcf, 'unit', 'normalized', 'position', [0.05,0.5,0.3,0.4]);
set(gca, 'YLim', [-100 0]);
ylabel('Power(dB)');
xlabel('Freq(MHz)');
grid minor;
hold on;