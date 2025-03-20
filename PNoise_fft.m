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
win = hann(length(phase_error));
N = 2^nextpow2(length(phase_error));
x = phase_error.*win;
y_fft = fft(phase_error, N);
y_power = (1/N^2) * abs(y_fft).^2;
enbw = freq_osc * sum(abs(win).^2)/abs(sum(win))^2;
y_psd = y_power / enbw;
y_psd_ssb = y_psd(1:N/2+1);
y_psd_ssb(2:end-1) = 2*y_psd_ssb(2:end-1);
PNxx = pow2db(y_psd_ssb);

% freq define
fres = freq_osc/N;
f = 0:fres:(N-1)*fres;
f = f(1:N/2+1);
fxx = f';

figure(1);
semilogx(fxx, PNxx, 'r','LineWidth',2);
grid minor;
hold on
axis([1e3 1e9 -180 -20]);
set(gcf,'unit','normalized','position',[0.3,0.3,0.5,0.5]);
set(gca,'YColor','k');
set(gca,'FontWeight','Bold');
set(gca,'FontSize',12);
xlabel('f_{offset}(Hz)');
ylabel('Phase Noise(dBc/Hz)');

% calc jitter
index_range = find((fxx>1e3)&(fxx<100e6));
frange = fxx(index_range);
PNrange = PNxx(index_range);
PNint = 2*trapz(frange,10.^(PNrange/10));
jitter = (1/(2*pi/freq_osc))*sqrt(PNint);

% figure: jitter vs BW
fprintf('PNint = %f dBc \n', 10*log10(PNint)); % SSB
fprintf('jitter = %f fs \n', jitter*1e15);

%% calc use pwelch
[freq,psd,freq_osc,jitter_total,jitter_random,spur_idx,spur_level] = pnpwelch(tcross, 'seg_num', 8, 'spur_fmax', 100e3, 'jitter_fmin', 1e3, 'jitter_fmax', 100e6, 'spur_sens', 6);

figure(2);
semilogx(freq, pow2db(psd), 'r','LineWidth',2);
grid minor;
hold on
axis([1e3 1e9 -180 -20]);
set(gcf,'unit','normalized','position',[0.3,0.3,0.5,0.5]);
set(gca,'YColor','k');
set(gca,'FontWeight','Bold');
set(gca,'FontSize',12);
xlabel('f_{offset}(Hz)');
ylabel('Phase Noise(dBc/Hz)');