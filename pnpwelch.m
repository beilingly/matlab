function [freq,psd,freq_osc,jitter_total,jitter_random,spur_idx,spur_level]=pnpwelch(xtime,varargin)
%PNPWELCH Calculate approximated power density spectrum of phase noise.
%
%   See also PWELCH.
%
%   For more information, please contact <a href="matlab: 
%   web('mailto:yjy22@mails.tsinghua.edu.cn')">Yin Junyang</a>.

% check validity of cross-zero timestamp sequence
if(length(size(xtime))>2)
    error('The dimension of crossing-zero timestamp is not correct.');
end

if(size(xtime,2)~=1)
    xtime=xtime';
end

% read-in extra parameters
p=inputParser;
addRequired(p,'xtime');
addParameter(p,'trim_idx',1);
addParameter(p,'seg_num',1);
addParameter(p,'spur_fmax',1e6);
addParameter(p,'spur_sens',6);
addParameter(p,'jitter_fmin',10e3);
addParameter(p,'jitter_fmax',10e6);
addParameter(p,'freq_assigned',[]);
parse(p,'trim_idx',varargin{:});
parse(p,'seg_num',varargin{:});
parse(p,'spur_fmax',varargin{:});
parse(p,'spur_sens',varargin{:});
parse(p,'jitter_fmin',varargin{:});
parse(p,'jitter_fmax',varargin{:});
parse(p,'freq_assigned',varargin{:});

trim_idx=p.Results.trim_idx;
seg_num=p.Results.seg_num;
spur_fmax=p.Results.spur_fmax;
spur_sens=p.Results.spur_sens;
jitter_fmin=p.Results.jitter_fmin;
jitter_fmax=p.Results.jitter_fmax;
freq_assigned=p.Results.freq_assigned;

% trim timestamp sequence
xtime=xtime(trim_idx:end);
xtime=xtime-xtime(1);
smp_num=length(xtime);

% use 1st-order linear fitting on decentralized data sequence
[pcoeff,~,mu]=polyfit((0:smp_num-1)',xtime,1);
pcoeff(2)=pcoeff(2)-pcoeff(1)*mu(1)/mu(2);
pcoeff(1)=pcoeff(1)/mu(2);

% derive basic infomation
if(isempty(freq_assigned))
    period_osc=pcoeff(1);
else
    period_osc=1/freq_assigned;
end
time_offset=pcoeff(2);
freq_osc=1/period_osc;

% derive noise information
xtime_ideal=(0:smp_num-1)'.*period_osc+time_offset;
xtime_error=xtime-xtime_ideal;
jitter_total=rms(xtime_error);
phase_error=xtime_error./period_osc.*(2*pi);

% Caution: To ensure accurate spectral analysis, it is highly necessary
% to remove the mean of the time-domain data once again before spectral
% analysis, especially for oscillators with significant 1/f noise. Due to
% the finite simulation duration, the integral of flicker noise will not
% be zero, implying a slow drift in the oscillator's center frequency.
% Failing to remove the DC offset of phase_error may result in spectral
% analysis results resembling a frequency-modulated signal.
phase_error=phase_error-mean(phase_error);

% calculate phase noise spectrum
if (seg_num==1)
    window_len=smp_num;
else
    window_len=floor(smp_num/(seg_num-0.5));
end
[psd,freq]=pwelch(phase_error,hann(window_len),[],[],freq_osc,'psd');
dc_idx=(freq==0);
psd(~dc_idx)=0.5*psd(~dc_idx);%convert single-sided spectrum to two-sided

% derive integrated jitter
jitter_fmin_idx=find(freq>jitter_fmin,1);
jitter_fmax_idx=find(freq>jitter_fmax,1);
if(isempty(jitter_fmin_idx))
    jitter_fmin_idx=2;
end
if(isempty(jitter_fmax_idx))
    jitter_fmax_idx=length(freq);
end
jitter_random=sqrt(2*trapz(freq(jitter_fmin_idx:jitter_fmax_idx),psd(jitter_fmin_idx:jitter_fmax_idx)))/(2*pi)*period_osc;

% detect spurious tone
psd_log=10*log10(psd);
[spur_psd,spur_idx]=findpeaks(psd_log,'Threshold',spur_sens);
spur_psd=spur_psd((freq(spur_idx)<=spur_fmax));
spur_idx=spur_idx((freq(spur_idx)<=spur_fmax));
spur_level=zeros(1,length(spur_idx));
if(~isempty(spur_idx))
    for i=1:length(spur_idx)
%         spur_level(i)=spur_psd(i)+10*log10((freq(spur_idx(i)+1)-freq(spur_idx(i)-1))/2);
        spur_level(i)=spur_psd(i);
    end
end

end