clear all
format long e
global Ts;              % s
% CONSTANTS
c_const = 3e8;          % speed of light

% NUMERICAL PARAMETERS
n = 1024;               % number of samples in a block
num = 41;               % observed time period number
N = 1024*num;
fm = 1e10;              % modulation frequency (Hz)
Ts = 1/fm/n;            % recalculate Ts so that Tm = N*Ts (s)
dzN = 15;               % distance stepsize (m)
nzN = 10;               % number of z-steps

% Wavelength parameters:
lamda0 = 1554e-9;           % m
lamdaS = 1540e-9;           % m
lamdaP = 1556e-9;           % m
deltavP = c_const/lamda0-c_const/lamdaP;	% Hz
Ip = 0.60;                                  % pump power (W)

gainpoints = 200;                           % number of gain points to be plotted
input_signal_power = 1e-9;                  % W

% CALCULATED QUANTITIES 
t = ((1:N)'-(N+1)/2)*Ts;                    % vector of t values (s)
v = [(0:N/2-1),(-N/2:-1)]'/(Ts*N);          % vector of v values (Hz)
vw = fftshift(c_const./(v + c_const/lamdaP))'*1e9;      % wavelengths (nm)
vs = fftshift(v);                           % swap halves for plotting 

% Frequency parameters:
[aux2a,aux2b] = min(abs(v - deltavP));
CarrierP = exp(-1j*2*pi*v(aux2b)*t);
Ep = sqrt(Ip)*CarrierP;
    
%STARTING FIELD
%------adding input pumps in time domain-----------------------------------------------
stepfreq = floor(N/(gainpoints + 1));       % determines the indices of the frequencies used for plot
Eo = zeros(N,1);
for ii = 1:gainpoints/2;
    Eo = Eo + sqrt(input_signal_power).*exp(j*2*pi*vs(N/2 + stepfreq*ii)*t);
    Eo = Eo + sqrt(input_signal_power).*exp(j*2*pi*vs(N/2 - stepfreq*ii + floor(stepfreq/2))*t);
end
%PROPAGATE
tx=cputime;
Eout = SSFM_FOPA(Eo,Ep,dzN,nzN);
fprintf(1,'\nCompleted in (%.2f s)\n',cputime-tx);

%VISUALIZATION
%----Plot Gain Using SSFM result---------------------------------------------
spectralpower0x = fftshift((abs(fft(Eo)).^2)/N^2);          % W
spectralpowerx = fftshift((abs(fft(Eout)).^2)/N^2);         % W

gains = spectralpowerx./spectralpower0x;
selectedgains = [];
selectedfreq = [];

for ii = 1:gainpoints/2;
    selectedfreq = [selectedfreq, vs(N/2 - stepfreq*ii + floor(stepfreq/2))];
    selectedgains = [selectedgains, gains(N/2 - stepfreq*ii + floor(stepfreq/2))];
end
selectedfreq = sort(selectedfreq);
selectedgains = selectedgains(length(selectedgains):-1:1);

for ii = 1:gainpoints/2;
    selectedfreq = [selectedfreq, vs(N/2 + stepfreq*ii)];
    selectedgains = [selectedgains, gains(N/2 + stepfreq*ii)];
end

vp = [c_const./ (selectedfreq + c_const/lamda0)]'.*1e9;     % wavelengths for plot (nm)
Gains_in_dB = 10*log10(selectedgains);

plot(vp,Gains_in_dB,'r');
axis([1523 1590 0 20]);
xlabel('Wavelength (nm)');
ylabel('Gain (dB)');