clear all
format long e
global v c_const Ts;
% Constant
c_const = 3e8;              % speed of light

% Time parameters:
n = 2048;                   % number of samples in a time period
num = 64;                   % observed time period number
N = n*num;                  % Total sampling points
fm = 1e10;                  % repetition rate (Hz)
Ts = 1/fm/n;                % sampling distance (s)
t = ((1:N)'-(N+1)/2)*Ts;    % vector of t values (s)

% Frequency parameters:
dv = 1/Ts/N;                                % frequency resolution (Hz)
v = [(0:N/2-1),(-N/2:-1)]'/(Ts*N);          % vector of v values (Hz)
vs = fftshift(v);                           % swap halves for plotting 

% Wavelength parameters:
lamda0 = 1554.8e-9;     	% m
lamdaS = 1543.3e-9;         % m
vw = fliplr(fftshift(c_const./(v + c_const/lamda0))'*1e9);	% wavelengths (nm)

% Generation of initial pulsed signal Es:
Is = 1;                     % W
Pulsewidth = 1e-12;                                         % OPA pump pulsewidth (s)
Base = exp(-2*log(2)*(t./Pulsewidth).^2);
Shape = Base;
for pp = 1:num-1
    Base = circshift(Base,n);
    Shape = Shape + Base;
end

% Optical carrier
deltavS = c_const/lamda0-c_const/lamdaS;                    % Hz
[~,aux1b] = min(abs(v - deltavS));
CarrierS = exp(-1j*2*pi*v(aux1b)*t);

Es = sqrt(Is)*Shape.*CarrierS;
spectralpower = fliplr(fftshift((abs(fft(Es)).^2)/N^2)');	% W
S_dBm = 10*log10(spectralpower);
[Wavel1,Spect1] = OSA_Resol(1535,1555,spectralpower,0.02,lamda0);
[Wavel2,Spect2] = OSA_Resol(1535,1555,spectralpower,0.04,lamda0);
[Wavel3,Spect3] = OSA_Resol(1535,1555,spectralpower,0.08,lamda0);

% Figures
figure(1)
subplot(2,2,1);
plot(vw,S_dBm,'b');
axis([1535 1555 -100 20]);
xlabel('Wavelength (nm)');
ylabel('Intensity (dBm)');
title('Full resolution');
subplot(2,2,2);
plot(Wavel1,Spect1,'b');
axis([1535 1555 -100 20]);
xlabel('Wavelength (nm)');
ylabel('Intensity (dBm)');
title('0.02-nm resolution');
subplot(2,2,3);
plot(Wavel2,Spect2,'b');
axis([1535 1555 -100 20]);
xlabel('Wavelength (nm)');
ylabel('Intensity (dBm)');
title('0.04-nm resolution');
subplot(2,2,4);
plot(Wavel3,Spect3,'b');
axis([1535 1555 -100 20]);
xlabel('Wavelength (nm)');
ylabel('Intensity (dBm)');
title('0.08-nm resolution');

