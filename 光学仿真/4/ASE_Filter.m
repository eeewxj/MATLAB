clear all
format long e
global Ts;
% Constant
c_const = 3e8;              % speed of light

% Time parameters:
n = 1024;                   % number of samples in a time period
num = 127;                  % observed time period number
N = n*num;                  % Total sampling points
fm = 1e10;                  % repetition rate (Hz)
Ts = 1/fm/n;                % sampling distance (s)
t = ((1:N)'-(N+1)/2)*Ts;    % vector of t values (s)

% Frequency parameters:
dv = 1/Ts/N;                                % frequency resolution (Hz)
v = [(0:N/2-1),(-N/2:-1)]'/(Ts*N);          % vector of v values (Hz)
vs = fftshift(v);                           % swap halves for plotting 

% Wavelength parameters:
lamda0 = 1554.8e-9;         % reference wavelength (m)
lamdaF = 1543.3e-9;         % carrier wavelength (m)
vw = fliplr(fftshift(c_const./(v + c_const/lamda0))'*1e9);  % wavelengths vector (nm)

% Filter bandwidth:
lamda3dB = 10e-9;          % m
f3dB = lamda3dB*c_const/lamdaF^2;

% Filter centra carrier
deltavF = c_const/lamdaF-c_const/lamda0;                    % Hz
[~,aux1b] = min(abs(v - deltavF));
CarrierF = exp(-1j*2*pi*v(aux1b)*t);

% White band noise
Es = wgn(N,1,-40,'complex');

% Filter @ lamdaF
Esf = conj(CarrierF).*Es;
Ssf = ifft(Esf,N)*N;
Ssf = filter_gaus(Ssf,f3dB,4);
%Ssf = filter_rect(Ssf,f3dB);
Eout = fft(Ssf)/N.*CarrierF;
Spectrum = fliplr(fftshift((abs(ifft(Eout)*N).^2)/N^2)');
S_dBm = 10*log10(Spectrum);

% Drawing
figure(1)
plot(vw,S_dBm);
axis([1520 1590 -inf 0]);
xlabel('Wavelength (nm)');
ylabel('Intensity (dBm)');
