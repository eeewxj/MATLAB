clear all
format long e
global Ts;
% Constant
c_const = 299792458;              % speed of light

% Fiber Length
dz = 15;                    % distance stepsize (m)
nz = 10;                    % number of z-steps

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
lamda0 = 1554.8e-9;         % zero-dispersion wavelength (m)
lamdaS = 1543.3e-9;         % signal wavelength (m)
lamdaP = 1556.5e-9;         % pump wavelength (m)
vw = fliplr(fftshift(c_const./(v + c_const/lamda0))'*1e9);  % wavelengths vector (nm)

% singal and pump power:
Is = 0.001;                	% signal power (W)
Ip = 1;                     % pump power (W)

% Optical carrier
deltavS = c_const/lamda0-c_const/lamdaS;                    % Hz
[~,aux1b] = min(abs(v - deltavS));
CarrierS = exp(-1j*2*pi*v(aux1b)*t);
deltavP = c_const/lamda0-c_const/lamdaP;                    % Hz
[~,aux2b] = min(abs(v - deltavP));
CarrierP = exp(-1j*2*pi*v(aux2b)*t);

% Coupled together
Es = sqrt(Is).*CarrierS;
Ep = sqrt(Ip).*CarrierP;
Ein = Es + Ep;

% Transmission in HNL-DSF
Eout = SSFM_Fiber(Ein,dz,nz);

spectraOUT = fliplr(fftshift((abs(fft(Eout)).^2)/N^2)');	% W
S_dBm = 10*log10(spectraOUT);

figure(1)
plot(vw,S_dBm,'b');
axis([1510 1600 -100 0]);
xlabel('Wavelength (nm)');
ylabel('Intensity (dBm)');
