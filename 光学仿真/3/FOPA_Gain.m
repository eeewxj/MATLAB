clear;clc;
format long e
% Constant
c_const = 299792458e-12;              % speed of light
tol = 1e-5;
% Fiber Length
dz = 1;                    % distance stepsize (m)
nz = 20;                    % number of z-steps

% Time parameters:
n = 2^8;                   % number of samples in a time period (/ps)
num = 2^8;                  % observed time period number (ps)
N = n*num;                  % Total sampling points
Ts = 1/n;                   % sampling distance (ps)
t = ((1:N)'-(N+1)/2)*Ts;    % vector of t values (ps)

% Frequency parameters:
dv = 1/Ts/N;                                % frequency resolution (THz)
v = [(0:N/2-1),(-N/2:-1)]'/(Ts*N);          % vector of v values (THz)
vs = fftshift(v);                           % swap halves for plotting 

% Wavelength parameters:
lamda0 = 1064.0e-9;         % zero-dispersion wavelength (m)
lamdaP = 1064.3e-9;         % pump wavelength (m)
vw = flipud(fftshift(c_const./(v + c_const/lamda0))*1e9);  % wavelengths vector (nm)

% singal and pump power:
Is = 1e-9;                	% small signal power (W)
Ip = 20;                     % pump power (W)

% generation of comb like signal:
gainpoints = 200;           % number of gain points to be plotted
stepfreq = floor(N/(gainpoints + 1));       % determines the indices of frequencies used for plot
Es = zeros(N,1);

for ii = 1:gainpoints;
    Es = Es + sqrt(Is).*exp(-1i*2*pi*vs(N/2 + stepfreq*(ii-gainpoints/2) - floor(stepfreq/2))*t);
end

% Optical carrier
deltavP = c_const/lamdaP-c_const/lamda0;	% Hz
[~,aux2b] = min(abs(v - deltavP));
CarrierP = exp(-1j*2*pi*v(aux2b)*t);
% Coupled together
u0 = sechpulse(t,0,15,1,8.96);
%Ep = sqrt(Ip).*u0.*CarrierP;
Ep = sqrt(Ip).*CarrierP;
Ein = Es + Ep;
spectralpower0x = fftshift((abs(ifft(Ein)*N).^2));         % W

% Transmission in HNL-DSF
alpha = 1e-6;
beta = [0 0 8.19e-6 6.17e-5 -7.41e-8 -1.54e-11];
gamma = 0.018;
Eout = ssrklem(Ein,Ts,dz,dz*nz,alpha,beta,gamma,lamda0/c_const,tol);
%Eout = SSFM_Fiber(Ein,dz,nz);
spectralpowerx = fftshift((abs(ifft(Eout)*N).^2));         % W
gains = spectralpowerx./spectralpower0x;

% Select those discrete gain points
selectedgains = zeros(gainpoints,1);
selectedfreq = zeros(gainpoints,1);
ii = 1:gainpoints;
selectedfreq(ii) = vs(N/2 + stepfreq*(ii-gainpoints/2) - floor(stepfreq/2));
selectedgains(ii) = gains(N/2 + stepfreq*(ii-gainpoints/2) - floor(stepfreq/2));


% Visualization
figure(1);
vp = flipud((c_const./(selectedfreq + c_const/lamda0)).*1e9);     % wavelengths for plot (nm)
Gains_in_dB = flipud(10*log10(selectedgains));

plot(vp,Gains_in_dB,'r');
axis([vp(1) vp(gainpoints) -30 100]);
xlabel('Wavelength (nm)');
ylabel('Gain (dB)');
title_string = ['\lambda_p=',num2str(lamdaP*1e9,'%5.1fnm'),...
    ',L=',num2str(dz*nz,'%3.0fm'),...
    ',P_p_e_a_k=',num2str(Ip,'%2.0fW'),...
    ',tol=',num2str(tol,'%1.0e')];
title (title_string);
figure(2);
plot(t,abs(Eout).^2)
title (title_string);


