
clear all
close all
global Ts;

% CONSTANTS
    c = 299792458;                          % speed of light (m/s)
    lamda_IN = 1550e-9;                     % m
    lamda0 = 1550e-9;                       % m

% Fiber parameters:
    dz = 50;                               % distance stepsize (m)
    nz = 20;                                % number of z-steps
    Pulsewidth = 6e-12;                     % s
    diffe = 20e-12;                         % s

% CALCULATED QUANTITIES 
    N = 2^17;                              % number of points in FFT
    Ts = 5e-14;                             % time step (s)
    T = N*Ts;                              % FFT window size (s)
    t = ((1:N)'-(N+1)/2)*Ts;              % vector of t values (s)
    v = [(0:N/2-1),(-N/2:-1)]'/T;         % vector of v values (Hz)
    
deltav = c/lamda0-c/lamda_IN;
[aux1a,aux1b] = min(abs(v - deltav));
Ein = exp(-1j*2*pi*v(aux1b)*t).*(exp(-2*log(2)*(t./Pulsewidth).^2) + 0.6.*exp(-2*log(2)*((t-diffe)./Pulsewidth).^2));
Eout = SSFM_SMF(Ein,dz,nz);

figure(1)
subplot(2,1,1)
    plot(t*1e12,abs(Ein));
    axis([-50 50 -inf inf]);
subplot(2,1,2)
    plot(t*1e12,abs(Eout));
    axis([-50 50 -inf inf]);



