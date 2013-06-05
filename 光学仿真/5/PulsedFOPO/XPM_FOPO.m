% Pulse pump FOPO with XPM effect
clear all
close all

global Ts;              % sampling period
global Fcar;            % carrier frequency (optical frequency)

c_const = 3e8;          % speed of light (m/s)

n = 1024;               % number of samples in a block
num = 1;               % observed time period number
N = 1024*num;

% Wavelength parameters:
    lamda0 = 1554e-9;           % m
    lamdaS = 1540e-9;           % m
    lamdaP = 1556e-9;           % m
    lamdaI = lamdaS*lamdaP/(2*lamdaS-lamdaP);   % m
    Fcar = c_const/lamdaS;                      % Hz
    deltavS = c_const/lamda0-c_const/lamdaS;	% Hz
    deltavP = c_const/lamda0-c_const/lamdaP;	% Hz
    deltavI = c_const/lamda0-c_const/lamdaI;	% Hz
    Ip = 1.5;                                   % pump power (W)
    
% Fiber parameters:
    dz = 2;             % distance stepsize (m)
    nz = 6;             % number of z-steps
    dzN = 15;           % distance stepsize (m)
    nzN = 10;           % number of z-steps
    dzD = 3;           % distance stepsize (m)
    nzD = 0;           % number of z-steps
    
% Filter bandwidth:
    lamda3dB = 1e-9;    % m
    f3dB = lamda3dB*c_const/lamdaS^2;
    WDMC = 15e-9;        % m
    Band = WDMC*c_const/lamdaS^2;
    lamda3dBout = 2.0e-9;  % m
    f3dBout = lamda3dBout*c_const/lamdaP^2;
    loss = 18;          % dB
    atten = 1/10^(loss/20);
    
% Detune parameters:
    fm = 1e10;          % modulation frequency (Hz)
    NHar = round((dz*nz+dzN*nzN+dzD*nzD)*fm/c_const);       % harmonic order
    Ts = 1/fm/n;        % recalculate Ts so that Tm = N*Ts (s)
    
% Frequency parameters:
    t = ((1:N)'-(N+1)/2)*Ts;                                % vector of t values (s)
    v = [(0:N/2-1),(-N/2:-1)]'/(Ts*N);                      % vector of v values (Hz)
    vws = fftshift(c_const./(v + c_const/lamdaS))'*1e9;     % wavelength (nm)
    vwp = fftshift(c_const./(v + c_const/lamdaP))'*1e9;     % wavelength (nm)
    [aux1a,aux1b] = min(abs(v - deltavS));
    CarrierS = exp(-1j*2*pi*v(aux1b)*t);
    [aux2a,aux2b] = min(abs(v - deltavP));
    CarrierP = exp(-1j*2*pi*v(aux2b)*t);

% Generation of initial pump Ep
    Pulsewidth = 5e-12;                                    % OPA pump pulsewidth (s)
    Base = exp(-2*log(2)*(t./Pulsewidth).^2);
    Shape = Base;
    for pp = 1:num-1
        Base = circshift(Base,n);
        Shape = Shape + Base;
    end
    Ep = sqrt(Ip)*Shape.*CarrierP;
    N_pass = 100;

for kk = 1:25
delay_per_pass = kk*5-60+21;
% Generation of initial signal Ein
Eo = wgn(N,1,-80,'complex');
for ii = 1:N_pass
    Eo = CarrierS.*Eo;
    Eo = SSFM_SMF(Eo,dz,nz);
    Eo = SSFM_FOPA(Eo,Ep,dzN,nzN);
    Eo = SSFM_DCF(Eo,dzD,nzD);
    Eoo = conj(CarrierS).*Eo;
    Eo = fft(Eoo,N);
    Eo = filter_gaus(Eo,f3dB,1);            % filter out signal noise
    Eo = filter_rect(Eo,Band);              % filter out the pump
    Eo = ifft(Eo);
    Eo = circshift(Eo,delay_per_pass);
    Eo = Eo*atten;
    Eo = Eo + wgn(N,1,-80,'complex');
end
Ealp = Eoo.*CarrierS./CarrierP;
Eall = fft(Ealp,N);
Iall(:,kk) = (10^kk)*Eall.*conj(Eall)/N^2;
end


figure(4);
plot(vwp,fftshift(10*log(abs(Iall))));
axis([1518 1595 -250 530]);
xlabel('Wavelength (nm)');
ylabel('Intensity (dBm)');

