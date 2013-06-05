% Pulse pump FOPA
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
    Ip = 18;                                   % pump power (W)
    Is = 0.01;                                  % signal power (W)
    
% Fiber parameters:
    dzN = 3;           % distance stepsize (m)
    nzN = 5;           % number of z-steps
    loss = 10;          % dB
    atten = 1/10^(loss/20);
    
% Detune parameters:
    fm = 1e10;          % modulation frequency (Hz)
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
    
for kk = 1:25
	delay_per_pass = kk*10-130;
    Eo = sqrt(Is)*Shape;
    Eo = circshift(Eo,delay_per_pass);
    Eo = CarrierS.*Eo;
    Eo = SSFM_FOPA(Eo,Ep,dzN,nzN);
    Eo = conj(CarrierS).*Eo;
    Eo = fft(Eo,N);
    Eo = ifft(Eo);
    Eo = Eo*atten;
    Ealp = Eo.*CarrierS./CarrierP;
    Eall = fft(Ealp,N);
    Iall(:,kk) = (10^kk)*Eall.*conj(Eall)/N^2;
end
figure(4);
plot(vwp,fftshift(10*log(Iall)));
axis([1517 1598 -280 530]);
xlabel('Wavelength (nm)');
ylabel('Intensity (dBm)');

