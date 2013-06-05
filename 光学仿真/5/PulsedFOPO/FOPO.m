% Pulse pump FOPO
clear all
close all

global Ts;              % sampling period
global Fcar;            % carrier frequency (optical frequency)

c_const = 3e8;          % speed of light (m/s)

n = 1024;               % number of samples in a block
num = 11;               % observed time period number
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
    Ip = 1.2;                                   % pump power (W)
    
% Fiber parameters:
    dz = 2;             % distance stepsize (m)
    nz = 6;             % number of z-steps
    dzN = 15;           % distance stepsize (m)
    nzN = 10;           % number of z-steps
    dzD = 6;           % distance stepsize (m)
    nzD = 0;           % number of z-steps
    
% Filter bandwidth:
    lamda3dB = 1e-9;    % m
    f3dB = lamda3dB*c_const/lamdaS^2;
    WDMC = 15e-9;        % m
    Band = WDMC*c_const/lamdaS^2;
    lamda3dBout = 1.0e-9;  % m
    f3dBout = lamda3dBout*c_const/lamdaP^2;
    loss = 10;          % dB
    atten = 1/10^(loss/20);
    
% Detune parameters:
    fm = 1e10;          % modulation frequency (Hz)
    NHar = round((dz*nz+dzN*nzN+dzD*nzD)*fm/c_const);       % harmonic order
    Ts = 1/fm/n;        % recalculate Ts so that Tm = N*Ts (s)
    f_detune = 3.8e4;   % detuned frequency (Hz)
    l_detune = c_const*NHar*f_detune/fm^2;                  % detuned cavity length (m)
    %l_detune = -0.055;
    delay_per_pass = round(l_detune/c_const/Ts);
    
% Frequency parameters:
    t = ((1:N)'-(N+1)/2)*Ts;                                % vector of t values (s)
    v = [(0:N/2-1),(-N/2:-1)]'/(Ts*N);                      % vector of v values (Hz)
    vws = fftshift(c_const./(v + c_const/lamdaS))'*1e9;     % wavelength (nm)
    vwp = fftshift(c_const./(v + c_const/lamdaP))'*1e9;     % wavelength (nm)
    vwi = fftshift(c_const./(v + c_const/lamdaI))'*1e9;     % wavelength (nm)
    [aux1a,aux1b] = min(abs(v - deltavS));
    CarrierS = exp(-1j*2*pi*v(aux1b)*t);
    [aux2a,aux2b] = min(abs(v - deltavP));
    CarrierP = exp(-1j*2*pi*v(aux2b)*t);
    [aux3a,aux3b] = min(abs(v - deltavI));
    CarrierI = exp(-1j*2*pi*v(aux3b)*t);

% Generation of initial pump Ep
    Pulsewidth = 20e-12;                                    % OPA pump pulsewidth (s)
    Base = exp(-2*log(2)*(t./Pulsewidth).^2);
    Shape = Base;
    for pp = 1:num-1
        Base = circshift(Base,n);
        Shape = Shape + Base;
    end
    Ep = sqrt(Ip)*Shape.*CarrierP;
    
% Generation of initial signal Ein
    Ein = wgn(N,1,-40,'complex');
    Eout = Ein;
    Eo = Ein;

    tx = cputime;
    N_pass = 100;
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
    % ------ End ------ output of the peak power
    if mod(ii,N_pass/50)==0
        Eout = [Eout , Eo];
    end
    Eo = Eo + wgn(N,1,-40,'complex');
end

Eout = Eout/atten;
Iout = Eout.*conj(Eout);
fprintf(1,'\nSimulation completed in (%.2f s)\n',cputime-tx);

%figure(1);
%mesh (Iout','meshstyle','row','facecolor','none');
%axis tight;
%xlabel('T (0.1ps)');
%ylabel('Pass number');
%zlabel('Intensity (W)');

figure(3);
% Signal
Eals = circshift(Eoo,0);
Esps = fft(Eals,N);
Esps = filter_gaus(Esps,f3dBout,1);
Eals = ifft(Esps);
Isps = Esps.*conj(Esps)/N^2;
subplot(2,3,4);
plot(vws,fftshift(10*log(Isps)));
axis([1537 1543 -250 -40]);
xlabel('Wavelength (nm)');
ylabel('Intensity (dBm)');
subplot(2,3,1);
Ptime = Ts*(1:n)*1e12;
Tsignal = Eals.*conj(Eals);
for jj = 1:num
    hold on;
    plot(Ptime,Tsignal(n*(jj-1)+1:n*jj));
end
axis([0 100 0 0.2]);
title('Signal');
xlabel('Time (ps)');
ylabel('Intensity (W)');

% Pump
Ealp = circshift(Eoo,0).*CarrierS./CarrierP;
Espp = fft(Ealp,N);
Espp = filter_gaus(Espp,f3dBout,1);
Ealp = ifft(Espp);
Ispp = Espp.*conj(Espp)/N^2;
subplot(2,3,5);
plot(vwp,fftshift(10*log(Ispp)));
axis([1553 1559 -250 -40]);
xlabel('Wavelength (nm)');
ylabel('Intensity (dBm)');
subplot(2,3,2);
Tpump = Ealp.*conj(Ealp);
for jj = 1:num
    hold on;
    plot(Ptime,Tpump(n*(jj-1)+1:n*jj));
end
axis([0 100 0 0.8]);
title('Pump');
xlabel('Time (ps)');
ylabel('Intensity (W)');

% Idler
Eali = circshift(Eoo,0).*CarrierS./CarrierI;
Espi = fft(Eali,N);
Espi = filter_gaus(Espi,f3dBout,1);
Eali = ifft(Espi);
Ispi = Espi.*conj(Espi)/N^2;
subplot(2,3,6);
plot(vwi,fftshift(10*log(Ispi)));
axis([1569 1575 -250 -40]);
xlabel('Wavelength (nm)');
ylabel('Intensity (dBm)');
subplot(2,3,3);
Tidle = Eali.*conj(Eali);
for jj = 1:num
    hold on;
    plot(Ptime,Tidle(n*(jj-1)+1:n*jj));
end
axis([0 100 0 0.2]);
title('Idler');
xlabel('Time (ps)');
ylabel('Intensity (W)');

figure(4);
Ealp = Eoo.*CarrierS./CarrierP;
Eall = fft(Ealp,N);
Iall = Eall.*conj(Eall)/N^2;
plot(vwp,fftshift(10*log(Iall)),'b');
axis([1523 1590 -170 -20]);
xlabel('Wavelength (nm)');
ylabel('Intensity (dBm)');

