% MLLsimple.m
% simplemode of mode locked laser

clear all
close all

global Ts;              % sampling period
global Fcar;            % carrier frequency (optical frequency)
c_const = 3e8;          % speed of light

n = 1024;               % number of samples in a block
num = 64;               % observed time period number
N = 1024*num;

% Wavelength parameters:
    lamda0 = 1550e-9;           % m
    lamdaC = 1550e-9;           % m
    Fcar = c_const/lamdaC;      % Hz
    deltav = c_const/lamda0-c_const/lamdaC;     % Hz
    
% Amplifier parameters:
    GssdB = 20;         % (dB)
    PoutsatdB = 10;     % (dBm)
    NF = 8;             % (dB)
    
% Fiber parameters:
    dz = 10;             % distance stepsize (m)
    nz = 3;            % number of z-steps
    
% Filter bandwidth:
    lamda3dB = 1.2e-9;  % m
    f3dB = lamda3dB*c_const/lamdaC^2;
    
% Modulator parameter:
    alpha = -0.07;
    epsilon = 10;       % (dB) extinction ratio
    
% Modulation parameters:
    m = 0.5;            % modulation index
    fm = 1e10;          % modulation frequency (Hz)
    NHar = 1000;        % harmonic order
    Ts = 1/fm/n;        % recalculate Ts so that Tm = N*Ts (s)
    f_detune = 0;   % detuned frequency (Hz)
    l_detune = c_const*NHar*f_detune/fm^2;      % detuned cavity length (m)
    delay_per_pass = round(l_detune/c_const/Ts);
    
% Frequency parameters:
    t = ((1:N)'-(N+1)/2)*Ts;                    % vector of t values (s)
    v = [(0:N/2-1),(-N/2:-1)]'/(Ts*N);          % vector of v values (Hz)
    [aux1a,aux1b] = min(abs(v - deltav));
    Carrier = exp(-1j*2*pi*v(aux1b)*t);
    
% Loss:
    loss = 10;          % dB
    atten = 1/10^(loss/20);

% generation an initial block of signal Ein
Ein = wgn(N,1,-40,'complex');
Eout = Ein;
Eo = Ein;

tx = cputime;
N_pass = 100;
for ii = 1:N_pass
    [Eo,G] = amp_simp1(Eo,GssdB,PoutsatdB,NF,num);
    Eo = Carrier.*Eo;
    Eo = SSFM_SMF(Eo,dz,nz);
    Eo = conj(Carrier).*Eo;
    Eo = fft(Eo,N);
    Eo = filter_gaus(Eo,f3dB,1);
    Eo = ifft(Eo);
    Eo = circshift(Eo,delay_per_pass);
    Eo = modInt(Eo(1:N),alpha,epsilon,m,fm,0.5);
    Eo = Eo*atten;
    % ------ End ------ output of the peak power
    if mod(ii,N_pass/50)==0
        Eout = [Eout , Eo];
    end
    Eo = Eo + wgn(N,1,-40,'complex');
end

Eout = Eout/atten;
close all
Iout = Eout.*conj(Eout);
fprintf(1,'\nSimulation completed in (%.2f s)\n',cputime-tx);

figure(1);
mesh (Iout','meshstyle','row','facecolor','none');
axis tight;
xlabel('T (0.1ps)');
ylabel('Pass number');
zlabel('Intensity (W)');

N1 = size(Eout,2);
Kmag = 1;
Nplot = 100;
Eoutfreq = fft(Eout(:,N1),N*Kmag);
Ioutfreq = Eoutfreq.*conj(Eoutfreq)/(N*Kmag)^2;
vw = fftshift(c_const./(v + c_const/lamdaC))'*1e9;  % wavelengths (nm)
figure(2);
subplot(2,1,2);
plot(vw,fftshift(10*log(Ioutfreq)));
axis([1549 1551 -200 -40]);
xlabel('Wavelength (nm)');
ylabel('Intensity (dBm)');

%figure(3);
subplot(2,1,1);
Trace = Iout(:,N1);
Ptime = Ts*(1:n)*1e12;
for jj = 1:num
    hold on;
    plot(Ptime,Trace(n*(jj-1)+1:n*jj));
end
axis([0 100 0 0.3]);
xlabel('Time (ps)');
ylabel('Intensity (W)');
