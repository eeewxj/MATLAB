% MLLsimple.m
% simplemode of mode locked laser

clear all
close all

global Ts;              % sampling period
global Fcar;            % carrier frequency (optical frequency)
c_const = 3e8;          % speed of light

lamda = 1550e-9;        % m
Fcar = c_const/lamda;   % Hz

Ts = 1e-13;             % s
N = 1024;               % number of samples in a block: Tblk = N * Ts = 102.4 ps

% Amplifier parameters:
    GssdB = 20;         % (dB)
    PoutsatdB = 10;     % (dBm)
    NF = 8;             % (dB)
    
% Filter bandwidth:
    lamda3dB = 1e-9;    % m
    f3dB = lamda3dB*c_const/lamda^2;
    
% Modulator parameter:
    alpha = -0.07;
    epsilon = 40;       % (dB) extinction ratio
    
% Modulation parameters:
    m = 0.5;            % modulation index
    fm = 1e10;          % modulation frequency
    
% Loss:
    loss = 10;          % dB
    atten = 1/10^(loss/20);

% generation an initial block of signal Ein
% Ein = 1e-3*gausswin(N,2);
Ein = wgn(N,1,-40,'complex');

Eout = Ein;
Eo = Ein;
N_pass = 500;
for ii = 1:N_pass
    [Eo,G] = amp_simp(Eo,GssdB,PoutsatdB,NF);
  % [Eo,G] = AmpSimpNoise(Eo,GssdB,PoutsatdB);    % no noise
    Eo = fft(Eo,N*4);
  % Eo = filter_bessel(Eo,f3dB,G);
  % Eo = filter_gaus1(Eo,f3dB);
    Eo = filter_gaus(Eo,f3dB,1);
    Eo = ifft(Eo);
    Eo = modInt(Eo(1:N),alpha,epsilon,m,fm,0.5);
  % Eo = modPhase(Eo(1:N),m,fm);
    Eo = Eo*atten;
    if mod(ii,10)==0
        Eout = [Eout , Eo];
    end
end

Eout = Eout/atten;
close all
% mesh (abs(Eout'),'edgecolor','black','meshstyle','row','facecolor','none');
Iout = Eout.*conj(Eout);
mesh (Iout','meshstyle','row','facecolor','none');
axis tight;
% set(gca,'XTick',tt_mark);
% set(gca,'XTickLable',tt_tick);
% set(gca,'XDir','reverse');
xlabel('T (0.1ps)');
% set(gca,'YTick',yy_mark);
% set(gca,'YTickLable',yy_tick);
ylabel('Pass number');
zlabel('intensity (W)');

N1 = size(Eout,2);
dPhi = angle(Eout(2:N,N1)) - angle(Eout(1:N-1,N1));
% figure (2);
% plot(dPhi);

Tp = fwhm(Iout(:,N1))*Ts;
pulse_alpha = 2*log(2)/(Tp^2);
pulse_beta = (dPhi(N/2+100) - dPhi(N/2 - 100))/200/Ts/Ts;
chirp = pulse_beta/pulse_alpha;

Kmag = 8;
Nplot = 100;
Eoutfreq = fft(Eout(:,N1),N*Kmag);
Ioutfreq = Eoutfreq.*conj(Eoutfreq)/(N*Kmag)^2;
figure(2);
ind = (- Nplot/2 : Nplot/2);
freq = ind/Ts/N/Kmag;
ind = mod((ind + N*Kmag),N*Kmag) + 1;
plot(freq,Ioutfreq(ind));

n = 1;
n = 2*n;
Tfil = exp(-log(2)*(2/f3dB*freq).^n);           % n order gaussian
% filter VPI
Tfil = Tfil*max(Ioutfreq(ind));
hold on
plot(freq,Tfil,'r');

% plot the gaussion fit curve
% gaussionFit(Iout(:,N1));

pulseBW = fwhm(Ioutfreq(ind))/Ts/N/Kmag;
Tp = fwhm(Iout(:,N1))*Ts;
TBP = pulseBW*Tp;

