% MLLsimpleDetune.m
% simple model of mode locked laser with detuning

clear all
close all

global Ts;              % sampling period
global Fcar;            % carrier frequency (optical frequency)
c_const = 3e8;          % speed of light

lamda = 1550e-9;        % m
Fcar = c_const/lamda;

Ts = 1e-13;
N = 1024;               % number of samples in a block: Tblk = N * Ts = 102.4 ps

% Amplifier parameters:
    GssdB = 20;         % (dB)
    PoutsatdB = 10;     % (dBm)
    NF = 8;             % (dB)
    
% Filter bandwidth:
    lamda3dB = 1.2e-9;  % m
    f3dB = lamda3dB*c_const/lamda^2;
    
% Modulator parameter:
    alpha = -0.007;
    epsilon = 40;       % (dB) extinction ratio
    
% Modulation parameters:
    m = 0.5;            % modulation index
    fm = 1e10;          % modulation frequency
    NHar = 1000;        % harmonic order
    Ts = 1/fm/N;        % recalculate Ts so that Tm = N*Ts
    f_detune = 0.5e5;   % detuned frequency
    delay_per_pass = round(NHar*f_detune/fm^2/Ts);
    
% Loss:
    loss = 10;          % dB
    atten = 1/10^(loss/20);
    % calculate the number of samples at the ends of the window must be
    % attenuatted
    p = floor((N*Ts*fm - 1)/2);
    if (p<0)
        Ntrunc = 1;
    else
        Ntrunc = round((N - (2*p+1)/fm/Ts)/2);
    end

% generation an initial block of signal Ein
% Ein = 1e-3*gausswin(N,2);
Ein = wgn(N,1,-40,'complex');

Eout = Ein;
Eo = Ein;
N_pass = 500;
for ii = 1:N_pass
  % [Eo,G] = amp_simp(Eo,GssdB,PoutsatdB,NF);
    [Eo,G] = AmpSimpNonoise(Eo,GssdB,PoutsatdB);    % no noise
    Eo = fft(Eo);
    Eo = filter_gaus(Eo,f3dB,1);
    Eo = ifft(Eo);
    
    if (delay_per_pass>0)
        Et = Eo(N - delay_per_pass + 1:N);
        Eo(delay_per_pass + 1:N) = Eo(1:N-delay_per_pass);
        Eo (1:delay_per_pass) = Et;
    else
        Et = Eo(1:-delay_per_pass);
        Eo(1:N + delay_per_pass) = Eo(1 - delay_per_pass:N);
        Eo(N+delay_per_pass+1:N) = Et;
    end
    
  % Eo = modInt(Eo(1:N),alpha,epsilon,m,fm,0.5);        % MZI modulator
    Eo = modInt_theory(Eo(1:N),m,fm);                   % Gaussian approximation modulator
    Eo = Eo*atten;
    
    
    
    [PeakVal MaxInd] = max(Eo.*conj(Eo));
    PulsePeak(ii) = PeakVal;
    PulsePeakPos(ii) = MaxInd;
    GOut(ii) = G;
    
    % ------ End ------ output of the peak power    
    if mod(ii,N_pass/50)==0
        Eout = [Eout , Eo];
    end
    if (ii>N_pass - 10)
        Eout = [Eout , Eo];
    end
end

% Eout = Eout/atten;
% ---------------- Display the results ------------------------------------
% mesh (abs(Eout'),'edgecolor','black','meshstyle','row','facecolor','none');
Iout = Eout.*conj(Eout);
mesh (Iout','meshstyle','row','facecolor','none');
axis tight;
% set(gca,'XTick',tt_mark);
% set(gca,'XTickLable',tt_tick);
% set(gca,'XDir','reverse');
xlabel('T (0.1ps)')
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
% chirp = pulse_beta/pulse_alpha;
% ------- plot spectrum of the lastest pulse ------------------------------
% Kmag = 8;
% Nplot = 100;
% Eoutfreq = fft(Eout(:,N1),N*Kmag);
% Ioutfreq = Eoutfreq.*conj(Eoutfreq)/(N*Kmag)^2;
% figure(2);
% ind = (- Nplot/2 : Nplot/2);
% freq = ind/Ts/N/Kmag;
% ind = mod((ind + N*Kmag),N*Kmag) + 1;
% plot(freq,Ioutfreq(ind));

% n = 1;
% n = 2*n;
% Tfil = exp(-log(2)*(2/f3dB*freq).^n);           % n order gaussian
% filter VPI
% Tfil = Tfil*max(Ioutfreq(ind));
% hold on
% plot(freq,Tfil,'r');
% ----------- End ------------ plot spectrum
% plot the gaussion fit curve
% gaussionFit(Iout(:,N1));
pulse_position = (find(Iout(:,N1)==max(Iout(:,N1))) - N/2)*Ts;
% calculation using the theoretical equations
wd = (2*pi*f3dB);
be = 2*log(2)/wd^2;
wm = 2*pi*fm;
al = 1/8*wm*(m*wm + sqrt(m^2*wm^2 + 2*m*wd^2/log(2)));
TpCal = sqrt(2*log(2)/al);
sprintf('Calculation using Li 2000 with correacted m: Tp = %e',TpCal);
td = NHar*f_detune/fm^2;
ts = td/4/al/be;
x = fm/f3dB;
delta_max_Li = sqrt(m*2*log(2))*x/4/NHar/(1+x*sqrt(m*2*log(2)))*fm;
Gss = 10^(GssdB/10);
% delta_max_Go = sqrt(2*log(2)*log(sqrt(Gss)*atten/sqrt(1+4*al*be)))*x/pi/NHar*fm;
% delta_max_2Tp = (0.589*sqrt(m)*x/NHar - 0.575*sqrt(sqrt(m))*x^1.5/NHar)*fm;
delta_max_Go = sqrt(4*be*log(sqrt(Gss)*atten/sqrt(1+4*al*be)))/NHar*fm;
delta_max_2Tp = 2*al*be/NHar/(1+4*al*be);

% plot the peak
figure (2)
subplot(2,1,1);
plot(PulsePeak);
subplot(2,1,2);
plot(PulsePeakPos - N/2);


