% MLLsimple.m
% simplemode of mode locked laser

clear all
close all

global Ts;              % sampling period
global Fcar;            % carrier frequency (optical frequency)
c_const = 3e8;          % speed of light

lamda = 1550e-9;        % m
Fcar = c_const/lamda;   % Hz

n = 1024;               % number of samples in a block
num = 11;               % observed time period number
N = 1024*num;

% Amplifier parameters:
    GssdB = 20;         % (dB)
    PoutsatdB = 10;     % (dBm)
    NF = 8;             % (dB)
    
% Filter bandwidth:
    lamda3dB = 0.8e-9;    % m
    f3dB = lamda3dB*c_const/lamda^2;
    
% Modulator parameter:
    alpha = -0.07;
    epsilon = 40;       % (dB) extinction ratio
    
% Modulation parameters:
    m = 0.5;            % modulation index
    fm = 1e10;          % modulation frequency
    NHar = 100;         % harmonic order
    Ts = 1/fm/n;        % recalculate Ts so that Tm = N*Ts
    f_detune = 1e6;     % detuned frequency
    delay_per_pass = round(NHar*f_detune/fm^2/Ts);
    
% Loss:
    loss = 10;          % dB
    atten = 1/10^(loss/20);

% generation an initial block of signal Ein
Ein = wgn(N,1,-40,'complex');
Eout = Ein;
Eo = Ein;

N_pass = 500;
for ii = 1:N_pass
    [Eo,G] = amp_simp1(Eo,GssdB,PoutsatdB,NF,num);
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

%figure(1);
%mesh (Iout','meshstyle','row','facecolor','none');
%axis tight;
%xlabel('T (0.1ps)');
%ylabel('Pass number');
%zlabel('intensity (W)');

N1 = size(Eout,2);
Kmag = 1;
Nplot = 100;
Eoutfreq = fft(Eout(:,N1),N*Kmag);
Ioutfreq = Eoutfreq.*conj(Eoutfreq)/(N*Kmag)^2;
figure(2);
plot(fftshift(10*log(Ioutfreq)));
%ind = (- Nplot/2 : Nplot/2);
%freq = ind/Ts/N/Kmag;
%ind = mod((ind + N*Kmag),N*Kmag) + 1;
%plot(freq,Ioutfreq(ind));

figure(3);
Trace = Iout(:,N1);
for jj = 1:num
    hold on;
    plot(Trace(n*(jj-1)+1:n*jj));
end
axis tight;
xlabel('T (0.1ps)');
ylabel('intensity (W)');
