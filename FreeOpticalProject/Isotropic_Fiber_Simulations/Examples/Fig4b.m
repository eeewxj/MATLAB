clear all
format long e
fprintf(1,'\n----------------------------------------------');
fprintf(1,'\nSimulating Figure 4b'); 

% CONSTANTS

c = 299792458;                          % speed of light (m/s)

% NUMERICAL PARAMETERS

nt = 2^14;                              % number of points in FFT
dt = 4e-2;                              % time step (ps)
dz = 0.4;                               % distance stepsize (m)
nz = 500;                               % number of z-steps

% OPTICAL PARAMETERS 

gamma = 1e-2;                           % W^-1 m^-1

beta1x = 0;                             % ps/m                                      
beta2x = 0;                             % ps^2/m
beta3x =  1e-4;                         % ps^3/m   
beta4x = 1e-7;                          % ps^4/m
alphax = 0;                             % 1/m

beta1y = 0;                             % ps/m   
beta2y = 0;                             % ps^2/m
beta3y =  1e-4;                         % ps^3/m
beta4y = 1e-7;                          % ps^4/m
alphay = 0;                             % 1/m

lamdapump1 = 1502.6;                    % nm
lamdapump2 = 1600.6;                    % nm
lamdazero = 1550.0;                     % nm
pump_power_per_pump_per_mode = 0.5;       % W

gainpoints = 200;                      % number of gain points to be plotted
input_signal_power = 1e-9;             % W

% CALCULATED QUANTITIES 

T = nt*dt;                              % FFT window size (ps)
t = ((1:nt)'-(nt+1)/2)*dt;              % vector of t values (ps)
w = 2*pi*[(0:nt/2-1),(-nt/2:-1)]'/T;    % vector of w values (rad/ps)
v = 1000*[(0:nt/2-1),(-nt/2:-1)]'/T;    % vector of v values (GHz)
vw = [c ./ (v + c/lamdazero)]';         % wavelengths (nm) %Note that LamdaZero is used as the simulation central frequecy
vs = fftshift(v);                       % swap halves for plotting 

betapx = [0 beta1x beta2x beta3x beta4x]';	    % polynomial beta coefs 
betapy = [0 beta1y beta2y beta3y beta4y]';	    % polynomial beta coefs


%STARTING FIELD
fprintf(1,'\nCreating Input Field... \n')

%------adding input pumps in time domain-----------------------------------------------
upump = sqrt(pump_power_per_pump_per_mode); 

signalandpumpfrequencies = c ./[lamdapump1 lamdapump2 lamdazero];

deltavi1 = signalandpumpfrequencies(3)  - signalandpumpfrequencies(1);
deltavi2 = signalandpumpfrequencies(3)-signalandpumpfrequencies (2);

[aux1a,aux1b] = min(abs(v - deltavi1));
[aux2a,aux2b] = min(abs(v - deltavi2));

u0x = upump* exp(-j*(pi/2))*exp(-j*2*pi*1e-3*v(aux1b)*t) + exp(-j*(pi/2))*upump*exp(-j*2*pi*1e-3*v(aux2b)*t);
u0y = upump*exp(-j*2*pi*1e-3*v(aux1b)*t) + upump*exp(-j*2*pi*1e-3*v(aux2b)*t);
 
%-------adding input signals at each half of the spectrum---------------------------
% Note that odd frequencies are chosen in half of the spectrum and even
% ones in the other half, so we can get the full gain-spectrum runing just
% one propagation

stepfreq = floor(nt/(gainpoints + 1)); % determines the indices of the frequencies used for plot
for ii = 1:gainpoints/2;
    u0x = u0x + sqrt(input_signal_power).*exp(j*2*pi*1e-3*vs(nt/2 + stepfreq*ii)*t);
    u0x = u0x + sqrt(input_signal_power).*exp(j*2*pi*1e-3*vs(nt/2 - stepfreq*ii + floor(stepfreq/2))*t);
    u0y = u0y + exp(-j*(pi/2))*sqrt(input_signal_power).*exp(j*2*pi*1e-3*vs(nt/2 + stepfreq*ii)*t);
    u0y = u0y + exp(-j*(pi/2))*sqrt(input_signal_power).*exp(j*2*pi*1e-3*vs(nt/2 - stepfreq*ii + floor(stepfreq/2))*t);
end

%--------------------------------------------------------------------------

%PROPAGATE
fprintf(1,'\n\nProapagation using the SSFM');
tx=cputime;
[ux,uy] = SSFM_for_CNLSE(u0x,u0y,dt,dz,nz,alphax,alphay,betapx,betapy,gamma);
fprintf(1,'\nSSFM completed in (%.2f s)\n',cputime-tx);
%--------------------------------------------------------------------------


%VISUALIZATION
fprintf(1,'\n\nVisualization...'); 

%----Plot Gain Using SSFM result---------------------------------------------
spectralpower0x = fftshift((abs(fft(u0x)).^2)/nt^2);     % Watts
spectralpower0y = fftshift((abs(fft(u0y)).^2)/nt^2);     % Watts
spectralpowerx = fftshift((abs(fft(ux)).^2)/nt^2);       % W
spectralpowery = fftshift((abs(fft(uy)).^2)/nt^2);       % W
spectralpower0 = spectralpower0x + spectralpower0y;      % W
spectralpower = spectralpowerx + spectralpowery;         % W

gains = spectralpower./spectralpower0;

selectedgains = [];
selectedfreq = [];

for ii = 1:gainpoints/2;
    selectedfreq = [selectedfreq, vs(nt/2 - stepfreq*ii + floor(stepfreq/2))];
    selectedgains = [selectedgains, gains(nt/2 - stepfreq*ii + floor(stepfreq/2))];
end
selectedfreq = sort(selectedfreq);
selectedgains = selectedgains(length(selectedgains):-1:1);

for ii = 1:gainpoints/2;
    selectedfreq = [selectedfreq, vs(nt/2 + stepfreq*ii)];
    selectedgains = [selectedgains, gains(nt/2 + stepfreq*ii)];
end


vw = [c ./ (selectedfreq + c/lamdazero)]';            % wavelengths for plot (nm)
Gains_in_dB = 10*log10(selectedgains);


% figure(1)
% plot(vw,10*log10(selectedgains),'.-' );
% grid on; 
% xlabel ('\lambda (nm)'); 
% ylabel ('Gain (dB)'); 
% title ('Gain Spectra');

%----Plot Analytical Result-----------------------------------------------
w1 = 2*pi*c /lamdapump1*1e-3;        % THz
w2 = 2*pi*c /lamdapump2*1e-3;        % THz
wc = (w1+w2)/2;                      % THz (central frequency betwen pumps)       
wz = 2*pi*c /lamdazero*1e-3;         % THz zero dispersion wavelength

wz_wc = wc-wz;                       % THz
beta2x = (beta4x/2) * (wz_wc)^2 + (beta3x) * (wz_wc);  %  ps^2/m % beta2 at each wavelength

selectedfreq = 2*pi*c./vw*1e-3;


DeltaBeta = beta2x*((selectedfreq-wc).^2 -((w1-w2)/2)^2) + (beta4x/12)*( (selectedfreq-wc).^4 -((w1-w2)/2)^4 );
P1= 2*pump_power_per_pump_per_mode;
P2 = 2*pump_power_per_pump_per_mode;
P0 = P1+P2;
L = nz*dz;
u = 2/3;
b = 0;
gmax = 2*sqrt(P1*P2)* gamma * b;
k = DeltaBeta + gamma*P0*u;
g = (gmax^2 - (k/2).^2).^0.5;
G = 1 + abs( (gmax./g) .* sinh(g*L)).^2;
G_Analytical_dB = 10*log10(G);

createfigure(vw,Gains_in_dB,G_Analytical_dB)

fprintf(1,'\n----------------------------------------------');
fprintf(1,'\n');
