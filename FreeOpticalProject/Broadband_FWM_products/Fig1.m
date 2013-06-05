clear all
format long e

fprintf(1,'\nStart "Simulation of broadband FWM cascade"'); 

% CONSTANTS

c = 299792458;                          % speed of light (m/s)

% NUMERICAL PARAMETERS

nt = 2^15;      % number of points in FFT
dt = 200/nt;    % time step (ps)
L = 3;          % fiber length (m) CHANGE HERE TO GENERATE THE SPECTRA IN FIG1a!!

% OPTICAL PARAMETERS
gamma = 10/1000; % W^-1 m^-1
S0 = 0.075/1000; % ps/nm^2/m    CHANGE HERE TO GENERATE THE SPECTRA IN FIG1b!!

lamdapump1 = 1555; % nm
lamdapump2 = 1563; % nm
lamdazero = 1550;  % nm
FrequencyZeroForSimulartion=(c/lamdapump1+c/lamdapump2)/2; % Ghz
LamdaZeroForSimulartion=c/FrequencyZeroForSimulartion;     % Ghz

beta1 = 0;                                   % ps/m                                      
beta2 = 0;                                   % ps^2/m
beta3 =  S0*lamdazero^4/(4*pi^2*(c/1000)^2); % ps^3/m   
beta4 = 1e-7;                                % ps^4/m
alpha = 0;                                   % 1/m

DeltaW = 2*pi*c*(lamdazero-LamdaZeroForSimulartion)/lamdazero/LamdaZeroForSimulartion; % GHz

beta2 = beta3*DeltaW/1000; % ps^2/m

% CALCULATED QUANTITIES

T = nt*dt;                              % FFT window size (ps)
t = ((1:nt)'-(nt+1)/2)*dt;              % vector of t values (ps)
w = 2*pi*[(0:nt/2-1),(-nt/2:-1)]'/T;    % vector of w values (rad/ps)
v = 1000*[(0:nt/2-1),(-nt/2:-1)]'/T;    % vector of v values (GHz)
vs = fftshift(v);                       % swap halves for plotting 
vw = [c ./ (v + c/LamdaZeroForSimulartion)]';  % wavelengths (nm) 

betap = [0 beta1 beta2 beta3 beta4]';	    % polynomial beta coefs 

%STARTING FIELD
fprintf(1,'\n\n\nConstructing input field'); 

%------pumps------------------------------------------------------------
power1 = 100; % W
power2 = 100; % W
up1 = (power1*ones(nt,1)).^0.5; % W^0.5
up2 = (power2*ones(nt,1)).^0.5; % W^0.5

%--------noise----------------------------------------------------------

noise_powerx= 1e-6;% in W (change here to change the noise power)
noise_resx= 12.5;% in GHz (change here to change the bandwidth at which the noise power is given)
noise = randn(length(t),1).*(noise_powerx*1000/dt/noise_resx).^0.5.*exp(i*2*pi*rand(length(t),1));


%-------Construct input field: noise and pumps-------------------------------------------

signalandpumpfrequencies = c ./[lamdapump1 lamdapump2 LamdaZeroForSimulartion];

deltavi1 = signalandpumpfrequencies(3)  - signalandpumpfrequencies(1);
deltavi2 = signalandpumpfrequencies(3)-signalandpumpfrequencies (2);

[aux1a,aux1b] = min(abs(v - deltavi1));
[aux2a,aux2b] = min(abs(v - deltavi2));

u0 = noise + up1.*exp(-j*2*pi*1e-3*v(aux1b)*t) + up2.*exp(-j*2*pi*1e-3*(v(aux1b)+ 1000)*t);  

%----------------------------------------------------------------------

%PROPAGATE
fprintf(1,'\n\nProapagation using the SSFM');    
tx=cputime;
u = UPM(u0,dt,L,1,alpha,betap,gamma,1e-4);
fprintf(1,'\n\nSSFM completed in (%.2f s)\n',cputime-tx);


%----OSA output -----------------------------------
res =12.5; %RB in GHz
[slotfreq,powerinres] = resolution(u,vs,res,nt,dt);

figure(1)
plotwav(powerinres.^0.5,LamdaZeroForSimulartion,slotfreq,1000,LamdaZeroForSimulartion)