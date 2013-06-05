clear all
format long e
fprintf(1,'\n----------------------------------------------');
fprintf(1,'\nSimulating Figure 16a of Rieznik´s thesis'); 

% CONSTANTS

c = 299792458;                          % speed of light (m/s)

% NUMERICAL PARAMETERS

nt = 2^15;                              % number of points in FFT
dt = 4e-2;                              % time step (ps)
dz = 0.5;                                 % distance stepsize (m)
nz = 400;                               % number of z-steps

% OPTICAL PARAMETERS

gamma = 10/1000;                        % W^-1 m^-1

beta1x = 0;                             % ps/m                                      
beta2x = 0;                             % ps^2/m
beta3x =  0.1/1000 ;                    % ps^3/m   
beta4x = 1e-7;                          % ps^4/m
alphax = 0;                             % 1/m

beta1y = 0;                             % ps/m   
beta2y = 0;                             % ps^2/m
beta3y =  0.1/1000;                     % ps^3/m
beta4y = 1e-7;                          % ps^4/m
alphay = 0;                             % 1/m


% CALCULATED QUANTITIES

T = nt*dt;                              % FFT window size (ps)
t = ((1:nt)'-(nt+1)/2)*dt;              % vector of t values (ps)
w = 2*pi*[(0:nt/2-1),(-nt/2:-1)]'/T;    % vector of w values (rad/ps)
v = 1000*[(0:nt/2-1),(-nt/2:-1)]'/T;    % vector of v values (GHz)
vs = fftshift(v);                       % swap halves for plotting 

betapx = [0 beta1x beta2x beta3x beta4x]';	    % polynomial beta coefs 
betapy = [0 beta1y beta2y beta3y beta4y]';	    % polynomial beta coefs

%STARTING FIELD
fprintf(1,'\nCreating Input Field.. \n');

%------pump------------------------------------------------------------

up = sqrt(0.5)*ones(length(t),1);


%-------Construct input field: noise and pumps-------------------------------------------
lamdapump1 = 1502.6;
lamdapump2 = 1600.6;
lamdazero = 1550;

signalandpumpfrequencies = c ./[lamdapump1 lamdapump2 lamdazero];

deltavi1 = signalandpumpfrequencies(3)  - signalandpumpfrequencies(1);
deltavi2 = signalandpumpfrequencies(3)-signalandpumpfrequencies (2);

[aux1a,aux1b] = min(abs(v - deltavi1));
[aux2a,aux2b] = min(abs(v - deltavi2));

u0x =  exp(-i*pi/4)*up.*exp(-j*2*pi*1e-3*v(aux1b)*t) + exp(i*pi/4)*up.*exp(-j*2*pi*1e-3*v(aux2b)*t); 
u0y =  exp(i*pi/4)*up.*exp(-j*2*pi*1e-3*v(aux1b)*t) + exp(-i*pi/4) * up.*exp(-j*2*pi*1e-3*v(aux2b)*t);

%-------Construct input field: signals at each half of the spectrum---------------------------

gainpoints = 150;
inputsignalpower = 1e-10; % W
stepfreq = floor(nt/(gainpoints + 1));
for ii = 1:gainpoints/2;
    u0x = u0x + exp(-i*pi/4)*sqrt(inputsignalpower).*exp(j*2*pi*1e-3*vs(nt/2 + stepfreq*ii)*t);
end
for iii = 1:gainpoints/2;
    u0x = u0x + exp(-i*pi/4)*sqrt(inputsignalpower).*exp(j*2*pi*1e-3*vs(nt/2 - stepfreq*iii + floor(stepfreq/2))*t);
end

for ii = 1:gainpoints/2;
    u0y = u0y + exp(i*pi/4)*sqrt(inputsignalpower).*exp(j*2*pi*1e-3*vs(nt/2 + stepfreq*ii)*t);
end
for iii = 1:gainpoints/2;
    u0y = u0y + exp(i*pi/4)*sqrt(inputsignalpower).*exp(j*2*pi*1e-3*vs(nt/2 - stepfreq*iii + floor(stepfreq/2))*t);
end
%----------Fixed Fiber random variables generation---------------------
correlationlength = 50;% tetas = [0]; beta0s = [0]; angles = [0];


load tetas400B
load beta0s400B

beta1s = 1000*beta0s/(2*pi*signalandpumpfrequencies(3));


correlationlength = dz;

angles = zeros(1,length(tetas));
 
%----------------------------------------------------------------------

%PROPAGATE
fprintf(1,'\n\nProapagation using the SSFM with ER terms started');
tx=cputime;
[ux,uy] = SSFM_for_CNLSE_withPMD(u0x,u0y,dt,dz,nz,alphax,alphay,betapx,betapy,gamma,correlationlength,angles,tetas,beta0s,beta1s);
fprintf(1,'\nSSFM completed in (%.2f s)\n',cputime-tx);

fprintf(1,'\n\nProapagation using the SSFM without ER terms started');
tx=cputime;
[ux2,uy2] = SSFM_for_CNLSE_withPMD_without_ER_terms(u0x,u0y,dt,dz,nz,alphax,alphay,betapx,betapy,gamma,correlationlength,angles,tetas,beta0s,beta1s);
fprintf(1,'\nSSFM completed in (%.2f s)\n',cputime-tx);

fprintf(1,'\n\nVisualization...'); 
%----Calculating gain with ER--------------------------------------------------------------------
spectralpower0x = fftshift((abs(fft(u0x)).^2)/nt^2);
spectralpower0y = fftshift((abs(fft(u0y)).^2)/nt^2);
spectralpowerx = fftshift((abs(fft(ux)).^2)/nt^2);
spectralpowery = fftshift((abs(fft(uy)).^2)/nt^2);
spectralpower0 = spectralpower0x + spectralpower0y;
spectralpower = spectralpowerx + spectralpowery;

gains = spectralpower./spectralpower0;

selectedgains = [];
selectedfreq = [];

for iii = 1:gainpoints/2;
    selectedfreq = [selectedfreq, vs(nt/2 - stepfreq*iii + floor(stepfreq/2))];
    selectedgains = [selectedgains, gains(nt/2 - stepfreq*iii + floor(stepfreq/2))];
end
selectedfreq = sort(selectedfreq);
selectedgains = selectedgains(length(selectedgains):-1:1);

for ii = 1:gainpoints/2;
    selectedfreq = [selectedfreq, vs(nt/2 + stepfreq*ii)];
    selectedgains = [selectedgains, gains(nt/2 + stepfreq*ii)];
end


%----Calculating gain without ER-----------------------------------------------------
spectralpowerx2 = fftshift((abs(fft(ux2)).^2)/nt^2);
spectralpowery2 = fftshift((abs(fft(uy2)).^2)/nt^2);
spectralpower2 = spectralpowerx2 + spectralpowery2;

gains2 = spectralpower2./spectralpower0;
selectedgains2 = [];
for iii = 1:gainpoints/2;
    selectedgains2 = [selectedgains2, gains2(nt/2 - stepfreq*iii + floor(stepfreq/2))];
end
selectedgains2 = selectedgains2(length(selectedgains2):-1:1);

for ii = 1:gainpoints/2;
    selectedgains2 = [selectedgains2, gains2(nt/2 + stepfreq*ii)];
end
%%%%%%%%%%   Vsualize   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
vw = [c ./ (selectedfreq + c/lamdazero)]';            % wavelengths (nm)
createFig3(vw,10*log10(selectedgains),10*log10(selectedgains2));

fprintf(1,'\n----------------------------------------------');
fprintf(1,'\n');