clear all
format long e
fprintf(1,'\n----------------------------------------------');
fprintf(1,'\nSimulation of 8-channels WDM System'); 

% CONSTANTS

c = 299792458;                          % speed of light (m/s)

% NUMERICAL PARAMETERS

bits = 2^5;                             % bits number
nperbit = 2^11;
nt = bits*nperbit;                      % number of points in FFT
dt = bits*100/nt;                       % time step (ps)
dz = 10;                              % distance stepsize (m)
nz = 10000/dz;                              % number of z-steps
tol = 1e-3;                             % error tolerance

% OPTICAL PARAMETERS

beta2 = -6e-006;            % ps^2/m
gamma = 2/1000;                            % W^-1 m^-1
alpha = 0;

% CALCULATED QUANTITIES

T = nt*dt;                              % FFT window size (ps)
t = ((1:nt)'-(nt+1)/2)*dt;              % vector of t values (ps)
w = 2*pi*[(0:nt/2-1),(-nt/2:-1)]'/T;    % vector of w values (rad/ps)
v = 1000*[(0:nt/2-1),(-nt/2:-1)]'/T;    % vector of v values (GHz)
vs = fftshift(v);                       % swap halves for plotting 
betap = [0 0 0 0.1/1000 -4e-7]';	% polynomial beta coefs


load input_8_channels;
u0 = input_8_channels;
  
%------------------PROPAG--------------------------------
%PROPAGATE

tx=cputime;
tol = 10^-0; %change here to change the achieved accuracy
fprintf(1,'\n\nProapagation using the SSFM started');
[u,number_of_FFTs] = UPM_SSFM(u0,dt,dz,nz,alpha,betap,gamma,tol);
% [u,number_of_FFTs] = LEM_SSFM(u0,dt,dz,nz,alpha,betap,gamma,tol);
% [u,number_of_FFTs] = WOM_SSFM(u0,dt,dz,nz,alpha,betap,gamma,tol);
fprintf(1,'\nSSFM completed in (%.2f s)\n',cputime-tx);

% giving the number of FFTs performed
number_of_FFTs
% Calculating global error
load analitical_8_channels_10km;
globalerror = sqrt(sum(abs(analitical_8_channels_10km-u).^2) / sum(abs(analitical_8_channels_10km).^2))

fprintf(1,'----------------------------------------------');
fprintf(1,'\n');