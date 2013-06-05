clear all
format long e
fprintf(1,'\n----------------------------------------------');
fprintf(1,'\nSimulation of Second-Order Soliton Propagation'); 

% PARAMETERS

A = 2;
ni = 0.44; % ps^-1
gamma = 2.2; % W^-1 * km^-1
duration = 1/ni; %ps
nt = 2^10; % = 1024
beta2 = -0.1; % ps^2/m
time = 50; % ps
dt = time/nt; % ps
t = -time/2:dt:(time/2-dt); % ps

LD = duration^2/-beta2; 

u0 = (A*ni*sqrt(abs(beta2)/gamma)*sech(ni*t))'; %initial field shape in W^0.5
 
alpha = 0;
betap = [0 0 beta2]';


% PR0PAGATE finding numerical solution
tx=cputime;
dz = 10;
tol= 10^-7; %change here to change the achieved accuracy
nz =  8.113617390469509e+001/dz;
fprintf(1,'\n\nProapagation using the SSFM started');
[u,number_of_FFTs] = UPM_SSFM(u0,dt,dz,nz,alpha,betap,gamma,tol);
% [u,number_of_FFTs] = LEM_SSFM(u0,dt,dz,nz,alpha,betap,gamma,tol);
% [u,number_of_FFTs] = NLPR_SSFM(u0,dt,dz,nz,alpha,betap,gamma,tol);
fprintf(1,'\nSSFM completed in (%.2f s)\n',cputime-tx);



% Calculating analitical solution
LD = ni^-2/-beta2;
period=pi*ni^-2/(2*-beta2);
length = period;
qsi = length/LD;
tao = t/ni^-1;

num = 4*(cosh(3*tao) + 3*exp(4*j*qsi)*cosh(tao)).*exp(j*qsi/2);
den = cosh(4*tao)+4*cosh(2*tao)+3*cos(4*qsi);
u_analitic = num./den;
u_analitic = u_analitic./(gamma*LD)^0.5;


% Calculating global error
globalerror = sqrt(sum(abs(u_analitic-u').^2) / sum(abs(u_analitic).^2))

% giving the number of FFTs performed
number_of_FFTs

fprintf(1,'----------------------------------------------');
fprintf(1,'\n');

% PLOT OUTPUT
figure(1);
plot (t,abs(u0).^2,'*-',t,abs(u).^2,'o-');
grid on;
xlabel ('t (ps)');
ylabel ('|u(z,t)|^2 (W)');
title ('Initial (blue) and Final (green) Pulse Shapes');




