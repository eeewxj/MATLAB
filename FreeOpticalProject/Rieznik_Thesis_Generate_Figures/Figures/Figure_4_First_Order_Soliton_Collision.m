clear all
format long e
fprintf(1,'\n----------------------------------------------');
fprintf(1,'\nSimulation of First Order Solitons Collision'); 

% PARAMETERS
A = 1;
ni = 0.44;
gamma = 2.2;
nt = 3072;
beta2 = -0.1;
dt = 400/nt;
t = -200:dt:(200-dt); %time in ps

u1 = (A*ni*sqrt(abs(beta2)/gamma)*sech(ni*t))';
u2 = (A*ni*sqrt(abs(beta2)/gamma)*sech(ni*(t-100)))';

u0 = u1 + u2.*exp(j*2*pi*1e-3*800*t');
 
alpha = 0;
betap = [0 0 beta2]';

% PROPAGATE
tx=cputime;
dz =10;
tol= (10^-5);
nz =  400/dz;
fprintf(1,'\n\nProapagation using the SSFM started');
[u,number_of_FFTs] = UPM_SSFM(u0,dt,dz,nz,alpha,betap,gamma,tol);
% [u,number_of_FFTs] = LEM_SSFM(u0,dt,dz,nz,alpha,betap,gamma,tol);
% [u,number_of_FFTs] = NLPR_SSFM(u0,dt,dz,nz,alpha,betap,gamma,tol);
fprintf(1,'\nSSFM completed in (%.2f s)\n',cputime-tx);

% Calculating global error
load analitical_collision;
globalerror = sqrt(sum(abs(analitical_collision-u).^2) / sum(abs(analitical_collision).^2))

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
title ('Initial (blue) and Final (red) Pulse Shapes');
