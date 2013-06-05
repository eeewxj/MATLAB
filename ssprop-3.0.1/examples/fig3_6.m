clear all;close all;clc;
T = 48;                                 % time window (period)
nt = 2^10;                              % number of points
dt = T/nt;                              % timestep (dt)
t = ((1:nt)'-(nt+1)/2)*dt;              % time vector

dz = 1;                            		% total distance per step
nz = 5;
u0 = gaussian(t,0,2*sqrt(log(2)));

betap = [0,0,0,1];						% dispersion polynomial
u1 = ssrklem(u0,dt,dz,nz,0,betap,0);
%u1 = ssfm(u0,dt,dz,nz,0,betap,0,3e-15);

betap = [0,0,1,1];						% dispersion polynomial
u2 = ssrklem(u0,dt,dz,nz,0,betap,0);
%u2 = ssfm(u0,dt,dz,nz,0,betap,0,3e-15);

plot(t,abs(u0).^2,':',...
     t,abs(u1).^2,'-',...
     t,abs(u2).^2,'--');
xlim([-12 12]);
xlabel ('(t-\beta_1z)/T_0');
ylabel ('|u(z,t)|^2/P_0');
legend('original pulse','\beta_2=0,\beta_3=1','\beta_2=1,\beta_3=1');
