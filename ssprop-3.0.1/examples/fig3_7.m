clear;close all;clc;
T = 200;                                % time window (period)
nt = 2^12;                              % number of points
dt = T/nt;                              % timestep (dt)
t = ((1:nt)'-(nt+1)/2)*dt;              % time vector

z = 6;                                  % total distance
nz = 50;                                 % total number of steps
dz = z/nz;                         		% total distance per step

betap = [0,0,0,1];						% dispersion polynomial

zv = (dz)*(0:nz);
u = zeros(length(t),nz+1);

u(:,1) = gaussian(t,0,2*sqrt(log(2)),1,3);
tic;
for ii = 1:nz,
  u(:,ii+1) = ssprop(u(:,ii),dt,dz,1,0,betap,0);
end
toc;
mesh(zv,t,abs(u).^2,...
     'MeshStyle', 'col', 'EdgeColor', 'black');
set(gca,'YDir', 'reverse');
hidden off;
xlim([0 6]);
ylim([-20 20]);
xlabel ('Z/L_d');
ylabel ('(t-\beta_1z)/T_0');
zlabel ('|u(z,t)|^2/P_0');
