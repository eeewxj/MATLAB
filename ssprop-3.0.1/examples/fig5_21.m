clc;
T = 100;                                % time window (period)
nt = 2^12;                              % number of points
dt = T/nt;                              % timestep (dt)
t = ((1:nt)'-(nt+1)/2)*dt;              % time vector
w = wspace(T,nt);                       % angular frequency vector
vs = fftshift(w/(2*pi));                % shifted for plotting

z = 3*(pi/2);                      		% total distance
nz = 5000;                              % total number of steps
nplot = 3;								% number of plots to make
n1 = round(nz/nplot);					% number of steps per plot
nz = n1*nplot;							% total number of steps (revised)
dz = z/nz;                              % step-size

beta2 = -1;
beta3 = .03;
betap = [0,0,beta2,beta3];				% dispersion polynomial
s = 0.05;
tr = 0.1;
N = 2;

u0 = solitonpulse(t);

zv = (z/nplot)*(0:nplot);
u = zeros(length(t),length(zv));
U = zeros(length(t),length(zv));

u(:,1) = u0;
U(:,1) = fftshift(abs(dt*fft(u(:,1))/sqrt(2*pi)).^2);

for ii = 1:nplot,
  u(:,ii+1) = sspropc(u(:,ii),dt,dz,n1,0,betap,N^2,tr,2*pi*s);
  U(:,ii+1) = fftshift(abs(dt*fft(u(:,ii+1))/sqrt(2*pi)).^2);
end

subplot(211);
mesh(zv/(pi/2),t,abs(u).^2, ...
     'MeshStyle', 'col', 'EdgeColor', 'black');
set(gca,'YDir', 'reverse');
hidden off;
xlim([0 3]);
ylim([-25 25]);
xlabel ('Z/Z_0');
ylabel ('(t-\beta_1z)/T_0');
zlabel ('|u(z,t)|^2/P_0');

subplot(212);
mesh(zv/(pi/2),vs,U, ...
     'MeshStyle', 'col', 'EdgeColor', 'black');
set(gca,'YDir', 'reverse');
hidden off;
xlim([0 3]);
ylim([-2 2]);
xlabel ('Z/L_d');
ylabel ('(\nu-\nu_0) T_0');
zlabel ('|U(z,\nu)|^2/P_0');
