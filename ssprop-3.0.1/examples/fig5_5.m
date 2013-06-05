T = 50;                                 % time window (period)
nt = 2^12;                              % number of points
dt = T/nt;                              % timestep (dt)
t = ((1:nt)'-(nt+1)/2)*dt;              % time vector
w = wspace(T,nt);                       % angular frequency vector
vs = fftshift(w/(2*pi));                % shifted for plotting

z = pi/2;                          		% total distance
nz = 500;                               % total number of steps
nplot = 10;	    						% number of plots to make
n1 = round(nz/nplot);					% number of steps per plot
nz = n1*nplot;							% total number of steps (revised)
dz = z/nz;                              % step-size

betap = [0,0,-1];						% dispersion polynomial

zv = (z/nplot)*(0:nplot);
u = zeros(length(t),length(zv));
U = zeros(length(t),length(zv));

u(:,1) = solitonpulse(t,0,1,3);         % input pulse
U(:,1) = fftshift(abs(dt*fft(u(:,1))/sqrt(2*pi)).^2);

for ii = 1:nplot,
  u(:,ii+1) = sspropc(u(:,ii),dt,dz,n1,0,betap,1);
  U(:,ii+1) = fftshift(abs(dt*fft(u(:,ii+1))/sqrt(2*pi)).^2);
end

mesh(zv/(pi/2),vs,U,...
     'MeshStyle', 'col', 'EdgeColor', 'black');
set(gca,'YDir','reverse');
hidden off;
xlim([0 1]);
ylim([-5 5]);
xlabel ('Z/Z_0');
ylabel ('(\nu-\nu_0) T_0');
zlabel ('|U(z,\nu)|^2/P_0');
