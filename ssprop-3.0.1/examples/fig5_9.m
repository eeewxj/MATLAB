T = 50;                                 % time window (period)
nt = 2^12;                              % number of points
dt = T/nt;                              % timestep (dt)
t = ((1:nt)'-(nt+1)/2)*dt;              % time vector
w = wspace(T,nt);                       % angular frequency vector
vs = fftshift(w/(2*pi));                % shifted for plotting

z = 2*(pi/2);                           % total distance
nz = 1000;                               % total number of steps
nplot = 30;	    						% number of plots to make
n1 = round(nz/nplot);					% number of steps per plot
nz = n1*nplot;							% total number of steps (revised)
dz = z/nz;                              % step-size

betap = [0,0,+1];						% dispersion polynomial

u0 = 3*tanh(t).*exp(1i*pi*t/T);

zv = (z/nplot)*(0:nplot);
u = zeros(length(t),length(zv));
U = zeros(length(t),length(zv));

u(:,1) = u0;
U(:,1) = fftshift(abs(dt*fft(u(:,1))/sqrt(2*pi)).^2);

for ii = 1:nplot,
  u(:,ii+1) = ssprop(u(:,ii),dt,dz,n1,0,betap,1);
  U(:,ii+1) = fftshift(abs(dt*fft(u(:,ii+1))/sqrt(2*pi)).^2);
end

mesh(zv/(pi/2),t,abs(u).^2, ...
     'MeshStyle', 'col', 'EdgeColor', 'black');
xlim([0 2]);
ylim([-15 15]);
xlabel ('Z/Z_0');
ylabel ('(t-\beta_1z)/T_0');
zlabel ('|u(z,t)|^2/P_0');
