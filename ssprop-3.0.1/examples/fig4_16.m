T = 16;              					% FFT window size / rep rate
nt = 2^10;                              % number of points in FFT
dt = T/nt;                              % normalized time step
z = 20;                           		% total distance
nz = 1000;                              % total number of steps
nplot = 2;								% number of plots to make
n1 = round(nz/nplot);					% number of steps per plot
nz = n1*nplot;							% total number of steps (revised)
dz = z/nz;                              % step-size

t = ((1:nt)'-(nt+1)/2)*dt;              % vector of t values

s = 0.01;								% toptical/(2*pi*T0)

zv = (z/nplot)*(0:nplot);
u = zeros(length(t),length(zv));
U = zeros(length(t),length(zv));

u(:,1) = gaussian(t,0,2*(sqrt(log(2))));
U(:,1) = fftshift(abs(dt*ifft(u(:,1))*nt/sqrt(2*pi)).^2);

for ii = 1:nplot,
  u(:,ii+1) = sspropc(u(:,ii),dt,dz,n1,0,0,1,0,2*pi*s);
  U(:,ii+1) = fftshift(abs(dt*ifft(u(:,ii+1))*nt/sqrt(2*pi)).^2);
end

plot (t,abs(u).^2);
xlim([-3 3]);
grid on;
xlabel ('(t-\beta_1z)/T_0');
ylabel ('|u(z,t)|^2/P_0');
