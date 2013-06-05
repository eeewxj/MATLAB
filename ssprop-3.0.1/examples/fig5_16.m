T = 16;              					% FFT window size / rep rate
nt = 2^10;                              % number of points in FFT
dt = T/nt;                              % normalized time step
z = 10;                           		% total distance
nz = 1000;                              % total number of steps
nplot = 2;								% number of plots to make
n1 = round(nz/nplot);					% number of steps per plot
nz = n1*nplot;							% total number of steps (revised)
dz = z/nz;                              % step-size

t = ((1:nt)'-(nt+1)/2)*dt;              % vector of t values

s = 0.2;								% toptical/(2*pi*T0)
betap = [0,0,-1];                       % dispersion

zv = (z/nplot)*(0:nplot);
u = zeros(length(t),length(zv));

u(:,1) = solitonpulse(t);

for ii = 1:nplot,
  u(:,ii+1) = sspropc(u(:,ii),dt,dz,n1,0,betap,1,0,2*pi*s);
end

plot (t,abs(u).^2);
xlim([-5 5]);
grid on;
xlabel ('(t-\beta_1z)/T_0');
ylabel ('|u(z,t)|^2/P_0');
