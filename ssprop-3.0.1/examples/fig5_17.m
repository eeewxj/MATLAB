T = 50;              					% FFT window size / rep rate
nt = 2^10;                              % number of points in FFT
dt = T/nt;                              % normalized time step
z = 5*(pi/2);                      		% total distance
nz = 1000;                              % total number of steps
nplot = 5;								% number of plots to make
n1 = round(nz/nplot);					% number of steps per plot
nz = n1*nplot;							% total number of steps (revised)
dz = z/nz;                              % step-size

t = ((1:nt)'-(nt+1)/2)*dt;              % vector of t values

s = 0.2;								% toptical/(2*pi*T0)
betap = [0,0,-1];                       % dispersion

zv = (z/nplot)*(0:nplot);
u = zeros(length(t),length(zv));

u(:,1) = solitonpulse(t,0,1,2);

for ii = 1:nplot,
  u(:,ii+1) = sspropc(u(:,ii),dt,dz,n1,0,betap,1,0,2*pi*s);
end

mesh(zv/(pi/2),t,abs(u).^2, ...
     'MeshStyle', 'col', 'EdgeColor', 'black');
set(gca,'YDir','reverse');
hidden off;
xlim([0 5]);
ylim([-16 16]);
xlabel ('Z/Z_0');
ylabel ('(t-\beta_1z)/T_0');
zlabel ('|u(z,t)|^2/P_0');
