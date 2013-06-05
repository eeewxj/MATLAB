T = 100;                                % time window (period)
nt = 2^12;                              % number of points
dt = T/nt;                              % timestep (dt)
t = ((1:nt)'-(nt+1)/2)*dt;              % time vector

z = 18;                         		% total distance
nz = 500;                               % total number of steps
nplot = 9;	    						% number of plots to make
n1 = round(nz/nplot);					% number of steps per plot
nz = n1*nplot;							% total number of steps (revised)
dz = z/nz;                              % step-size

betap = [0,0,-1];						% dispersion polynomial

zv = (z/nplot)*(0:nplot);
u = zeros(length(t),length(zv));
U = zeros(length(t),length(zv));

u(:,1) = sechpulse(t,0,2*acosh(sqrt(2)),1,-0.5);

for ii = 1:nplot,
  u(:,ii+1) = ssprop(u(:,ii),dt,dz,n1,0,betap,1);
end

h = mesh(zv,t,abs(u).^2,...
         'MeshStyle', 'col', 'EdgeColor', 'black');
set(gca,'YDir','reverse');
hidden off;
xlim([0 18]);
ylim([-16 16]);
xlabel ('Z/Z_0');
ylabel ('(t-\beta_1z)/T_0');
zlabel ('|u(z,t)|^2/P_0');
