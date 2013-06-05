T = 30;                                 % time window (period)
nt = 2^10;                              % number of points
dt = T/nt;                              % timestep (dt)
t = ((1:nt)'-(nt+1)/2)*dt;              % time vector

z = 90;                      		    % total distance
nz = 2000;                              % total number of steps
nplot = 45;	    						% number of plots to make
n1 = round(nz/nplot);					% number of steps per plot
nz = n1*nplot;							% total number of steps (revised)
dz = z/nz;                              % step-size

betap = [0,0,-1];						% dispersion polynomial
q0 = 3.5;                               % separation between pulses

subplot(221);

theta = 0;
r = 1;
u0 = solitonpulse(t,-q0) + r*solitonpulse(t,+q0)*exp(i*theta);

zv = (z/nplot)*(0:nplot);
u = zeros(length(t),length(zv));

u(:,1) = u0;

for ii = 1:nplot,
  u(:,ii+1) = sspropc(u(:,ii),dt,dz,n1,0,betap,1);
end

mesh(zv,t,abs(u).^2, ...
     'MeshStyle', 'col', 'EdgeColor', 'black');
xlim([0 90]);
ylim([-8 8]);
set(gca,'YDir', 'reverse');
view(-75,60);
xlabel ('Z/Z_0');
ylabel ('(t-\beta_1z)/T_0');
zlabel ('|u(z,t)|^2/P_0');

subplot(222);

theta = pi/4;
r = 1;
u0 = solitonpulse(t,-q0) + r*solitonpulse(t,+q0)*exp(i*theta);

zv = (z/nplot)*(0:nplot);
u = zeros(length(t),length(zv));

u(:,1) = u0;

for ii = 1:nplot,
  u(:,ii+1) = sspropc(u(:,ii),dt,dz,n1,0,betap,1);
end

mesh(zv,t,abs(u).^2, ...
     'MeshStyle', 'col', 'EdgeColor', 'black');
xlim([0 90]);
ylim([-8 8]);
set(gca,'YDir', 'reverse');
view(-75,60);
xlabel ('Z/Z_0');
ylabel ('(t-\beta_1z)/T_0');
zlabel ('|u(z,t)|^2/P_0');

subplot(223);

theta = pi/2;
r = 1;
u0 = solitonpulse(t,-q0) + r*solitonpulse(t,+q0)*exp(i*theta);

zv = (z/nplot)*(0:nplot);
u = zeros(length(t),length(zv));

u(:,1) = u0;

for ii = 1:nplot,
  u(:,ii+1) = sspropc(u(:,ii),dt,dz,n1,0,betap,1);
end

mesh(zv,t,abs(u).^2, ...
     'MeshStyle', 'col', 'EdgeColor', 'black');
xlim([0 90]);
ylim([-8 8]);
set(gca,'YDir', 'reverse');
view(-75,60);
xlabel ('Z/Z_0');
ylabel ('(t-\beta_1z)/T_0');
zlabel ('|u(z,t)|^2/P_0');

subplot(224);

theta = 0;
r = 1.1;
u0 = solitonpulse(t,-q0) + r*solitonpulse(t,+q0)*exp(i*theta);

zv = (z/nplot)*(0:nplot);
u = zeros(length(t),length(zv));

u(:,1) = u0;

for ii = 1:nplot,
  u(:,ii+1) = sspropc(u(:,ii),dt,dz,n1,0,betap,1);
end

mesh(zv,t,abs(u).^2, ...
     'MeshStyle', 'col', 'EdgeColor', 'black');
xlim([0 90]);
ylim([-8 8]);
set(gca,'YDir', 'reverse');
view(-75,60);
xlabel ('Z/Z_0');
ylabel ('(t-\beta_1z)/T_0');
zlabel ('|u(z,t)|^2/P_0');
