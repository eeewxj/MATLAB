T =200;              					% FFT window size / rep rate
nt = 2^12;                              % number of points in FFT
dt = T/nt;                              % normalized time step
z = 5*(pi/2);                      		% total distance
nz = 500;                               % total number of steps
dz = z/nz;                              % step-size

t = ((1:nt)'-(nt+1)/2)*dt;              % vector of t values
w = wspace(t);                          % vector of w values
vs = fftshift(w/(2*pi));                % used for plotting

tr = 0.01;								% toptical/(2*pi*T0)
betap = [0,0,-1];                       % dispersion

u0 = solitonpulse(t,0,1,2);
U0 = fftshift(abs(dt*fft(u0)/sqrt(2*pi)).^2);

u = sspropc(u0,dt,dz,nz,0,betap,1,tr);
U = fftshift(abs(dt*fft(u)/sqrt(2*pi)).^2);

plot (vs,U0,vs,U);
xlim([-2 2]);
grid on;
xlabel ('(\nu-\nu_0) T_0');
ylabel ('|U(z,\nu)|^2/P_0');
