T = 60;              					% FFT window size / rep rate
nt = 2^12;                              % number of points in FFT
dt = T/nt;                              % normalized time step
z = 0.2;                           		% total distance
nz = 1000;                              % total number of steps
dz = z/nz;                              % step-size

t = ((1:nt)'-(nt+1)/2)*dt;              % vector of t values
w = wspace(t);                          % vector of w values
vs = fftshift(w/(2*pi));                % used for plotting

s = 0.01;								% toptical/(2*pi*T0)
betap = [0,0,1];                        % dispersion 

u0 = gaussian(t,0,2*(sqrt(log(2))));

u1 = sspropc(u0,dt,dz,nz,0,betap,10^2,0,2*pi*s);
U1 = fftshift(abs(dt*ifft(u1)*nt/sqrt(2*pi)).^2);

u2 = ssprop(u1,dt,dz,nz,0,betap,10^2,0,2*pi*s);
U2 = fftshift(abs(dt*ifft(u2)*nt/sqrt(2*pi)).^2);

subplot(221);
plot (t,abs(u1).^2);
xlim([-6 6]);
xlabel ('(t-\beta_1z)/T_0');
ylabel ('|u(z,t)|^2/P_0');

subplot(222);
plot (vs,U1);
xlim([-4 4]);
xlabel ('(t-\beta_1z)/T_0');
ylabel ('|u(z,t)|^2/P_0');

subplot(223);
plot (t,abs(u2).^2);
xlim([-8 8]);
xlabel ('(\nu-\nu_0) T_0');
ylabel ('|U(z,\nu)|^2/P_0');

subplot(224);
plot (vs,U2);
xlim([-4 4]);
xlabel ('(\nu-\nu_0) T_0');
ylabel ('|U(z,\nu)|^2/P_0');
