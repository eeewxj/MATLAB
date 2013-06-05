T = 200;
nt = 2^12;
dt = T/nt;
t = ((1:nt)'-(nt+1)/2)*dt;
w = wspace(T,nt);
vs = fftshift(w/(2*pi));
z = 3;
nz = 200;
dz = z/nz;

betap = [0,0,0,1];

u0 = solitonpulse(t,0,1);
U0 = fftshift(abs(dt*fft(u0)/sqrt(2*pi)).^2);

u = sspropc(u0,dt,dz,nz,0,betap,2^2);
U = fftshift(abs(dt*fft(u)/sqrt(2*pi)).^2);

subplot(121);
plot (t,abs(u0).^2,t,abs(u).^2);
xlim([-10 20]);
grid on;
xlabel ('(t-\beta_1z)/T_0');
ylabel ('|u(z,t)|^2/P_0');

subplot(122);
plot (vs,U0,vs,U);
xlim([-0.8 0.8]);
grid on;
xlabel ('(\nu-\nu_0) T_0');
ylabel ('|U(z,\nu)|^2/P_0');
