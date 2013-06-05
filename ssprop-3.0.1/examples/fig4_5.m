T = 200;                                % time window (period)
nt = 2^14;                              % number of points
dt = T/nt;                              % timestep (dt)
t = ((1:nt)'-(nt+1)/2)*dt;              % time vector
w = wspace(T,nt);                       % angular frequency vector
vs = fftshift(w/(2*pi));                % shifted for plotting

z = 4.5*pi;                             % distance

subplot(121);

u0 = gaussian(t,0,2*sqrt(log(2)),1,1,10);

u = ssrklem(u0,dt,z/10,z,0,0,1);         % perform calculations
U = fftshift(abs(dt*ifft(u)*nt/sqrt(2*pi)).^2);

plot (vs,U);
xlim([-5 5]);
xlabel ('(\nu-\nu_0) T_0');
ylabel ('|U(z,\nu)|^2/P_0');
title ('SPM Pulse Spectrum, Gaussian Pulse (m=1,C=+10,\phi = 4.5\pi)');

subplot(122);

u0 = gaussian(t,0,2*sqrt(log(2)),1,1,-10);

u = ssrklem(u0,dt,z/10,z,0,0,1);         % perform calculations
U = fftshift(abs(dt*ifft(u)*nt/sqrt(2*pi)).^2);

plot (vs,U);
xlim([-6 6]);
xlabel ('(\nu-\nu_0) T_0');
ylabel ('|U(z,\nu)|^2/P_0');
title ('SPM Pulse Spectrum, Gaussian Pulse (m=1,C=-10,\phi = 4.5\pi)');
