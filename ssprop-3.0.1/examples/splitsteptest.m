% This m-file is designed as a test of the splitstep routine.
% When it is working correctly, this should reproduce the results
% of Agrawal, Fig. 4.14.  Specifically, this routine simulates
% the case of pulse propagation at the zero-dispersion point, in
% the presence of nonlinearity.  This test case was chosen
% because both the output waveform and the output spectrum are
% asymmetric, which allows one to identify sign errors in the
% algorithms.

T = 100;                     % time window
nt = 2^12;                   % number of points
dt = T/nt;                   % timestep
t = ((1:nt)'-(nt+1)/2)*dt;   % time vector
w = wspace(T,nt);            % angular frequency vector
vs = fftshift(w/(2*pi));     % frequency (shifted for plotting)

z = 5;                       % total distance
nz = 2000;                   % number of steps
dz = z/nz;                   % stepsize
alpha = 0;                   % loss coefficient
betap = [0,0,0,1];           % dispersion polynomial
gamma = 1;                   % nonlinearity coefficient
u0 = exp(-t.^2);             % input pulse (Gaussian)

tic;
u = sspropc(u0,dt,dz,nz,alpha,betap,gamma);   % propagate signal
toc;

subplot(211);
plot (t,abs(u0).^2,t,abs(u).^2);
xlim([-4 16]);
grid on;
xlabel ('(t-\beta_1z)/T_0');
ylabel ('|u(z,t)|^2/P_0');
title ('Initial and Final Pulse Shapes');

subplot(212);
U0 = fftshift(abs(dt*fft(u0)/sqrt(2*pi)).^2); % input power spectrum
U = fftshift(abs(dt*fft(u)/sqrt(2*pi)).^2);   % output power spectrum
plot (vs,U0,vs,U);
xlim([-0.6 0.6]);
grid on;
xlabel ('(\nu-\nu_0) T_0');
ylabel ('|U(z,\nu)|^2/P_0');
title ('Initial and Final Pulse Spectra');
