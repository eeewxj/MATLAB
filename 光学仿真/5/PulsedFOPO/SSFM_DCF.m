function Eout = SSFM_DCF(Eo,dz,nz)

global Ts;              % s
gamma = 1.4e-3;         % W^-1 m^-1
beta1x = 0;             % ps/m
beta2x =  2.13e-1;      % ps^2/m
beta3x = -1.27e-3;      % ps^3/m
beta4x =  5.26e-6;      % ps^4/m
betap = [0 beta1x beta2x beta3x beta4x]';	    % @ 1550nm
alpha = 1e-4;           % 1/m

N = length(Eo);         % the number of point
w = 2*pi*[(0:N/2-1),(-N/2:-1)]'/(Ts*N*1e12);         %constructing used frequencies (rad.THz)

linearoperatorx = -alpha/2;
for ii = 0:length(betap)-1;
  linearoperatorx = linearoperatorx - 1j*betap(ii+1)*(w).^ii/factorial(ii);     % (rad/m)
end
halfstepx = exp(linearoperatorx*dz/2);

ufft = fft(Eo); 
for i=1:nz,
  %------linear propagation until dz/2------------
  uhalf = ifft(halfstepx.*ufft);
  u1 = uhalf.*exp(-1j*gamma*(abs(uhalf).^2)*dz); % non linear propagation in CP basis
  %------linear propagation from dz/2 to dz------------
  ufft = halfstepx.*fft(u1); % 
end
Eout = ifft(ufft);
