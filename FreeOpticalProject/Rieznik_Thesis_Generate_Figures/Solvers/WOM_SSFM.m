function [u1,number_of_FFTs] = ssprop(u0,dt,dz,nz,alpha,betap,gamma,tol);


if (nargin<8)
  tol = 1e-3;
end

propagedlength = 0;
fiberlength = dz*nz;

nt = length(u0);
w = 2*pi*[(0:nt/2-1),(-nt/2:-1)]'/(dt*nt);

halfstep = -alpha/2;
for ii = 0:length(betap)-1;
  halfstep = halfstep - j*betap(ii+1)*(w).^ii/factorial(ii);
end
halfstep = exp(halfstep*dz/2);

nf = 0;
ufft = fft(u0); nf = nf + 1;

fprintf(1, '\nSimulation running...      ');
for iz = 1:nz,
  uhalf = ifft(halfstep.*ufft); nf = nf + 1;
  uv = uhalf .* exp(-j*gamma*(abs(uhalf).^2)*dz);
  ufft = halfstep.*fft(uv); nf = nf + 1;
  propagedlength = propagedlength + dz;
  fprintf(1, '\b\b\b\b\b\b%5.1f%%', propagedlength * 100.0 /fiberlength );
end

if propagedlength < fiberlength,
  dz = fiberlength - propagedlength;
  uhalf = ifft(halfstep.*ufft); nf = nf + 1;
  uv = uhalf .* exp(-j*gamma*(abs(uhalf).^2)*dz);
  ufft = halfstep.*fft(uv); nf = nf + 1;
  propagedlength = propagedlength + dz;
  fprintf(1, '\b\b\b\b\b\b%5.1f%%', propagedlength * 100.0 /fiberlength );
end

u1 = ifft(ufft);
number_of_FFTs = nf;

% load propaged2;
% globalerror = sqrt(sum(abs(propaged2-ux).^2) / sum(abs(propaged2).^2))
% u1=ux;
