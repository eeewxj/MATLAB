function [u1,number_of_FFTs] = LEM_SSFM(u0,dt,dz,nz,alpha,betap,gamma,tol)

% This function solves the nonlinear Schrodinger equation for
% pulse propagation in an optical fiber using the split-step
% Fourier method described in: 
% 
% Oleg V. Sinkin, Ronald Holzlner, John Zweck, and Curtis Menyuk, "Optimization
% of the Split-Step Fourier Method in Modeling Optical-Fiber Communications Systems",
% IEEE J. of Lightwave Technol. 21, 61-68 (2003).
%
% USAGE
%
% u1 = ssprop(u0,dt,dz,nz,alpha,betap,gamma);
% u1 = ssprop(u0,dt,dz,nz,alpha,betap,gamma,tol);
%
% INPUT
%
% u0 - starting field amplitude (vector)
% dt - time step
% dz - initial propagation stepsize
% nz - number of steps to take, ie, ztotal = dz*nz
% alpha - power loss coefficient, ie, P=P0*exp(-alpha*z)
% betap - dispersion polynomial coefs, [beta_0 ... beta_m]
% gamma - nonlinearity coefficient
% maxiter - max number of iterations (default = 4)
% tol - local error method tolerance (default = 1e-5)
%
% OUTPUT
%
% u1 - field at the output
% number_of_FFTs - number of Fast Fourier Transforms performed during the
% propagation
%
% NOTES  The dimensions of the input and output quantities can
% be anything, as long as they are self consistent.  E.g., if
% |u|^2 has dimensions of Watts and dz has dimensions of
% meters, then gamma should be specified in W^-1*m^-1.
% Similarly, if dt is given in picoseconds, and dz is given in
% meters, then beta(n) should have dimensions of ps^(n-1)/m.


if (nargin<8)
  tol = 1e-3;
end

nt = length(u0);
w = 2*pi*[(0:nt/2-1),(-nt/2:-1)]'/(dt*nt);
fiberlength = dz*nz;

propagedlength = 0;
nf = 0; % parameter to save the number of FFTs

% constructing linear operator
LinearOperator = -alpha/2;
for ii = 0:length(betap)-1;
  LinearOperator = LinearOperator - 1i*betap(ii+1)*(w).^ii/factorial(ii); 
end

ufft = fft(u0); nf = nf +1;

% Performig the SSFM according to the LEM spatial-step size
fprintf(1, '\nSimulation running...      ');
while (propagedlength < fiberlength),
  if (dz + propagedlength) > fiberlength,
      dz = fiberlength - propagedlength;
  end
  halfstep = exp(LinearOperator*dz/2);
  quarterstep=exp(LinearOperator*dz/4);
  uhalf = ifft(halfstep.*ufft); nf = nf +1;
  uquarter=ifft(quarterstep.*ufft); nf = nf +1;
 
  uc = uhalf .* exp(-1i*gamma*(abs(uhalf.^2 )*dz));
  ufftc = halfstep.*fft(uc); nf = nf +1;
  uc = ifft(ufftc); nf = nf +1;
    
  uf = uquarter.*exp(-1i*gamma*(abs(uquarter).^2 )*dz/2);
  ufftf = quarterstep.*fft(uf); nf = nf +1;
  uf = ifft(ufftf); nf = nf +1;
  uquarter=ifft(quarterstep.*ufftf); nf = nf +1;
  uf2 = uquarter.*exp(-1i*gamma*(abs(uquarter).^2 )*dz/2);
  ufftf = quarterstep.*fft(uf2); nf = nf +1;
  uf = ifft(ufftf); nf = nf +1;
    
  delta = sqrt(sum((abs(uf-uc)).^2))/sqrt(sum((abs(uf)).^2));
        
  if (delta < (tol/2))
      u1 = (4/3)*uf-(1/3)*uc;
      ufft = fft(u1); nf = nf +1;
      propagedlength = propagedlength + dz;
      dz = dz*(2^(1/3));
     
  elseif ( tol <= delta <= (2*tol))
     u1 = (4/3)*uf-(1/3)*uc;
     ufft = fft(u1); nf = nf +1;
     propagedlength = propagedlength + dz;
     dz = dz/(2^(1/3));
    
 elseif delta > (2*tol)
    dz = dz/2;
    halfstep = exp(LinearOperator*dz/2);
    quarterstep=exp(LinearOperator*dz/4);
    uhalf = ifft(halfstep.*ufft); nf = nf +1;
    uquarter=ifft(quarterstep.*ufft); nf = nf +1;
 else
    u1 = (4/3)*uf-(1/3)*uc;
    ufft = fft(u1); nf = nf +1;
    propagedlength = propagedlength + dz;
  end
  fprintf(1, '\b\b\b\b\b\b%5.1f%%', propagedlength * 100.0 /fiberlength );
end

% giving output parameters
u1 = ifft(ufft);
number_of_FFTs = nf;