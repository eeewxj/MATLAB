function [u1x,u1y] = SSFM_for_CNLSE(u0x,u0y,dt,dz,nz,alphax,alphay,betapx,betapy,gamma);

nt = length(u0x); % nt is the number of point
w = 2*pi*[(0:nt/2-1),(-nt/2:-1)]'/(dt*nt); %constructing used frequencies

linearoperatorx = -alphax/2; % in the lines above are constructed the linear operator in X and Y. Observe that for the 
for ii = 0:length(betapx)-1; % first propagation until correlation lenght beta0x=beta0y=0.
  linearoperatorx = linearoperatorx - j*betapx(ii+1)*(w).^ii/factorial(ii);
end

linearoperatory = -alphay/2;
for ii = 0:length(betapy)-1;
  linearoperatory = linearoperatory - j*betapy(ii+1)*(w).^ii/factorial(ii);
end

u1x = u0x;
ufftx = fft(u0x); 

u1y = u0y;
uffty = fft(u0y);
  
fiberlength = nz*dz;
propagedlength = 0;


fprintf(1, '\nSimulation running...      ');
for i=1:nz,
    
  %Performing the SSFM with the given dz
  %------linear propagation until dz/2------------
  halfstepx = exp(linearoperatorx*dz/2);
  uhalfx = ifft(halfstepx.*ufftx); 
  
  halfstepy = exp(linearoperatory*dz/2);
  uhalfy = ifft(halfstepy.*uffty); 
          
  uhalfxaux = uhalfx;
  uhalfx = (uhalfxaux + j*uhalfy)/sqrt(2); % transforming into a CP basis
  uhalfy = (uhalfxaux - j*uhalfy)/sqrt(2); % transforming into a CP basis

  
  u1x = uhalfx .* exp(-j*gamma*(2/3)*( (abs(uhalfx).^2) + (2*abs(uhalfy).^2) )*dz); % non linear propagation in CP basis
  u1y = uhalfy .* exp(-j*gamma*(2/3)*( (abs(uhalfy).^2) + (2*abs(uhalfx).^2) )*dz); % non linear propagation in CP basis
    
  u1xaux = u1x; 
  u1x = (u1xaux + u1y)/sqrt(2); % transforming into a LP basis
  u1y = -j*(u1xaux - u1y)/sqrt(2); % transforming into a LP basis
  
  %------linear propagation from dz/2 to dz------------
  ufftx = halfstepx.*fft(u1x); % 
  uffty = halfstepy.*fft(u1y);
  
  propagedlength = propagedlength + dz;
  propagatededlength = propagedlength ;
   fprintf(1, '\b\b\b\b\b\b%5.1f%%', propagedlength * 100.0 /fiberlength );
  
end

u1x = ifft(ufftx);
u1y = ifft(uffty);

