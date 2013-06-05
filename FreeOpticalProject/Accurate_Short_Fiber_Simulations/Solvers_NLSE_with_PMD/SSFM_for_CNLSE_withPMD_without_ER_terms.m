function [u1x,u1y] = UPMPol(u0x,u0y,dt,dz,nz,alphax,alphay,betapx,betapy,gamma,correlationlength,angles,tetas,beta0s,beta1s);

nt = length(u0x);
w = 2*pi*[(0:nt/2-1),(-nt/2:-1)]'/(dt*nt);

linearoperatorx = -alphax/2;
for ii = 0:length(betapx)-1;
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
propagedlength =0;
auxlength = 0;

jj=1;

fprintf(1, '\nSimulation running...      ');
while propagedlength < fiberlength,
      % the next "if" and the variable "auxlenght" are used to check
      % if the correlation length was riched
      if (auxlength >= correlationlength | auxlength ==0),
      auxlength = 0;
      

      angle = angles(jj);  
      teta = tetas(jj);
      ufftxaux = ufftx * cos(teta) - uffty * sin(teta)*exp(-i*angle);
      ufftyaux = ufftx * sin(teta)*exp(i*angle) + uffty * cos(teta);
      ufftx = ufftxaux;
      uffty = ufftyaux;
            
      betapx(1) = -beta0s(jj); % observe that b0x and b0y are defined such that b0x = -b0y = -beta0s(1)
      betapy(1) = beta0s(jj);
      betapx(2) = -beta1s(jj); % tha same as above for beta1x and beta1y
      betapy(2) = beta1s(jj);
      
      jj=jj+1; % I increase the value of jj so in the next time the correlation length is reached, different pseudo-random variables are taken
      
      linearoperatorx = -alphax/2; % in the lines bellow I define the linear operators again using the new beta0s and beta1s
      for ii = 0:length(betapx)-1;
        linearoperatorx = linearoperatorx - j*betapx(ii+1)*(w).^ii/factorial(ii);
      end
      linearoperatory = -alphay/2;
      for ii = 0:length(betapy)-1;
        linearoperatory = linearoperatory - j*betapy(ii+1)*(w).^ii/factorial(ii);
      end  
  end % here ends the "if" part that is done just at each correlation length.
  
   %Performing the SSFM with the given dz
  %------linear propagation until dz/2------------
  halfstepx = exp(linearoperatorx*dz/2);
  uhalfx = ifft(halfstepx.*ufftx); 
  
  halfstepy = exp(linearoperatory*dz/2);
  uhalfy = ifft(halfstepy.*uffty); 
          
   %------non linear propagation in dz------------
  u1x = uhalfx .* exp(-j*gamma*( (abs(uhalfx).^2) + ((2/3)*abs(uhalfy).^2) )*dz);
  u1y = uhalfy .* exp(-j*gamma*( (abs(uhalfy).^2) + ((2/3)*abs(uhalfx).^2) )*dz);
  
  %------linear propagation from dz until dz------------
  ufftx = halfstepx.*fft(u1x);
  uffty = halfstepy.*fft(u1y);
  %-----------------------------------------------------
  auxlength = auxlength + dz;
  propagedlength = propagedlength + dz;
  fprintf(1, '\b\b\b\b\b\b%5.1f%%', propagedlength * 100.0 /fiberlength );
end

u1x = ifft(ufftx);
u1y = ifft(uffty);

