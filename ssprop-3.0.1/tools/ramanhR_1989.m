function hR =ramanhR_1989(t)
% This function computes the raman response function
%
% USAGE:
% 
% u = ramanhR (t);t unit /ps

tau1 = 12.2e-3;     %ps
tau2 = 32e-3;       %ps      
   
tres = t-t(1);  
hR = ((tau1^2+tau2^2)/(tau1*tau2^2))*exp(-tres/tau2).*sin(tres/tau1);
end