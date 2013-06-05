function hR =ramanhR(t)
% This function computes the raman response function
%
% USAGE:
% 
% u = ramanhR (t);t unit /ps

tau1 = 12.2e-3;     %ps
tau2 = 32e-3;       %ps      
taub = 96e-3;       %ps      
fb = 0.21;
tres = t-t(1);  
ha = ((tau1^2+tau2^2)/(tau1*tau2^2))*exp(-tres/tau2).*sin(tres/tau1);
hb=(2*taub-tres)/taub^2 .* exp(-tres/taub);
hR=(1-fb)*ha + fb * hb;  

end