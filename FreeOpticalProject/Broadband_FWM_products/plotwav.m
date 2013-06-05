function  f = plotwav(u,lambda,v,range,vwcenterplot)

% This function graphs u^2 (in dBm) vs wavelenght for a given 
% value of wavelenght for the zero spectral component
% of u (lambda). The frequencies and v are required.
%
% INPUTS
% 
% u - input field (W^0.5)
% lambda - wavelength of the zero frequency (nm)
% v - frequencies at which u is given (GHz)
% range - range for plot (nm)
% vwcenterplot - center wavelength for plot (nm)
%
% OUTPUTS
% Plot wavelength vs u^2 (nm vs dBm)
% Exported files of input and output signal 

c = 299792458;                          % speed of light (m/s) 
vw = [c ./ (v + c/lambda)]';            % wavelengths (nm) 
nt = length(u);

aux = vw - vwcenterplot; 
iv = find( abs(aux) < range); 
plot(vw(iv),10*log10( 1000 * u(iv).^2),'.-' ); 
grid on; 
xlabel ('\lambda (nm)'); 
ylabel ('|U(z,\nu)|^2 (dBm)'); 
title ('Spectra');