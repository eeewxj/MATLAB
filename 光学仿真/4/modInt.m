function Eout = modInt(Ein,alpha,epsilon,m,fm,bias)
% modInt(Ein,alpha,epsilon,m,fm)
% Amplitude modulator model

% modulator parameters
%   chirp alpha factor: alpha
%   extinction ratio: esilon (dB)

% modulation parameters
%   modulation index: m
%   modulation frequency: fm
%   bias: 0:1

global Ts;

N = size(Ein,1);
k = (1:N)';
Vm = m*cos(2*pi*fm*Ts*(k-N/2));
ext = 1 - 4/pi*atan(1/sqrt(10^(epsilon/10)));
delta_phi = pi/4*(2-bias*2-ext*Vm);
Eout = Ein.*cos(delta_phi).*exp(1i*alpha*delta_phi);




