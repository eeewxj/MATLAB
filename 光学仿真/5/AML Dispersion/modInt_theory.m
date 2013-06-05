function Eout = modInt_theory(Ein,m,fm)
% modInt_theory(Ein,m,fm)
% intensity modulator model

% modulator parameters
%   chirp alpha factor: alpha
%   extinction ratio: esilon (dB)
% modulation parameters
%   modulation index: m
%   modulation frequency: fm

global Ts;

N = size(Ein,1);
k = (1:N)';
Eout = Ein.*exp(-m/4*(2*pi*fm*Ts)^2*(k-N/2).*(k-N/2));



