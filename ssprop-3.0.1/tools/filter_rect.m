function Eout = filter_rect(Ein,f3dB)
% filter_rect(Ein,f3dB)
% filter the input signal with a rect function
global Ts;

N = size(Ein,1);
% the k element in the Ein corresponds to the frequency of (k-N/2)/Ts/N
k = wspace(2*pi,N);
w = k/Ts/N;

Eout = Ein.*(abs(w)<(f3dB/2));



