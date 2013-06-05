function Eout = filter_gaus(Ein,f3dB,n)
% filter_gaus(Ein,f3dB,n)
% 利用n阶高斯滤波器在频域滤波，注意输入的Ein应为频域形式
% T(f) = exp(-log(sqrt(2))*(2/f3dB/Ts/N)^2n*(k.^2n));

global Ts;
N = size(Ein,1);
% the k element in the Ein corresponds to the frequency of (k-N/2)/Ts/N
k = wspace(2*pi,N);
% Eout = Ein.*exp(-log(sqrt(2))*(2/f3dB/Ts/N)^n*(k.^(2*n)));
% n order of gaussian filter
temp = log(sqrt(2))*(2./f3dB./Ts./N)^(2*n);
Eout = Ein.*exp(-temp*(k.^(2*n)));      % n order gaussian filter VPI



