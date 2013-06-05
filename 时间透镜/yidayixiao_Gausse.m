%This code 考虑了高斯脉冲的传输并计算了其放大倍数
% This code solves the NLS equation with the split-step method
% idu/dz-(beta2/2)*d^2u/d(tau)^2+gama*|u|^2*u=0
%Written by Govind P.Agrawal in March 2005 for the NLFO book

close all;%
clear all;%
clc;
%--Specify input parameters
beta21=input('Enter Second dispersion value=');
% P0=input('Enter pulse peak value=');
T01=input('Enter pulse initial width value=');  %Gauss pulse half width
gama=input('Enter NL value=');

% Ld1=T01^2/abs(beta21);  %dispersion length
% Lnl=1/(gama*P0);  %nonlinear length

%--Specify input parameters
distance1=11;%
D1=beta21*distance1;
% mshape=input('m = 0 for sech,m > 0 for super-Gaussian =');
% chirp0=0;  %input pulse chirp(default value)

%---set simulation parameters
nt=4096; Tmax=100;   %FFT points and period size
num=5;          %num. of pulses
N=nt*num;        %total sampling points
step_num1=round(20*distance1*5);  %No. of z steps to
deltaz1=distance1/step_num1;  %step size in z
dtau=Tmax/nt;   % step size in tau

%---tau and omega arrays
tau=((1:N)'-(N+1)/2)*dtau;   % temporal grid 此处不能写成 tau=(-N/2:N/2-1)*dtau  不清楚什么原因
v=[(0:N/2-1),(-N/2:-1)]'/(Tmax*N);   %frequency  grid
vs=fftshift(v);       %swap halves for plotting
omega =2*pi*vs;       %Angular frequency grid

%--Input Field profile here only considering m>0
% if mshape == 0
%    uu = sech(tau).*exp(-0.5i*chirp0*tau.^2);  % soliton
% else    %super-Gaussian
base= exp(-0.5*tau.^2/T01^2);
shape=base;
% for pp=1:2:num-1
    shape= shape+0.5*circshift(base,nt);
%     shape=shape+base;
% end
% for pp=2:2:num-1
%     base= circshift(base,nt);
%     shape=shape+base;
% end
uu=shape;
    % end 
%---plus Gausse pulse

%---Plot input pulse shape and spectrume
temp = fftshift(fft(uu));  %spectrum
figure; subplot(2,2,1);
plot(tau,abs(uu).^2,'c');
hold on;
xlabel('Normalized Time');
ylabel('Normalized Power');
title('Input and Output Pulse Shape');
subplot(2,2,2);
plot(vs,abs(temp).^2,'c');
hold on;
xlabel('Normalized Frequency');
ylabel('Spectral Power');
title('Input and Output Pulse Spectrum');
%--store dispersive phase shifts to speedup code
dispersion1 = exp(1i*0.5*beta21*omega.^2*deltaz1);   %phase factor
hhz1 = 1i*gama*deltaz1; % nonlinear phase factor

%****************[Beginning of MAIN Loop]********************
%scheme:1/2N->D->1/2N;first half step nonlinear
temp = uu.*exp(abs(uu).^2.*hhz1/2);   %note hhz/2
% L1=length(temp)
% L2=length(dispersion1)
% L3=fftshift(fft(temp))
for n=1:step_num1
    f_temp = fftshift(fft(temp)).*dispersion1;
    uu = ifft(fftshift(f_temp));
    temp = uu.*exp(abs(uu).^2.*hhz1);
%     L4=length(temp)
end
uu = temp.*exp(-abs(uu).^2.*hhz1/2);  %Final field

temp=fftshift(fft(uu));  %Final spectrum
%***************[End of MAIN Loop]*************************
%---Plot output pulse shape and spectrum
subplot(2,2,1)
plot(tau, abs(uu).^2,'k');
subplot(2,2,2)
plot(vs,abs(temp).^2,'k');

%****************[Time Lens]******************
% A=50;Wm=10;   %frequency(GHz)
Phi = 0.5*0.005*(tau).^2;
uu1 = uu.*exp(-1i*Phi);
subplot(2,2,3);
plot(tau,abs(uu1).^2,'c')
xlabel('Normalized Time');
ylabel('Spectral Power');
title('After time-lens Pulse Shape');
temp1=fftshift(fft(uu1));
subplot(2,2,4);
plot(vs,abs(temp1).^2,'k');
xlabel('Normalized Frequency');
ylabel('Spectral Power');
title('After time-lens Pulse Spectrum');
%plot(tau,exp(Phi),'-k')

%*******************[pulse width after time lens]********************
% x1=find(abs(uu1).^2>=1/exp(1));     %find uu1 max. position
% x2=length(x1);       %find uu1 half-width position
% T02=x2/2*dtau;         %get uu1 FWHM


%******************[Second dispersion]****************
%input parameters
beta22 = -20;
distance2 = 110;
D2=beta22*distance2;
% Ld2=T02^2/abs(beta22);  %dispersion length

M=-D2/D1;    %magnification factor

step_num2=round(20*distance2*5);  
deltaz2=distance2/step_num2;

%--store dispersive phase shifts to speedup code
dispersion2 = exp(1i*0.5*beta22*omega.^2*deltaz2);   %phase factor
hhz2 = 1i*gama*deltaz2; % nonlinear phase factor

%****************[Beginning of MAIN Loop]********************
%scheme:1/2N->D->1/2N;first half step nonlinear
temp2 = uu1.*exp(abs(uu1).^2.*hhz2/2);   %note hhz/2
for n=1:step_num2
    f_temp = fftshift(fft(temp2)).*dispersion2;
    uu1 = ifft(fftshift(f_temp));
    temp2 = uu1.*exp(abs(uu1).^2.*hhz2);
end
uu2 = temp2.*exp(-abs(uu1).^2.*hhz2/2);  %Final field
temp2=fftshift(fft(uu2));  %Final spectrum
%***************[End of MAIN Loop]*************************
%---Plot output pulse shape and spectrum
figure;
subplot(2,1,1)
plot(tau.*M, abs(uu2).^2,'k')
xlabel('Normalized Time');
ylabel('Normalized Power');
title('Final Pulse Shape');
subplot(2,1,2)
plot(vs,abs(temp2).^2,'k')
xlabel('Normalized Frequency');
ylabel('Normalized Power');
title('Final Pulse Spectrum');

