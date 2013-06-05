clear all;clc;
format long e
fprintf(1,'\n\n\n----------------------------------------------');
fprintf(1,'\nSimulating Supercontinuum Generation in PCF');


% Constants or Parameters *****************************************************
c = 299792.458;                         %speed of ligth nm/ps
tol = 1e-4;                             % photon error number, or local error, depending on the used method. 

% Input Field Paramenters
tfwhm = 30e-3;                          % ps
ni = 1/tfwhm;                           % ps^-1
lamda_central = 835;                    %nm
fo=c/lamda_central;                     % central pulse frequency (THz)


% Fiber Parameters
gamma = 0.110;                          % W^-1 * m^-1
alpha = 0;                              % atenuation coef. (m^-1)
L = 1e-0;                               % fiber length (m)
% beta coefficients (ps^n/m)
betaw = [0 0 -11.830 8.1038e-2 -9.5205e-5 2.0737e-7 -5.3943e-10 1.3486e-12 -2.5495e-15 3.0524e-18 -1.714e-21]*1e-3; 


% Numerical Parameters
nt = 2^15;                              % number of spectral points`
time = 2^6;                             % ps
dt = time/nt;                           % ps
t = -time/2:dt:(time/2-dt);             % ps
dz = 1e-4;                              % initial longitudinal step (m)
v = [(0:nt/2-1),(-nt/2:-1)]/(dt*nt);    % frequencies frequencies (THz)

% INPUT FIELD ***********************************************************
lamdaP = 700;
deltavP = c/lamdaP-c/lamda_central;     % THz
[~,aux2b] = min(abs(v - deltavP));
CarrierP = exp(-1j*2*pi*v(aux2b)*t);
PeakPower = 10000; % W, change here the Input Power!
u0 = sqrt(PeakPower)*sechpulse(t,0,tfwhm).*CarrierP; %initial field shape in W^0.5



% PR0PAGATE finding numerical solution **********************************
%************************************************************************
fprintf(1,'\n\nInteraction Picture Method started    ');

tic

nplot = 10;                              % must be times of 5
up = zeros(nt,nplot+1);
ufourier = zeros(nt,nplot+1); 
up(:,1) = u0.';
for ii = 1:nplot
    fprintf('\b\b\b\b%4.0f', ii);
    ufourier(:,ii) = (fftshift(abs(dt*(ifft(up(:,ii))*nt).^2/sqrt(2*pi))));   
    up(:,ii+1) = ssrklem(up(:,ii),dt,dz,L/nplot,alpha,betaw,gamma,1/fo,tol);
end
ufourier(:,nplot+1) = (fftshift(abs(dt*(ifft(up(:,nplot+1))*nt).^2/sqrt(2*pi))));
plotstep = nplot/5;
vw = (c./(fftshift(v) + fo));
figure(1);
hold on;
for i = 1:5
    fprintf('\b\b\b\b%4.0f', i);
    surfc((L/nplot)*(plotstep*(i-1):i*plotstep),(vw),(10*log10(ufourier(:,plotstep*(i-1)+1:i*plotstep+1)/(max(max(ufourier))))),...
        'MeshStyle', 'col', 'EdgeColor', 'none');
end
set(gca,'YDir','reverse');              %Y坐标轴反转
camup([1 0 0]);                         %设置观察者的仰视方向向量
campos([L/2 lamda_central 0]);          %设置观察者的位置
camtarget([L/2 lamda_central -20]);     %设置观察目标的位置
grid off; 
set(gca,'TickDir','out');               %设置刻度线的朝向
hidden off;                             %设置隐藏线消除模式，即是否用虚线显示被遮挡的部分
set(gca,'CLim',[-40 0]);
xlim([0 L]);
ylim([500 1400]);
xlabel ('z/m','Rotation',90,'Position',[L/2 1550 0]);
ylabel ('\lambda/nm','Rotation',0,'Position',[-L/10 800 0]);
u = up(:,nplot+1);

tx = toc;

fprintf(1, '\n\nSimulation lasted (s) = ');
fprintf(1, '%5.2f%', tx );

%  ---------plot output spectrum------------------
specnorm = fftshift(abs(ifft(u)).^2);
specnorm = specnorm/max(specnorm);
figure(2)
hold on
plot(c./(fftshift(v) + fo),10*log10(specnorm));
grid on;
xlabel ('lambda (nm)');
ylabel ('Normalized Spectrum (a.u.)');
title ('Output Spectrum');
axis([500 1400 -70 1])
fprintf(1,'\n----------------------------------------------\n\n');
