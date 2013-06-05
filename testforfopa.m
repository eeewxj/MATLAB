clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1);
% surfc(L/nplot*(0:nplot),t,(10*log10(abs(up)/(max(max(abs(up)))))),...
%         'MeshStyle', 'col', 'EdgeColor', 'none');
%     
% set(gca,'YDir','reverse');                 %Y坐标轴反转
% camup([1 0 0]);                            %设置观察者的仰视方向向量
% campos([L/2 0 5]);                         %设置观察者的位置
% camtarget([L/2 0 -20]);                    %设置观察目标的位置
% grid off; 
% set(gca,'TickDir','out');                  %设置刻度线的朝向
% hidden off;                                %设置隐藏线消除模式，即是否用虚线显示被遮挡的部分
% set(gca,'CLim',[-40 0]);
% xlim([0 L]);
% ylim([-5 5]);
% xlabel ('z/m','Rotation',90,'Position',[L/2 6 0]);
% ylabel ('t/ps','Rotation',0,'Position',[-L/10 0 0]);
%%% make movie of pulse's spectrum of stft
% mov(1:nplot+1) = struct('cdata', [],...
%      'colormap', []);
% 
% framelength = 2^8;  
% d = framelength/4;
% k = framelength - d;
% fftnt = floor((nt-d)/k);
% win = hamming(framelength);
% for j=1:nplot+1
%  splitufft = zeros(framelength,fftnt);
%  for n = 1:fftnt
%     splitufft(:,n) = fftshift(abs((ifft(up((k*(n-1)+1:k*(n-1)+framelength),j).*win)...
%         *framelength).^2/sqrt(2*pi)));
%  end
%  surfc(t((0:fftnt-1)*k+framelength/2),vw((0.5:framelength-0.5)*(nt/framelength)),(10*log10(abs(splitufft)/(max(max(abs(splitufft)))))),...
%          'MeshStyle', 'col', 'EdgeColor', 'none');
%  set(gca,'YDir','reverse');                 %Y坐标轴反转
%  camup([1 0 0]);                            %设置观察者的仰视方向向量
%  campos([0 lamda_central 5]);                %设置观察者的位置
%  camtarget([0 lamda_central -20]);          %设置观察目标的位置
%  grid off; 
%  set(gca,'TickDir','out');                  %设置刻度线的朝向
%  hidden off;                                %设置隐藏线消除模式，即是否用虚线显示被遮挡的部分
%  set(gca,'CLim',[-40 0]);
%  xlim([-5,30]);
%  ylim([500 1400]);
%  
%  mov(:,j)=getframe; 
% end
% movie2avi(mov,'pulse_spectrual.avi','compression', 'None');


%%%% make movie of pulse's shape
% maxup = max(max(abs(up)));
% mov(1:nplot+1) = struct('cdata', [],...
%                         'colormap', []);
% hold off;
% for j=1:nplot+1
% plot(t,abs(up(:,j))/maxup);
% xlim([-5 40]);
% ylim([0,1])
% mov(:,j)=getframe; 
% end
% movie2avi(mov,'pulse_temporal.avi','compression', 'None');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure;
% plot(vw,S_dBm,'b');
% axis([vw(1) vw(N) -100 -30]);
% xlabel('Wavelength (nm)');
% ylabel('Intensity (dBm)');
% figure;
% L = dz*nz;
% lamda_central = vw(N/2);
% surfc((L/nplot)*(0:nplot),vw,(10*log10(ufourier/(max(max(ufourier))))),...
%         'MeshStyle', 'col', 'EdgeColor', 'none');
% 
% set(gca,'YDir','reverse');                      %Y坐标轴反转
% camup([1 0 0]);                                 %设置观察者的仰视方向向量
% campos([L/2 lamda_central 5]);                         %设置观察者的位置
% camtarget([L/2 lamda_central -20]);                    %设置观察目标的位置
% grid off; 
% set(gca,'TickDir','out');                       %设置刻度线的朝向
% hidden off;                                     %设置隐藏线消除模式，即是否用虚线显示被遮挡的部分
% set(gca,'CLim',[-60 0]);
% xlabel ('z/m','Rotation',90,'Position',[L/2 vw(N) 0]);
% ylabel ('\lambda/nm','Rotation',0,'Position',[-L/10 lamda_central 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spectraOUT1 = flipud((abs(fftshift(ifft(sechpulse(t,0,1e-12))*N)).^2)/N^2);	% W
% S_dBm1 = 10*log10(spectraOUT1);
% 
% figure;
% plot(vw,S_dBm1,'b');
% axis([vw(1) vw(N) -100 10]);
% xlabel('Wavelength (nm)');
% ylabel('Intensity (dBm)');


% figure(1);
% plot(vp,Gains_in_dB,'r');
% axis([vp(1) vp(gainpoints) -30 100]);
% xlabel('Wavelength (nm)');
% ylabel('Gain (dB)');
% title_string = ['\lambda_p=',num2str(lamdaP*1e9,'%5.1fnm'),...
%     ',L=',num2str(dz*nz,'%3.0fm'),...
%     ',P_p_e_a_k=',num2str(Ip,'%2.0fW'),...
%     ',tol=',num2str(tol,'%1.0e')];
% title (title_string);
% figure(2);
% plot(t,abs(Eout).^2)
% title (title_string);
% figure(3);
% plot(t,abs(Ein).^2)
% title (title_string);
% figure(1);
% plot(vw,S_dBm,'b');
% axis([vw(1) vw(N) -100 0]);
% xlabel('Wavelength (nm)');
% ylabel('Intensity (dBm)');
% title_string = ['\lambda_p=',num2str(lamdaP*1e9,'%5.1fnm'),...
%     ',L=',num2str(dz*nz,'%3.0fm'),...
%     ',P_p_e_a_k=',num2str(Ip,'%2.0fW'),...
%     ',tol=',num2str(tol,'%1.0e')];
% title (title_string);
% figure(2);
% plot(t,abs(Eout).^2)
% title (title_string);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
n = 2^7;                   % number of samples in a time period (/ps)
num = 2^8;                  % observed time period number (ps)
N = n*num;                  % Total sampling points
Ts = 1/n;                   % sampling distance (ps)
t = ((1:N)'-(N+1)/2)*Ts;    % vector of t values (ps)
v = [(0:N/2-1),(-N/2:-1)]'/(Ts*N);          % vector of v values (THz)
vs = fftshift(v);                           % swap halves for plotting 
u = sechpulse(t,0,15,1,8.95);
uu = abs(u).^2;
uufft = abs(fftshift(ifft(u))).^2;
%plot(vs,uufft);
disp(fwhm(t,uu));
disp(fwhm(vs,uufft));


