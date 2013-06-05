clear all;close all;clc;
a1 = imageprocess('1.png',1020,1060,-73.2,6.8,68,433);
a2 = imageprocess('2.png',1020,1060,-72.4,7.6,68,433);
a3 = imageprocess('3.png',1020,1060,-71.6,8.4,68,433);
a4 = imageprocess('4.png',1020,1060,-68.8,11.2,68,433);
a5 = imageprocess('5.png',1010,1070,-73.0,7.0,68,433);
a6 = imageprocess('6.png',1020,1060,-70.8,9.2,68,433);
a7 = imageprocess('7.png',1010,1070,-73.7,6.3,68,433);
a8 = imageprocess('8.png',1010,1070,-73.1,6.9,68,433);

figure;hold on;
set(gca, 'FontName', 'Times New Roman','Box','on','linewidth',2);
plot(a1(:,1),a1(:,2),'b--');
plot(a2(:,1),a2(:,2),'g-.');
plot(a3(:,1),a3(:,2),'r:');
plot(a4(:,1),a4(:,2),'k');

xlabel('wavelength/nm');
ylabel('power/dBm');
legend('P_p_u_m_p=10.1dBm','P_p_u_m_p=10.5dBm',...
    'P_p_u_m_p=10.9dBm','P_p_u_m_p=12.6dBm');

figure;hold on;
set(gca, 'FontName', 'Times New Roman','Box','on','linewidth',2);
plot(a5(:,1),a5(:,2),'b--');
plot(a6(:,1),a6(:,2),'g-.');
plot(a7(:,1),a7(:,2),'r:');
plot(a8(:,1),a8(:,2),'k');

xlabel('wavelength/nm');
ylabel('power/dBm');
legend('\lambda_p_u_m_p=1041.2nm','\lambda_p_u_m_p=1042.2nm',...
    '\lambda_p_u_m_p=1044.2nm','\lambda_p_u_m_p=1045.4nm');

