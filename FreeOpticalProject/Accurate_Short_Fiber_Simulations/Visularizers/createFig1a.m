function createfigure(x1, y1, y2)
%CREATEFIGURE(X1,Y1,Y2)
%  X1:  vector of x data
%  Y1:  vector of y data
%  Y2:  vector of y data
 
%  Auto-generated by MATLAB on 18-Oct-2007 10:54:52
 
%% Create figure
figure1 = figure('FileName','C:\Documents and Settings\pessoal\Meus documentos\ssf\IncludingPolarization\Short_Fiber_OPAs\Article\Definitivos\Fig1a.fig');
 
%% Create axes
axes1 = axes(...
  'FontSize',12,...
  'FontWeight','bold',...
  'XGrid','on',...
  'XTick',[0 40 80 120 160 200],...
  'YGrid','on',...
  'Parent',figure1);
xlabel(axes1,'fiber length (m)');
ylabel(axes1,'angle rotation (rad)');
hold(axes1,'all');
 
%% Create plot
plot1 = plot(x1,y1,'LineWidth',3);
 
%% Create plot
plot2 = plot(...
  x1,y2,...
  'Color',[0.8314 0.8157 0.7843],...
  'LineWidth',3);
 
