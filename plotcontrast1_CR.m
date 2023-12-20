%% 绘图脚本
% 本代码最初是为了绘制算法的参数敏感度，大家可以简单借鉴。作者：蒋林志
% 输入数据
beta_all = [0.225317239,1.618056317,16.10197733,161.0094647,1610.091186;
    0.03275159,0.03287445,0.040408391,0.233027839,2.301260148];
sel = '1/16';
cost_plot(sel,beta_all,'R'); % 调用函数

%% 绘图函数
function cost_plot(sel,cur_data,parameter)
%% 维数选择
% MNIST to USPS
Dim = [0:size(cur_data,2)-1];% x轴大小

t = tiledlayout(2,2,'TileSpacing','Compact'); % 设置图例位置
set(figure(1),'Position',[200,200,1000,260]);
% 设置颜色
c1 = [0.254901960784314,0.411764705882353,0.882352941176471]; % 蓝色
c2 = [1,0,0]; % 红色

% line1
Line1 = line(Dim(1:1:5),log10(cur_data(1,1:5)));
set(Line1, 'LineStyle', '--','LineWidth', 1,  'Color', c1);
set(Line1, 'Marker', 'o','MarkerEdgeColor',c1,'MarkerFaceColor','w','MarkerSize',6);
% line2
Line2 = line(Dim(1:1:5),log10(cur_data(2,1:5)));
set(Line2, 'LineStyle', '--','LineWidth', 1,  'Color', c2);
set(Line2, 'Marker', '+','MarkerEdgeColor',c2,'MarkerFaceColor','w','MarkerSize',6);


set(gca, 'Box', 'off', ...                                % 边框设置
         'XGrid', 'off', 'YGrid', 'on', ...               % 网格线调整
         'TickDir', 'out', 'TickLength', [.01 .01], ...   % 刻度调整
         'XMinorTick', 'off', 'YMinorTick', 'off', ...    
         'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1],...
         'FontName', 'Times', 'FontSize', 9);  % 坐标轴颜色

set(gca,'XTick',Dim) %x轴范围1-8，间隔1
set(gca,'XTicklabel',{'1','10','10^2','10^3','10^4'});% 设置x轴刻度

% set(gca,'YTick',0:100:10000)
set(gca,'YTicklabel',{'10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}', '10^{3}', '10^{4}'});% 设置x轴刻度

legend('Pressure','Velocity','Location','NorthOutside',...
    'Orientation','horizon','FontName','Helvetica','FontSize',12,'FontWeight','normal'); %右上角标注
title(['Error variation when h=',sel],'FontName','Helvetica','FontSize',10,'FontWeight','normal') % 标题（可去掉）
xlabel(parameter, 'FontName','Helvetica','FontSize',10,'FontWeight','normal')  %x轴坐标描述
ylabel('error', 'FontName','Helvetica','FontSize',10,'FontWeight','normal') %y轴坐标描述

end