%% 绘图脚本
% 本代码最初是为了绘制算法的参数敏感度，大家可以简单借鉴。作者：蒋林志
% 输入数据
beta_all = [0.133460191,0.065774866,0.03275159,0.016356342,0.008175684;
            0.142557501,0.066783498,0.03287445,0.016371683,0.008177614;
            0.480453105,0.119914123,0.040408391,0.017358269,0.008302896;
            4.569047658,0.985802434,0.233027839,0.058819601,0.016263084;
            45.62185218,9.816792749,2.301260148,0.563569148,0.140393226];
sel = 'error';
cost_plot(sel,beta_all,'h'); % 调用函数

%% 绘图函数
function cost_plot(sel,cur_data,parameter)
%% 维数选择
% MNIST to USPS
Dim = [0:size(cur_data,2)-1];% x轴大小

t = tiledlayout(2,2,'TileSpacing','Compact'); % 设置图例位置
set(figure(1),'Position',[200,200,1000,260]);
% 设置颜色
c1 = [0,0,0]; % 黑色
c2 = [1,0,0]; % 红色
c3 = [0.105882352941176,0.619607843137255,0.466666666666667];% 绿色
c4 = [0.545098039215686,0,0.545098039215686];% 紫色
c5 = [0.254901960784314,0.411764705882353,0.882352941176471]; % 蓝色

% line1
Line1 = line(Dim(1:1:5),log10(cur_data(1,1:5)));
set(Line1, 'LineStyle', '--','LineWidth', 1,  'Color', c1);
set(Line1, 'Marker', 'o','MarkerEdgeColor',c1,'MarkerFaceColor','w','MarkerSize',6);
% line2
Line2 = line(Dim(1:1:5),log10(cur_data(2,1:5)));
set(Line2, 'LineStyle', '--','LineWidth', 1,  'Color', c2);
set(Line2, 'Marker', '+','MarkerEdgeColor',c2,'MarkerFaceColor','w','MarkerSize',6);
% line3
Line3 = line(Dim(1:1:5),log10(cur_data(3,1:5)));
set(Line3, 'LineStyle', '--','LineWidth', 1,  'Color', c3);
set(Line3, 'Marker', '^','MarkerEdgeColor',c3,'MarkerFaceColor','w','MarkerSize',6);
% line4
Line4 = line(Dim(1:1:5),log10(cur_data(4,1:5)));
set(Line4, 'LineStyle', '--','LineWidth', 1,  'Color', c4);
set(Line4, 'Marker', '*','MarkerEdgeColor',c4,'MarkerFaceColor','w','MarkerSize',6);
% line5
Line5 = line(Dim(1:1:5),log10(cur_data(5,1:5)));
set(Line5, 'LineStyle', '--','LineWidth', 1,  'Color', c5);
set(Line5, 'Marker', 'd','MarkerEdgeColor',c5,'MarkerFaceColor','w','MarkerSize',6);


set(gca, 'Box', 'off', ...                                % 边框设置
         'XGrid', 'off', 'YGrid', 'on', ...               % 网格线调整
         'TickDir', 'out', 'TickLength', [.01 .01], ...   % 刻度调整
         'XMinorTick', 'off', 'YMinorTick', 'off', ...    
         'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1],...
         'FontName', 'Times', 'FontSize', 12);  % 坐标轴颜色

set(gca,'XTick',Dim) %x轴范围1-8，间隔1
set(gca,'XTicklabel',{'1/4','1/8','1/16','1/32','1/64'});% 设置x轴刻度

% set(gca,'YTick',0:100:10000)
set(gca,'YTicklabel',{'10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}', '10^{3}', '10^{4}'});% 设置x轴刻度

legend('R=1','R=10','R=10^2','R=10^3','R=10^4','Location','NorthOutside',...
    'Orientation','horizon','FontName','Helvetica','FontSize',12,'FontWeight','normal'); %右上角标注
title(['Velocity ',sel],'FontName','Helvetica','FontSize',16,'FontWeight','normal') % 标题（可去掉）
xlabel(parameter, 'FontName','Helvetica','FontSize',15,'FontWeight','normal')  %x轴坐标描述
ylabel('e_{0}(u)', 'FontName','Helvetica','FontSize',15,'FontWeight','normal') %y轴坐标描述

end