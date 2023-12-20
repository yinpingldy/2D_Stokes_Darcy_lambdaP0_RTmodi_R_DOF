clc
clear

% [x,y] = meshgrid(0:0.05:1,0:0.05:1);
% if x<=1/2
%     u = exp(x.*y);
%     v = exp(x).*exp(2.*y);
% else
%     u = exp(x.*y);
%     v = 3.*x+y.^2;
% end

% if x<=1/2
%     u = sin(y.^2 + 6.*x);
%     v = cos(4.*(x.^2).*y);
% else
%     u = sin(y.^2 + 6.*x);
%     v = sin(2.*x).*cos(3.*y);
% end

[x,y] = meshgrid(0.01:0.05:0.5,0.01:0.05:1);
u = exp(x.*y);
v = exp(x).*exp(2.*y);


figure
quiver(x,y,u,v)
hold on

[x,y] = meshgrid(0.51:0.05:1,0.01:0.05:1);
u = exp(x.*y);
v = 3.*x+y.^2;
quiver(x,y,u,v)
axis([0,1,0,1]);
grid on
set(gca,'xtick',0:0.1:1);
set(gca,'ytick',0:0.1:1);
ax = gca;
ax.GridColor = [0 0 0];
ax.GridLineStyle = '--';
ax.GridAlpha = 0.5;
ax.Layer = 'top';

hold off