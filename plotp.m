clc
clear

t=0:0.02:1;
[x,y]=meshgrid(t); 

z=(12.*x.^2.*exp(y)).*(x<0.5)+(16.*x.*y.^3-exp(1)-2).*(x>0.5); 
% z=(2.* (y-1).* (cos(x)).^2).*(x<0.5)+(y.* (cos(y).^2)+4.*x-5/2).*(x>0.5); 

i = find(z==0);    %找出值为0的点
z(i)=NaN;  %将值为0的点赋值NaN

surf(x,y,z)

xlabel('x-axis');
ylabel('y-axis');
zlabel('pressure');

% subplot(2,2,1)
% 
% surf(xx,yy,z);title('Surfplot'); %子图1，绘制三维图形
% 
% subplot(2,2,2)
% 
% mesh(xx,yy,z);title('Meshplot'); %子图2，绘制三维曲面
% 
% subplot(2,2,3)
% 
% surf(xx,yy,z);title('Surplot with shading interp'); %子图3，绘制三维曲面，表面为光滑
% 
% shading interp;
% 
% subplot(2,2,4)
% 
% contour(xx,yy,z);title('Meshplot'); %子图4，绘制等高曲线