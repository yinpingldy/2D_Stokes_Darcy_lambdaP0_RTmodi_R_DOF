clc
clear

syms x y
%% f_s
u1 = ;
u2 = ;
ps = ;

u1_xx = diff(u1,x,2);
u1_yy = diff(u1,y,2);
u2_xy = diff( diff(u2,x),y);
ps_x = diff(ps,x);

u1_yx = diff( diff(u1,y),x);
u2_xx = diff(u2,x,2);
u2_yy = diff(u2,y,2);
ps_y = diff(ps,y);
%% f_d
pd = ;
pd_x = diff(pd,x);
pd_y = diff(pd,y);

%% yita
u1_x = diff(u1,x);
u1_y = diff(u1,y);
u2_x = diff(u2,x);