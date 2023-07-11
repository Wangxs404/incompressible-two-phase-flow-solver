clc
clear
global imin imax jmin jmax nx ny 
load('D2Data_Multy');%FVM

%% 绘制流线图
subplot(221)
x=linspace (0 ,Lx, nx+1);
y=linspace (0 ,Ly, ny+1);
uu=u(4:imax+1,4:imax+1)';
vv=v(4:imax+1,4:imax+1)';
quiver(x,y,uu,vv);
startx = linspace(0,1,50);
starty = ones(size(startx));
streamline(x,y,uu,vv,startx,starty);
streamslice(x,y,uu,vv);
title("Re=3200 streamline");
axis tight;
%% 绘制压力分布
subplot(222)
plot(p(imin:imax,jmax));
title("顶盖压力分布曲线");
set(gcf,'position',[351 199 758 667]);

%% 绘制中轴线速度分布U
subplot(223)
u_inner=u((imax+imin-1)/2,jmin-1:jmax+2);
u_real=zeros(nx+1);
for i=1:nx+1
    u_real(i)=(u_inner(i+1)+u_inner(i))/2;
end
x=linspace(0,Lx,nx+1);
y=u_real(1:nx+1);
plot(x,y,'r-');   %present
hold on;

addpath("E:\Wangxs Terminal\Data base\mats")
axis([0,1,-0.5,1]) ;
load('U3200.mat');
e=U3200(2,:);
f=U3200(1,:);
plot(e,f,'bo');

title(" Vertical line, u");
legend('Re=3200 Present','Re=3200 Ghia et.al.')

%% 绘制中轴线速度分布V
subplot(224)
v_inner=v(imin-1:imax+2,(jmax+jmin-1)/2);
v_real=zeros(nx+1);
for i=1:nx+1
    v_real(i)=(v_inner(i+1)+v_inner(i))/2;
end
x=linspace(0,Lx,nx+1);
y=v_real(1:nx+1);
plot(x,y,'r-');
hold on;

load('V3200.mat');
e=V3200(1,:);
f=V3200(2,:);
plot(e,f,'bo');
title(" Horizontal line, v");
legend('Re=3200 Present','Re=3200 Ghia et.al.')



