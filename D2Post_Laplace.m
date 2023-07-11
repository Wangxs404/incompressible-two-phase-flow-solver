clc
clear
global imin imax jmin jmax nx ny sigma0 r
load('D2Data_Multy');%FVM

%绘制流线图
% figure(1)
subplot(121)
x=linspace (0 ,Lx, nx+1);
y=linspace (0 ,Ly, ny+1);
uu=u(4:imax+1,4:imax+1)';
vv=v(4:imax+1,4:imax+1)';
quiver(x,y,uu,vv);
startx = linspace(0,1,50);
starty = ones(size(startx));
streamline(x,y,uu,vv,startx,starty);
streamslice(x,y,uu,vv);
title("streamline");
axis equal;
set(gca,'xtick',[],'ytick',[]) ; 

%绘制压力分布
% figure(2)
subplot(122)
plot(p(imin+nx/2,jmin:jmax));

V=sqrt(u.^2+v.^2);
normV_2=norm(V(imin:imax,jmin:jmax),2);
normV_INF=norm(V(imin:imax,jmin:jmax),inf);
CaNUM=normV_INF*mu_Light/sigma0;
fprintf("V-2范数=%s \nV-∞范数=%s \nCa= %s \n",num2str(normV_2),num2str(normV_INF),num2str(CaNUM));

% normU=norm(u(imin:imax,jmin:jmax),inf);
% normV=norm(v(imin:imax,jmin:jmax),inf);
sigmaCalcu=(max(p(imin+nx/2,jmin:jmax))-min(p(imin+nx/2,jmin:jmax)))*r;
error=(sigma0-sigmaCalcu)/sigma0;
fprintf("Laplace Law误差=%s",strcat( num2str(error*100),"%"));

%% 绘制三个La数的压力分布曲线，并计算Laplace误差

load("La120.mat")
plot(linspace(0.5/128,1-0.5/128,64),p(imin+nx/2,jmin:2:jmax)-min(p(imin+nx/2,jmin:jmax)),'r--','MarkerSize',5);
hold on

load("La1200.mat")
plot(linspace(0.5/128,1-0.5/128,64),p(imin+nx/2,jmin+1:2:jmax+1)-min(p(imin+nx/2,jmin:jmax)),'g-.','MarkerSize',5);

load("La12000.mat")
plot(linspace(0.5/128,1-0.5/128,128),p(imin+nx/2,jmin:jmax)-min(p(imin+nx/2,jmin:jmax)),'k-');

plot(linspace(0,0.5,20),5*ones(20,1),'b--');

legend('La=120','La=1200','La=12000','Theory resolution','location','northeast');
axis tight;
xlabel('Horizontal coordinates x');ylabel('ΔP');
axis([0,0.5,0,7])
% linspace(0.5/128)
