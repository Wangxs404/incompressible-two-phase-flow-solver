% 绘制气泡定量对比结果

%%
% 绘制气泡定量对比结果=============================================================================================
clc
clear
run D2DataQuantify

addpath("E:\Wangxs Terminal\Data base\mats")
load 'Data10PY'
load 'Data10PY1'
load 'Hysing'

subplot(222)
%%  VC
% Time=Time(1:end);
% Vc=Vc(1:end);
% Yc=Yc(1:end);
figure(1)
plot(Time,Vc,'sk');hold on
plot(Data10PY1(:,1),Data10PY1(:,2),'--b');
plot(Hysing(:,1),Hysing(:,2),'r^');
legend('WENO-FVM','SEM,Cn=0.02','Hysing','location','southeast');
title("Rising Bubble Ratio=10 ");
axis([0,3,0,0.25]) ;

subplot(224)
%%  YC
% load 'Data_Multy'
plot(Time,Yc,'-k');hold on
load('Data10PYyc.mat')
plot(PY10Yc(:,1),PY10Yc(:,2),'--b');
load 'HysingYC.mat'
scatter(HysingYC1(:,1),HysingYC1(:,2),'r^');
legend('WENO-FVM','SEM,Cn=0.02','Hysing','location','northwest');
axis([0,3,0.5,1.2]) ;


subplot(2,2,[1,3])
%%  Shape
load('D2Data_Multy');%FVM
contour(Coord_x(4:1+nx),Coord_y(4:1+ny),phi(4:1+nx,4:1+ny)',1,'k') %等高线图
hold on
addpath 'E:\Wangxs Terminal\BDF V8 test\BFA V8 rb ratio10-5'
load('DataPY.mat');%谱元
plot(DataPY(:,1)+0.5,DataPY(:,2)+1,'b--');
load('HysingShape10.mat');
scatter(HysingShape(:,1),HysingShape(:,2),'r^');
hold on
legend('WENO-FVM','SEM,Cn=0.02','Hysing','location','northwest');
set(gca,'XLim',[0.1 0.9]);set(gca,'YLim',[0.5 1.8]);%限定坐标范围
axis equal;








% title('phi contour ');axis equal;

% %%
% %定量对比
% %T=0-3，Vc
% clc
% clear
% load 'Data_Multy'
% load 'Data10PY'
% load 'Data10PY1'
% load 'Hysing'
% Time=Time(1:end);
% Vc=Vc(1:end);
% Yc=Yc(1:end);
% figure(1)
% plot(Time,Vc,'-k');hold on
% % plot(Data10PY(:,1),Data10PY(:,2),'--r');
% plot(Data10PY1(:,1),Data10PY1(:,2),'--b');
% plot(Hysing(:,1),Hysing(:,2),'^k');
% % legend('WENO-FVM,Cn=0.02','SEM,Cn=0.02','SEM,Cn=0.01','Hysing','location','southeast');
% legend('WENO-FVM,Cn=0.01','SEM,Cn=0.01','Hysing','location','southeast');
% title("Rinsing bubble Ratio=10 ");
% axis([0,3,0,0.25]) ;
%
% %%
% %定量对比T=3s 时的，气泡现状（ratio=10）
% addpath 'E:\Wangxs Terminal\BDF V8 test\BFA V8 rb ratio10-5'
% load('DataPY.mat');%谱元
% plot(DataPY(:,1)+0.5,DataPY(:,2)+1,'b--');
% hold on
%
% load('HysingShape10.mat');%谱元
% scatter(HysingShape(:,1),HysingShape(:,2),'r^');
% hold on
%
% load('Data_Multy.mat');%FVM
% [xx,yy]=find(phi<0.3 & phi>-0.3);%返回相界面下标
% xx=(xx-3)/nx*Lx-0.5;yy=(yy-3)/ny*Ly-1;% 将下标还原到界面坐标（并移动与谱元对比）
% plot(xx+0.5,yy+1,'.k');
% set(gca,'XLim',[0.1 0.9]);set(gca,'YLim',[0.5 1.8]);%限定坐标范围
% legend('SEM,Cn=0.02','Hysing','WENO-FVM,Cn=0.01','location','northwest');
% %%
% %Yc
% load 'Data_Multy'
% plot(Time,Yc,'-k');hold on
% load('Data10PYyc.mat')
% plot(PY10Yc(:,1),PY10Yc(:,2),'--b');
% load 'HysingYC.mat'
% scatter(HysingYC1(:,1),HysingYC1(:,2),'^k');
% legend('WENO-FVM,Cn=0.01','SEM,Cn=0.02','Hysing','location','northwest');
% axis([0,3,0.5,1.2]) ;
%
%
% %% 质量泄漏
% load('Data_Initial.mat')
% phiV=(1+phi)/2  ;
% Mass_Initial =sum(sum(rho(imin:imax,jmin:jmax).*phiV(imin:imax,jmin:jmax)));
%
% load('Data10000.mat')
% phiV=(1+phi)/2  ;
% Mass_10000 =sum(sum(rho(imin:imax,jmin:jmax).*phiV(imin:imax,jmin:jmax)));
%
% load('Data20000.mat')
% phiV=(1+phi)/2  ;
% Mass_20000 =sum(sum(rho(imin:imax,jmin:jmax).*phiV(imin:imax,jmin:jmax)));
%
% load('Data30000.mat')
% phiV=(1+phi)/2  ;
% Mass_30000 =sum(sum(rho(imin:imax,jmin:jmax).*phiV(imin:imax,jmin:jmax)));
%
% load('Data40000.mat')
% phiV=(1+phi)/2  ;
% Mass_40000 =sum(sum(rho(imin:imax,jmin:jmax).*phiV(imin:imax,jmin:jmax)));
%
% load('Data50000.mat')
% phiV=(1+phi)/2  ;
% Mass_50000 =sum(sum(rho(imin:imax,jmin:jmax).*phiV(imin:imax,jmin:jmax)));
%
% load('Data60000.mat')
% phiV=(1+phi)/2  ;
% Mass_60000 =sum(sum(rho(imin:imax,jmin:jmax).*phiV(imin:imax,jmin:jmax)));
%
% a0=Mass_Initial/Mass_Initial;
% a1=Mass_10000/Mass_Initial;
% a2=Mass_20000/Mass_Initial;
% a3=Mass_30000/Mass_Initial;
% a4=Mass_40000/Mass_Initial;
% a5=Mass_50000/Mass_Initial;
% a6=Mass_60000/Mass_Initial;
% A=[a0 a1 a2 a3 a4 a5 a6];B=[1 1 1 1 1 1 1];
% plot(B,'--r');hold on;
% plot(A,'-b');axis([1,6,0.5,1.5]) ;
% legend('rho initial','rho now','location','northwest');




