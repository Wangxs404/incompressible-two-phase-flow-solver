clc
clear

run D2DataQuantify.m

addpath("D:\Wangxs\Data base\mats")
load('RTdata.mat')

subplot(2,3,[1,4])
contourf(Coord_x(4:1+nx),Coord_y(4:1+ny),phi(4:1+nx,4:1+ny)',5) 
title('phi contour ');axis equal;

subplot(2,3,[2,3,5,6])
plot(Time,RT_wall-2,'k-');hold on;
plot(RTdata(:,1),RTdata(:,2),'ro');hold on;
plot(Time,RT_midd-2,'k-');
axis([0,2.7,-1.5,1.5]) ;
legend('Present','Huang et.al','location','northwest');
title(" Location of the interfaceof the Rayleigh-Taylor instability","Re=256,At=0.5");
xlabel("t");ylabel("Location");


