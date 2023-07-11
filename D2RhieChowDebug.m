%%
pname=strcat('RC_P',string(RhieChow_or_Not));
save(pname,'p')

% close all;plot(p(imin:imax,jmax)-p(imax,jmax),'r-') ;
% close all
% load('RC_P1.mat')
% plot(p(imin:imax,jmax)-p(imax,jmax),'r-') ; hold on
% % % load('RC_P2.mat')
% % % plot(p(imin:imax,jmax)-p(imax,jmax),'b') ; hold on
% load('RC_P3.mat')
% plot(p(imin:imax,jmax)-p(imax,jmax),'kx') ; hold on
% % load('RC_P4.mat')
% % plot(p(imin:imax,jmax)-p(imax,jmax)) ; hold on
% load('RC_P5.mat')
% plot(p(imin:imax,jmax)-p(imax,jmax),'bo') ; hold on
% load('RC_P6.mat')
% plot(p(imin:imax,jmax)-p(imax,jmax),'r.') ; hold on
% % load('RC_P7.mat')
% % plot(p(imin:imax,jmax)-p(imax,jmax)) ; hold on
% 
% legend('p1','p3','p5','p6')


us=norm(u(imin:imax,jmin:jmax),2);
fprintf("伪势速度为：%s",num2str(us))
