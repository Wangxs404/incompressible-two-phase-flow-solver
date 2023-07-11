function [N] = D2set_BCPhiM(N,H0)
global imax imin jmax jmin dy
%设置边界条件,由于需对N使用WENO，故外三层都需填充
% PhiE_N,PhiE_S
% N(:,jmin-1)=N(:,jmin);%S
% N(:,jmin-2)=N(:,jmin);
% N(:,jmin-3)=N(:,jmin);
% 
% N(:,jmax+1)=N(:,jmax);%N
% N(:,jmax+2)=N(:,jmax);
% N(:,jmax+3)=N(:,jmax);
% 

N(:,jmin-1)=N(:,jmin)-dy*H0;
N(:,jmin-2)=N(:,jmin)-2*dy*H0;
N(:,jmin-3)=N(:,jmin)-3*dy*H0;

N(:,jmax+1)=N(:,jmax)+dy*H0;
N(:,jmax+2)=N(:,jmax)+2*dy*H0;
N(:,jmax+3)=N(:,jmax)+3*dy*H0;


N(imin-1,:)=N(imin,:);
N(imin-2,:)=N(imin,:);
N(imin-3,:)=N(imin,:);

N(imax+1,:)=N(imax,:);
N(imax+2,:)=N(imax,:);
N(imax+3,:)=N(imax,:);

end

