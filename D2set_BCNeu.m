function [N] = D2set_BCNeu(N)
global imax imin jmax jmin
%设置边界条件,由于需对N使用WENO，故外三层都需填充
N(:,jmin-1)=N(:,jmin);
N(:,jmin-2)=N(:,jmin);
N(:,jmin-3)=N(:,jmin);

N(:,jmax+1)=N(:,jmax);
N(:,jmax+2)=N(:,jmax);
N(:,jmax+3)=N(:,jmax);

N(imin-1,:)=N(imin,:);
N(imin-2,:)=N(imin,:);
N(imin-3,:)=N(imin,:);

N(imax+1,:)=N(imax,:);
N(imax+2,:)=N(imax,:);
N(imax+3,:)=N(imax,:);

end

