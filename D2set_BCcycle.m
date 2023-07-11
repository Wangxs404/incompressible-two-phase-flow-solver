function [F] = D2set_BCcycle(F)
global imin imax jmin jmax
%设置边界条件,由于需对F使用WENO，故外三层都需填充

F(:,jmax+1)=F(:,jmin+1);
F(:,jmax+2)=F(:,jmin+2);
F(:,jmax+3)=F(:,jmin+3);
F(:,jmin-1)=F(:,jmax-1);
F(:,jmin-2)=F(:,jmax-2);
F(:,jmin-3)=F(:,jmax-3);


F(imax+1,:)=F(imin+1,:);
F(imax+2,:)=F(imin+2,:);
F(imax+3,:)=F(imin+3,:);
F(imin-1,:)=F(imax-1,:);
F(imin-2,:)=F(imax-2,:);
F(imin-3,:)=F(imax-3,:);
end

