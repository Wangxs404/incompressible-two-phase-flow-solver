function [T] = D2set_BCT(T)
global imax imin jmax jmin
% NS
        % T ：无滑移 固壁  绝热
        % T ：无滑移 固壁  绝热
        T(:,jmin-1)=T(:,jmin);
        T(:,jmin-2)=T(:,jmin);
        T(:,jmin-3)=T(:,jmin);
        
        T(:,jmax+1)=T(:,jmax);
        T(:,jmax+2)=T(:,jmax);
        T(:,jmax+3)=T(:,jmax);
% % EW
%         % T ：周期 
%         % T ：周期 
%         T(imin-1,:)=T(imax,:);
%         T(imin-2,:)=T(imax-1,:);
%         T(imin-3,:)=T(imax-2,:);
%         
%         T(imax+1,:)=T(imin,:);
%         T(imax+2,:)=T(imin+1,:);
%         T(imax+3,:)=T(imin+2,:);

% EW
        % T ：开放
        % T ：开放
        T(imin-1,:)=T(imin,:);
        T(imin-2,:)=T(imin,:);
        T(imin-3,:)=T(imin,:);
        
        T(imax+1,:)=T(imax,:);
        T(imax+2,:)=T(imax,:);
        T(imax+3,:)=T(imax,:);
end

