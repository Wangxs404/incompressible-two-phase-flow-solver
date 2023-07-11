function [phi] = D2set_BCWet(phi,phiL,sigma)
global imax imin jmax jmin   theta dy dx epsilon
%设置边界条件,由于需对phi使用WEphiO，故外三层都需填充
lamda=(3*sigma(imin:imax,jmin)*epsilon) / (2*sqrt(2));
phiS=2*phi(imin:imax,jmin)-phiL(imin:imax,jmin);
a= 0.5*sigma(imin:imax,jmin)./lamda*cos(theta).*(pi/2*cos(pi/2*phiS));
%S
phi(imin:imax,jmin-1)=phi(imin:imax, jmin ) + a*dy;
phi(imin:imax,jmin-2)=phi(imin:imax,jmin-1) + a*dy;
phi(imin:imax,jmin-3)=phi(imin:imax,jmin-2) + a*dy;

%%
%phi
phi(:,jmax+1)=phi(:,jmax);
phi(:,jmax+2)=phi(:,jmax);
phi(:,jmax+3)=phi(:,jmax);

%W
phi(imin-1,:)=phi(imin,:);
phi(imin-2,:)=phi(imin,:);
phi(imin-3,:)=phi(imin,:);

%E
phi(imax+1,:)=phi(imax,:);
phi(imax+2,:)=phi(imax,:);
phi(imax+3,:)=phi(imax,:);

end

