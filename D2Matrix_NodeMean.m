function [FacefinalX,FacefinalY] = D2Matrix_NodeMean(Fv,Fu)  
global imin imax jmin jmax dxi dyi

%mu * DeltaVx
FaceX=0.5*(Fv + circshift(Fv,1,1));
FaceXX=0.5*(FaceX + circshift(FaceX,1,2));
FaceXXX=FaceXX(imin:imax+1,jmin:jmax+1);%Node's V
FaceXXXX=dxi*(circshift(FaceXXX,-1,1)-FaceXXX);
FacefinalX=FaceXXXX(1:end-1,:);
%mu * DeltaUy
FaceY=0.5*(Fu + circshift(Fu,1,2));
FaceYY=0.5*(FaceY + circshift(FaceY,1,1));
FaceYYY=FaceYY(imin:imax+1,jmin:jmax+1);%Node's U
FaceYYYY=dyi*(circshift(FaceYYY,-1,2)-FaceYYY);
FacefinalY=FaceYYYY(:,1:end-1);


end
