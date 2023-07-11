function [Dtx,Dty]=D2Diver_Dtrans(uR,vR,mu)  %黏性项2的离散
global dxi dyi
[muX,muY] = D2Matrix_FaceMeanF(mu);

muGux =muX .* D2Gradx_face(uR);
Dtxx=dxi*(circshift(muGux,-1,1)-muGux);  Dtx1=Dtxx(1:end-1 ,:);%Dx1

muGvy =muY .* D2Grady_face(vR);
Dtyy=dyi*(circshift(muGvy,-1,2)-muGvy);  Dty1=Dtyy(:,1:end-1 );%Dy1

[vx,uy] = D2Matrix_NodeMean(vR,uR) ;
D2v=muY .* vx;
D2u=muX .* uy;

[Dtx2] = D2Grady_face2(D2v) ;                                    %Dx2
[Dty2] = D2Gradx_face2(D2u) ;                                    %Dx2

Dtx= Dtx1+Dtx2;
Dty= Dty1+Dty2;

end

