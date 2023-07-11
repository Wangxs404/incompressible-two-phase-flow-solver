function [Fe,PhiE,ChargeQ]= D2Eq_Electronic(phi)
global EK1 EK2 Eepsilon1 Eepsilon2  PhiE_N  PhiE_S
global D2GradX D2GradY Reconstrct
%重构物性参数===================
EK=Reconstrct(EK1,EK2,phi); Eepsilon=Reconstrct(Eepsilon1,Eepsilon2,phi);

%更新电势分布PhiE================
[PhiE] = D2Electronic_PhiE(EK,PhiE_N,PhiE_S);

%更新电荷密度q===================
[PhiE] = D2set_BCPhiE(PhiE,PhiE_N,PhiE_S);%PhiE边条
Difx= Eepsilon.*D2GradX(PhiE);Dify= Eepsilon.*D2GradY(PhiE);
ChargeQ= -(D2GradX(Difx)+D2GradY(Dify));    %电荷密度

%更新电场强度E===================
E.X=-D2GradX(PhiE);     E.Y=-D2GradY(PhiE);

%计算质点所受电场力==============
Emodulus=E.X.*E.X+E.Y.*E.Y;
Fe.X= ChargeQ.*E.X -0.5*Emodulus*D2GradX(Eepsilon);
Fe.Y= ChargeQ.*E.Y -0.5*Emodulus*D2GradY(Eepsilon);

end


%Show Electronic
%     contourf(Coord_x(imin:imax),Coord_y(jmin:jmax),phi(imin:imax,jmin:jmax)',3)
%     hold on;
%     contour(Coord_x(4:1+nx),Coord_y(jmin:jmax),PhiE(imin:imax,jmin:jmax)',30)
%     hold on;
%     contour(Coord_x(imin:imax),Coord_y(jmin:jmax),ChargeQ(imin:imax,jmin:jmax)',20)
%     hold on;