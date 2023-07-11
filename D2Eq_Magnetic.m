function [Fm,PhiM,H,Emodulus]= D2Eq_Magnetic(phi)
global miuM1 miuM2 Reconstrct D2GradX D2GradY H0
%重构物性参数===================
miuM=Reconstrct(miuM1,miuM2,phi);
%更新磁势分布PhiM================
[PhiM] = D2Magnetic_PhiM(miuM,H0);
[PhiM] = D2set_BCPhiM(PhiM,H0);
%更新磁通量======================
H.X=-D2GradX(PhiM);     H.Y=-D2GradY(PhiM);
%计算磁应力=====================
%Emodulus=H.X.*H.X+H.Y.*H.Y; kai0=2.2; miu0=4*pi*1e-7;
Emodulus=H.X.*H.X+H.Y.*H.Y; kai0=2.2; miu0=4*pi*1e-7;
Fm.X= 0.5*miu0*kai0*D2GradX(Emodulus);
Fm.Y= 0.5*miu0*kai0*D2GradY(Emodulus);
end

%Show Electronic
%     contour(Coord_x(imin:imax),Coord_y(jmin:jmax),Emodulus(imin:imax,jmin:jmax)',1)
%     hold on;
%     contour(Coord_x(jmin:jmax),Coord_y(jmin:jmax),PhiM(imin:imax,jmin:jmax)',30)
%     hold on;
%    quiver(Coord_x(jmin:jmax),Coord_y(jmin:jmax),H.X(imin:imax,jmin:jmax)',H.Y(imin:imax,jmin:jmax)',1)
% %     contour(Coord_x(imin:imax),Coord_y(jmin:jmax),ChargeQ(imin:imax,jmin:jmax)',20)
%     hold on;