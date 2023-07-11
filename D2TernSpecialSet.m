%  load('Data46000');    % 续算时开启
%  Res_Freq=50  ;  dat_Freq=1000;

switch initial_type   
    case 3 %Za
        [u,v]  =D2InitialUVza(Coord_x,Coord_y,imax,imin,jmax,jmin);  %Zalesak's Disk时开启,该算例通用性不强，参考改写E:\Wangxs\Multiphase Flow\phase move CH\phase move CH2 zalesak
        [uL,vL]=D2InitialUVza(Coord_x,Coord_y,imax,imin,jmax,jmin);  %Zalesak's Disk时开启
    case 12
        [u,v]  =D2InitialUVdb(Coord_x,Coord_y,imax,imin,jmax,jmin);  %Drop breakup时开启
        [uL,vL]=D2InitialUVdb(Coord_x,Coord_y,imax,imin,jmax,jmin);  %Drop breakup时开启
    case 6 %自动开启 singel drop 定量对比的向量初始化
        
        [u,v]  =D2InitialUVdb(Coord_x,Coord_y,imax,imin,jmax,jmin);  %Drop breakup时开启
        [uL,vL]=D2InitialUVdb(Coord_x,Coord_y,imax,imin,jmax,jmin);  %Drop breakup时开启       
end


if Thermal ==  "ON";  T=D2Initial_T(T);end
 


% 初始化表面活性剂浓度场==================================

csi=0.8;   %强度
cs = csi * (1-abs(phi.A))+0.01;     % Concentrations of Surfactant Bulk浓度不为0，需考虑附加小值
cpsi=0*phi.A;%便于未考虑表面活性剂计算时，export2可以运行





%Initial Mass
MassLeaf=0*speye(istep_max,1);
% rho=rho_Light*phi+rho_Heavy*(1-phi);   %rho_Light 对应 phi=1
% rho=(rho_Light+rho_Heavy)/2+(rho_Light-rho_Heavy)*phi/2;
% Mass_Initial = sum(sum(rho(imin:imax,jmin:jmax).*phi(imin:imax,jmin:jmax)));

% phiVol=phi ; phiVol(phiVol<0)=0; 
% VolumnIni=sum(sum(phiVol))*(dx*dy);
