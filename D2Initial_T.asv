function T=D2Initial_T(T)
global ny nx imin imax jmin jmax
T_High=16;    T_Low=0;%对应 W-E 壁面
Tvec=linspace(T_Low,T_High,nx)';Tin=repmat(Tvec,1,ny);
T(imin:imax,jmin:jmax)=Tin;
end
