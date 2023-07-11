function [w2,w1,px,e1,e2,e3]= D2weno_eleX(F)
global imin imax jmin jmax

px=F(imin-1:imax , jmin:jmax);
w1=F(imin-2:imax-1 , jmin:jmax);
w2=F(imin-3:imax-2 , jmin:jmax);
e1=F(imin:imax+1 , jmin:jmax);
e2=F(imin+1:imax+2 , jmin:jmax);
e3=F(imin+2:imax+3 , jmin:jmax);

end