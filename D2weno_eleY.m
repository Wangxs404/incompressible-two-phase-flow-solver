function [s2,s1,py,n1,n2,n3]= D2weno_eleY(F)
global imin imax jmin jmax

py=F(imin:imax , jmin-1:jmax);
s1=F(imin:imax , jmin-2:jmax-1);
s2=F(imin:imax , jmin-3:jmax-2);
n1=F(imin:imax , jmin:jmax+1);
n2=F(imin:imax , jmin+1:jmax+2);
n3=F(imin:imax , jmin+2:jmax+3);

end