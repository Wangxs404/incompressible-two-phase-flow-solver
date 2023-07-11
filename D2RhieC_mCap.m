function [mcapXP,mcapXE,mcapXf,mcapYP,mcapYN,mcapYf] = D2RhieC_mCap(u,v,uL,vL,uf,vf,ufL,vfL,rhoL,rhoLL)
global imin imax jmin jmax
mcapXP=   2 *rhoL(imin-1:imax,jmin:jmax).*u(imin-1:imax,jmin:jmax) ...
        - 0.5*rhoLL(imin-1:imax,jmin:jmax).*uL(imin-1:imax,jmin:jmax)  ;
mcapXE=   2 *rhoL(imin:imax+1,jmin:jmax).*u(imin:imax+1,jmin:jmax) ...
        - 0.5*rhoLL(imin:imax+1,jmin:jmax).*uL(imin:imax+1,jmin:jmax)  ;
%=================== 
RhoLx_f=2 ./ ( 1./rhoL(imin-1:imax , jmin:jmax) + 1./rhoL(imin:imax+1 , jmin:jmax) );
RhoLLx_f=2 ./ ( 1./rhoLL(imin-1:imax , jmin:jmax) + 1./rhoLL(imin:imax+1 , jmin:jmax) );

mcapXf=  2*RhoLx_f.*uf - 0.5*RhoLLx_f.*ufL ;
%%
mcapYP=   2 *rhoL(imin:imax,jmin-1:jmax).*v(imin:imax,jmin-1:jmax) ...
        - 0.5*rhoLL(imin:imax,jmin-1:jmax).*vL(imin:imax,jmin-1:jmax)  ;
mcapYN=   2 *rhoL(imin:imax,jmin:jmax+1).*v(imin:imax,jmin:jmax+1) ...
        - 0.5*rhoLL(imin:imax,jmin:jmax+1).*vL(imin:imax,jmin:jmax+1)  ;
    
RhoLy_f=2 ./ ( 1./rhoL(imin:imax , jmin-1:jmax) + 1./rhoL(imin:imax , jmin:jmax+1) );
RhoLLy_f=2 ./ ( 1./rhoLL(imin:imax , jmin-1:jmax) + 1./rhoLL(imin:imax , jmin:jmax+1) );

mcapYf=  2*RhoLy_f.*vf - 0.5*RhoLLy_f.*vfL ;
% k=1


end