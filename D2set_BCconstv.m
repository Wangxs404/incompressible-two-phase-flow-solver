function [bu_constv] = D2set_BCconstv(av,v)
global D2Loc imin imax jmin jmax nx ny
global v_bottom  v_top  v_left v_right
%New
%b-below
% vs=[v(imin:imax,jmin-1);zeros(nx*(ny-1),1) ];
% vn=[zeros(nx*(ny-1),1) ;v(imin:imax,jmax+1) ];
% vw=repmat(v(imin-1,jmin:jmax),nx,1);vw=reshape(vw,nx*ny,1);
% ve=repmat(v(imax+1,jmin:jmax),nx,1);ve=reshape(ve,nx*ny,1);
% 
% bc_s= vs  .*av.s.* D2Loc.s ;
% bc_n= vn  .*av.n.* D2Loc.n;
% bc_w= vw  .*av.w.* D2Loc.w;
% bc_e= ve  .*av.e.* D2Loc.e;
% bu_constv=bc_s+bc_w+bc_e+bc_n;

%Old
        bc_s=2*v_bottom *av.s.* D2Loc.s ;
        bc_n=2*v_top    *av.n.* D2Loc.n;
        bc_w=0  *av.w.* D2Loc.w;
        bc_e=0  *av.e.* D2Loc.e;
        bu_constv=bc_s+bc_w+bc_e+bc_n;

end
