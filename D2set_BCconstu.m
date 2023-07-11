function [bu_constu] = D2set_BCconstu(au,u)
global D2Loc imin imax jmin jmax nx ny
global u_bottom u_top  u_left  u_right
%New
% b-below
% us=[u(imin:imax,jmin-1);zeros(nx*(ny-1),1) ];
% un=[zeros(nx*(ny-1),1) ;u(imin:imax,jmax+1) ];
% uw=repmat(u(imin-1,jmin:jmax),nx,1);uw=reshape(uw,nx*ny,1);
% ue=repmat(u(imax+1,jmin:jmax),nx,1);ue=reshape(ue,nx*ny,1);

% bc_s= us  .*au.s.* D2Loc.s ;
% bc_n= un  .*au.n.* D2Loc.n;
% bc_w= uw  .*au.w.* D2Loc.w;
% bc_e= ue  .*au.e.* D2Loc.e;
% bu_constu=bc_s+bc_w+bc_e+bc_n;

%Old
bc_s=2*u_bottom *au.s.* D2Loc.s ;
bc_n=2*u_top    *au.n.* D2Loc.n;
bc_w=2*u_left   *au.w.* D2Loc.w;
bc_e=2*u_right  *au.e.* D2Loc.e;
bu_constu=bc_s+bc_w+bc_e+bc_n;

end
