function [imin,imax,jmin,jmax,Coord_x,Coord_y,dx,dy,dxi,dyi] = D2Mesh(Lx,Ly)
global nx ny ds
imin=4; imax=imin+nx-1;    %imin,imax分别对应内网格始末节点
jmin=4; jmax=jmin+ny-1;
x(imin : imax+1)=linspace (0 ,Lx, nx+1);
y(jmin : jmax+1)=linspace (0 ,Ly, ny+1);
Coord_x(imin : imax)=0.5*(x(imin : imax)+x(imin+1:imax+1));  %coordinate of grid nodeX
Coord_y(jmin : jmax)=0.5*(y(jmin : jmax)+y(jmin+1:jmax+1));
dx=x(imin+1)-x(imin );    % Create mesh sizes
dy=y(jmin+1)-y(jmin );
dxi=1/dx ;
dyi=1/dy ;
ds=dx*dy;
end



% dx=2/nx;
% Coord_x2=dx/2+(I-1)*dx;