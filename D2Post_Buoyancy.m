

load('D2Data_Multy');%FVM

%等高线图
contour(Coord_x(4:1+nx),Coord_y(4:1+ny),phi(4:1+nx,4:1+ny)',1,'k') 
hold on

% %流线图
% x=linspace (0 ,Lx, nx+1);
% y=linspace (0 ,Ly, ny+1);
% uu=u(4:imax+1,4:imax+1)';
% vv=v(4:imax+1,4:imax+1)';
% quiver(x,y,uu,vv);
% startx = linspace(0,1,50);
% starty = ones(size(startx));
% streamline(x,y,uu,vv,startx,starty);
% streamslice(x,y,uu,vv);
% title("streamline");
% axis equal;
% set(gca,'xtick',[],'ytick',[]) ; 


% [p]=D2set_BCNeu(p);  
% LaP=D2La_Oper(p);
% error=norm(LaP(imin:imax,jmin:jmax),2);
% fprintf("Laplace P=%s",num2str(error));

V=sqrt(u.^2+v.^2);
normV_2=norm(V(imin:imax,jmin:jmax),2);
normV_INF=norm(V(imin:imax,jmin:jmax),inf);
CaNUM=normV_INF*mu_Light/sigma0;
fprintf("Ca= %s \n",num2str(CaNUM));

contour(Coord_x(4:1+nx),Coord_y(4:1+ny),p(4:1+nx,4:1+ny)',8,'r') 

