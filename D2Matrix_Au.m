function [A,au]= D2Matrix_Au(rho,mu)
global imax imin jmax jmin nx ny dx dy dt 
global D2LocR AX AY D2Loc
%动量方程的系数矩阵A  （AU ~=  AV）
mue=circshift(mu,-1,1);muw=circshift(mu,1,1);mun=circshift(mu,-1,2);mus=circshift(mu,1,2);
Ap=1.5*(dx*dy)*rho + 0.5*dt*(dy/dx)*(mue+2*mu+muw) +0.5*dt*(dx/dy)*(mun+2*mu+mus) ;
Ae=-0.5*(dy/dx)*dt*(mue+mu);  Aw=-0.5*(dy/dx)*dt*(muw+mu);  An=-0.5*(dx/dy)*dt*(mun+mu); As=-0.5*(dx/dy)*dt*(mus+mu);

ap=reshape(Ap(imin:imax,jmin:jmax),nx*ny,1);
au.e=reshape(Ae(imin:imax,jmin:jmax),nx*ny,1);
au.w=reshape(Aw(imin:imax,jmin:jmax),nx*ny,1);
au.n=reshape(An(imin:imax,jmin:jmax),nx*ny,1);
au.s=reshape(As(imin:imax,jmin:jmax),nx*ny,1);

% NewBc
% app=ap;




%OldBc
a_bc_s=-au.s .* D2Loc.s;%  A对角线中引入边条
a_bc_n=-au.n .* D2Loc.n; %（系数=边条向量*（+1或-1）），常量引入B_final
a_bc_w=-au.w .* D2Loc.w;
a_bc_e=-au.e .* D2Loc.e;
app=ap+a_bc_s+a_bc_w+a_bc_e+a_bc_n;
        
        
% A的5对角稀疏矩阵
% 扣除无效零边： 系数向量 .* 去边界—位置向量----------------------
aee=au.e.*D2LocR.e; aww=au.w.*D2LocR.w;
ass=au.s.*D2LocR.s; ann=au.n.*D2LocR.n;
% 由于填充规则，将BWS掐头，而ENT掐尾,只取所需段并拼接乘Aele向量-------------------
aww= aww(2:end) ; aee=aee(1:end-1);
ass= ass(nx+1:end); ann=ann(1:end-nx);
% 矩阵构造sparse（元素的X坐标向量，元素的Y坐标向量，元素向量，矩阵维度）
Aele=[ass;aww;app;aee;ann];
A=sparse(AX,AY,Aele,nx*ny,nx*ny);

end



