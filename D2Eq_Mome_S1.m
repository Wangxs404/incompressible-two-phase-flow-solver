function [u_star,v_star]=D2Eq_Mome_S1(u,v,uL,vL,uf,vf,ufL,vfL,p,mu,rho,rhoL,rhoLL, Mccx, Mccy,psi,phi,sigma,Fe,Fm)
global imax imin jmax jmin gx gy  nx ny  Swith_SM upsilon0
%%
%==========================初始化和准备===================================
u_star=zeros(imax+3,jmax+3);v_star=zeros(imax+3,jmax+3);
uR=2*u-uL;    vR=2*v-vL;

if Swith_SM==1 
    mX=Mccx;   mY=Mccy;
else
    mX=rhoL.*u;   mY=rhoL.*v;
end

[mX] = D2set_BCNeu(mX); [mY] = D2set_BCNeu(mY);  %内外动量暂用纽曼边条，至少不等同速度边条，也不应直接计算后不设置边条，可能出现错误的动量穿越
%=============计算系数矩阵=========
% %系数矩阵A
% [Au,au]= D2Matrix_Au(rho,mu);
% [Av,av]= D2Matrix_Av(rho,mu);
%常系数矩阵A
[Au,au]= D2Matrix_Au(ones(nx+6,ny+6),upsilon0*ones(nx+6,ny+6));
[Av,av]= D2Matrix_Av(ones(nx+6,ny+6),upsilon0*ones(nx+6,ny+6));
%-----------------------
[Dtx,Dty]=D2Diver_Dtrans(uR,vR,mu);%黏性项2
[bu]=D2Calcu_b(u,uL,mX,gx,@D2Gradx_Matrix,Dtx,rho,rhoL,rhoLL,uf,vf,ufL,vfL,psi,phi,p,@D2Grady_Matrix,sigma,Fe.X,Fm.X);
[bv]=D2Calcu_b(v,vL,mY,gy,@D2Grady_Matrix,Dty,rho,rhoL,rhoLL,uf,vf,ufL,vfL,psi,phi,p,@D2Gradx_Matrix,sigma,Fe.Y,Fm.Y);
%================================计算bu====================================
bu_final=bu -D2set_BCconstu(au);  %Ghost cell的常量部分，引入bfinal
Lu = ichol(Au) ;      %一些情况下，不进行预处理反而算的更快
[uv,~]= pcg (Au,bu_final,1e-6,20,Lu,Lu');%同等情况下pcg优于cgs
u_star(imin : imax , jmin : jmax)=reshape(uv,nx,ny);  %赋值给U_start
%================================计算bv====================================
bv_final=bv - D2set_BCconstv(av);
Lv = ichol(Av);
[vv,~]= pcg (Av,bv_final,1e-6,20,Lv,Lv');
v_star(imin : imax , jmin : jmax)=reshape(vv,nx,ny);  %赋值给V_start

end
