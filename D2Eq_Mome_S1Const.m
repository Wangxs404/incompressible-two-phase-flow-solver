function [u_star,v_star]=D2Eq_Mome_S1Const(u,v,uL,vL,uf,vf,ufL,vfL,p,pL,mu,rho,rhoL,rhoLL, Mccx, Mccy,psi,phi,sigma,Fe,Fm)
global imax imin jmax jmin gx gy  nx ny 
global Au Av au av  LAu LAuT LAv LAvT Swith_Mcc
%%
%==========================初始化和准备===================================
u_star=zeros(imax+3,jmax+3);v_star=zeros(imax+3,jmax+3);
if Swith_Mcc=="ON"    % 守恒通量
    mX=Mccx;   mY=Mccy;
else
    mX=rhoL.*u;   mY=rhoL.*v;
%         mX=2*rhoL.*u-rhoLL.*uL;   mY=2*rhoL.*v-rhoLL.*vL;
end
[mX] = D2set_BCNeu(mX); [mY] = D2set_BCNeu(mY);
%转置黏性项
[Dtx,Dty]=D2Diver_Dtrans(2*u-uL,2*v-vL,mu);
Dtx=Dtx./rho(imin:imax,jmin:jmax);
Dty=Dty./rho(imin:imax,jmin:jmax);

%常系数b
[bu]=D2Calcu_bConst(u,uL,mX,gx,@D2Gradx_Matrix,rho,rhoL,rhoLL,mu,uf,vf,ufL,vfL,psi,phi,p,pL,@D2Grady_Matrix,sigma,Dtx,Fe.X,Fm.X);
[bv]=D2Calcu_bConst(v,vL,mY,gy,@D2Grady_Matrix,rho,rhoL,rhoLL,mu,uf,vf,ufL,vfL,psi,phi,p,pL,@D2Gradx_Matrix,sigma,Dty,Fe.Y,Fm.Y);
%================================常系数算法计算u_star====================================
bu_final=bu -D2set_BCconstu(au,u);  %Ghost cell的常量部分，引入bfinal
[uv,~]= pcg(Au,bu_final,1e-6,20,LAu,LAuT);%pcg 效率=5倍的chol
u_star(imin : imax , jmin : jmax)=reshape(uv,nx,ny);  %赋值给U_start
%================================常系数算法计算v_star====================================
bv_final=bv - D2set_BCconstv(av,v);
[vv,~]= pcg (Av,bv_final,1e-6,20,LAv,LAvT);
v_star(imin : imax , jmin : jmax)=reshape(vv,nx,ny);

end
