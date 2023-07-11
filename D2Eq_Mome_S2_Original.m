function [ppie]=D2Eq_Mome_S2_Original(u_star,v_star,rho)
global  imax imin jmax jmin dt nx ny dxi dyi 
% 中心差分计算通量
%%
%==================系数矩阵：▽·（1/ρ*▽P'）===========
ppie=zeros(imax+3,jmax+3);
[L] = D2Matrix_Laplace(rho);

[FaceUX,~] = D2Matrix_FaceMeanF(u_star) ;
[~,FaceVY] = D2Matrix_FaceMeanF(v_star) ;
adv_Ustar= dxi*( FaceUX(2:end,:)-FaceUX(1:end-1,:))...
    +dyi*(FaceVY(:,2:end)-FaceVY(:,1:end-1))  ;

Bp1=1.5/dt * adv_Ustar; b_poison=reshape(Bp1,nx*ny,1);
%%
%===================Directly Solve==============================
pv=sparse(L)\sparse(b_poison);   % Poison 方程非正定，不能使用共轭梯度法
ppie(imin : imax , jmin : jmax)=reshape(pv,nx,ny); %将解得的P填入网格

end
