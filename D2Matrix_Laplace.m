function [L,Le,Lw,Ln,Ls] = D2Matrix_Laplace(rho)     %此为Poisson方程全隐式求解，需Mome Step2 配合引入边界条件
global nx ny  dxi dyi refere_Pressure 
global D2Loc D2LocR AX AY
  
[rhoe,rhow,rhon,rhos] = D2Matrix_FaceMeanC(rho);  %use the rho(n+1)    %不必空间积分，多余了
rhoeT=1./rhoe;   rhowT=1./rhow;   rhonT=1./rhon;   rhosT=1./rhos;   %rho 源于Φ，已施加边条，故不需要再次外插

Lp= -dxi^2 * (rhoeT+rhowT) - dyi^2  * (rhonT+rhosT);

Le=  dxi^2 * rhoeT;
Lw=  dxi^2 * rhowT;
Ln=  dyi^2 * rhonT;
Ls=  dyi^2 * rhosT;

Lp=reshape(Lp,nx*ny,1);
Le=reshape(Le,nx*ny,1);
Lw=reshape(Lw,nx*ny,1);
Ln=reshape(Ln,nx*ny,1);
Ls=reshape(Ls,nx*ny,1);

p_bc_s=Ls .* D2Loc.s ;%p对角线中引入nuemann压力边条
p_bc_n=Ln .* D2Loc.n; % 仅适用于δP=0，等于给Lp对应点系数+1
p_bc_w=Lw .* D2Loc.w;
p_bc_e=Le .* D2Loc.e;
Lp=Lp+p_bc_s+p_bc_w+p_bc_e+p_bc_n;

%设置压力零点
switch refere_Pressure
    case 1
        Lp(1)=1e5;            %左下
    case 2
        Lp(nx)=1e5;           %右下
    case 3
        Lp(nx*ny-nx+1)=1e5;   %左上
    case 4
        Lp(nx*ny)=1e5;        %右上
    case 5
        Lp(nx*ny/2-nx/2)=1e5; %体心
    case 6
        Lp(nx*ny-nx/2)=1e5;   %顶中点
    case 7
        Lp(nx/2)=1e5;         %底中点
end

% L============== 扣除边界零元，构造7对角稀疏矩阵  
% A的5对角稀疏矩阵
% 扣除无效零边： 系数向量 .* 去边界—位置向量----------------------
Lpp=Lp(:);
Lee=Le.*D2LocR.e; Lww=Lw.*D2LocR.w;
Lnn=Ln.*D2LocR.n; Lss=Ls.*D2LocR.s;
% 由于填充规则，将BWS掐头，而ENT掐尾,只取所需段并拼接乘AeLe向量-------------------
Lww= Lww(2:end) ; Lee=Lee(1:end-1);
Lss= Lss(nx+1:end); Lnn=Lnn(1:end-nx);
% 矩阵构造sparse（元素的X坐标向量，元素的Y坐标向量，元素向量，矩阵维度）
Lele=[Lss;Lww;Lpp;Lee;Lnn];
L=sparse(AX,AY,Lele,nx*ny,nx*ny);

end
