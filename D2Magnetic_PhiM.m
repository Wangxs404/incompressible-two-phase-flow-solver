function [PhiM] = D2Magnetic_PhiM(miuM,H0)     %求解电场的电势分布
global nx ny  dxi dyi imin imax jmin jmax dy
global D2Loc D2LocR AX AY
PhiM=zeros(nx+6,ny+6);%初始化
[miuMe,miuMw,miuMn,miuMs] = D2Matrix_FaceMeanC(miuM);  %use the miuM(n+1)   
 %miuM 源于Φ，已施加边条，故不需要再次外插

Lp= -dxi^2 * (miuMe+miuMw) - dyi^2  * (miuMn+miuMs);

Le=  dxi^2 * miuMe;
Lw=  dxi^2 * miuMw;
Ln=  dyi^2 * miuMn;
Ls=  dyi^2 * miuMs;

Lp=reshape(Lp,nx*ny,1);
Le=reshape(Le,nx*ny,1);
Lw=reshape(Lw,nx*ny,1);
Ln=reshape(Ln,nx*ny,1);
Ls=reshape(Ls,nx*ny,1);

% NS固壁，WE自由出流
p_bc_s=  Ls .* D2Loc.s; % p对角线中引入nuemann压力边条
p_bc_n=  Ln .* D2Loc.n; % (k,b)=(1,±dy*H0)
p_bc_w=  Lw .* D2Loc.w;
p_bc_e=  Le .* D2Loc.e;
Lp=Lp+p_bc_s+p_bc_w+p_bc_e+p_bc_n;

% L============== 扣除边界零元，构造5对角稀疏矩阵  
% A的5对角稀疏矩阵
% 扣除无效零边： 系数向量 .* 去边界—位置向量----------------------
Lpp=Lp(:);
Lee=Le.*D2LocR.e; Lww=Lw.*D2LocR.w;
Lnn=Ln.*D2LocR.n; Lss=Ls.*D2LocR.s;
% 由于填充规则，将BWS掐头，而EN掐尾,只取所需段并拼接乘AeLe向量-------------------
Lww= Lww(2:end) ; Lee=Lee(1:end-1);
Lss= Lss(nx+1:end); Lnn=Lnn(1:end-nx);
% 矩阵构造sparse（元素的X坐标向量，元素的Y坐标向量，元素向量，矩阵维度）
Lele=[Lss;Lww;Lpp;Lee;Lnn];
L=sparse(AX,AY,Lele,nx*ny,nx*ny);


%============b

bc_s=  (dy*H0)*Ls.* D2Loc.s ;
bc_n= -(dy*H0)*Ln.* D2Loc.n;% (k,b)=(1,±dy*H0)
bc_w= 0*0*Lw.* D2Loc.w;
bc_e= 0*0*Le.* D2Loc.e;
bPhiM_constu=bc_s+bc_w+bc_e+bc_n;
bM=0-bPhiM_constu;

%============Solve
PhiMV=L\bM;
PhiM(imin : imax , jmin : jmax)=reshape(PhiMV,nx,ny);

end


% function [bPhiE_constu] = D2ElecBC(au,PhiE_N,PhiE_S)
% 
% 
% bc_s=2*PhiE_S *au.s.* D2Loc.s ;
% bc_n=2*PhiE_N *au.n.* D2Loc.n;
% bc_w=0*au.w.* D2Loc.w;
% bc_e=0*au.e.* D2Loc.e;
% bPhiE_constu=bc_s+bc_w+bc_e+bc_n;
% 
% end
