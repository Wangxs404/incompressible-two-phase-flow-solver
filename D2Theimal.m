function [sigma,Tn]=D2Theimal(phi,rho,uf,vf,T)
global dt D2fXmean D2fYmean sigmafX sigmafY
%Property====
Cp1=1 ; Cp2=1;%比热
K1=0.2 ; K2=0.2;%导热系数
T0=8 ; sigma0=1e-2;%参考温度及其表面张力系数
sigmaT=-1e-3;%sigma-T线性系数

%重构Cp&K==========
Cp=(Cp1+Cp2)/2+(Cp1-Cp2)*phi/2;
K=(K1+K2)/2+(K1-K2)*phi/2;

%求解温度场方程==========
Q1=TRHS(uf,vf,T,K./(rho.*Cp));
Q2=TRHS(uf,vf,T,K./(rho.*Cp));
Tn= T + 1/3*dt*(Q1+2*Q2);

%更新sigma=======================
sigma=sigma0 + sigmaT*(T-T0);
%==
sigmafX=D2fXmean(sigma);sigmafY=D2fYmean(sigma);
end

function [RHS] = TRHS(uf,vf,T,Cof)
global D2GradX D2GradY D2adv
Difx=Cof.*D2GradX(T);Dify=Cof.*D2GradY(T);
CLa_T=D2GradX(Difx)+D2GradY(Dify);    %M先内乘，再求散度

FT=D2adv(uf,vf,uf,vf,T,T);%恢复低阶时间格式，但用二阶时间推进
RHS= -FT + CLa_T;
end

