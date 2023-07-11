function [ phin,phiLn,psi,M_change] = D2Eq_PhaseLiuTVD3(uf,vf,phi)
global  M epsilon  dt  How_Many_Phase
%==============================离散求解CH方程===============================
psi=phi.^3-phi-epsilon^2*D2La_Oper(phi);
[psi]=D2set_BCNeu(psi); 
M_change=M*abs((1+phi).*(1-phi));



%% 更新相分数
%===============
phiLn=phi; %截留phiL
switch How_Many_Phase
    case 1
        phin=phi;  %增加该步恢复至单相流
    case 2
        phiA=phi+dt*CHRHS(uf,vf,phi,psi,M_change);
        phiB=0.75*phi + 0.25*phiA + 0.25*dt*CHRHS(uf,vf,phiA,psi,M_change);
        phiC=(1/3)*phi + (2/3)*phiB + (2/3)*dt*CHRHS(uf,vf,phiB,psi,M_change);
        phin=phiC;
end

end

function [RHS] = CHRHS(uf,vf,phi,psi,M_change)
global D2adv   
MLa_Psi=D2DiverGrad(M_change,psi);%M先内乘，再求散度
Fphi=D2adv(uf,vf,uf,vf,phi,phi);%恢复低阶时间格式，但用二阶时间推进
RHS= -Fphi + MLa_Psi;
end

