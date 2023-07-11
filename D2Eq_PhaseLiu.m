function [ phin,phiLn,psi,M_change] = D2Eq_PhaseLiu(uf,vf,phi)
global  M epsilon  dt  How_Many_Phase
%==============================离散求解CH方程===============================
psi=phi.^3-phi-epsilon^2*D2La_Oper(phi);
[psi]=D2set_BCNeu(psi); 
M_change=M*abs((1+phi).*(1-phi));

K1=CHRHS(uf,vf,phi,psi,M_change);
K2=CHRHS(uf,vf,phi+0.75*dt*K1,psi,M_change);

%% 更新相分数
%===============
switch How_Many_Phase
    case 1
        phin=phi; phiLn=phi;  %增加该步恢复至单相流
    case 2
        phiLn=phi; %截留phiL
        phin=phi + 1/3*dt*(K1+2*K2);
end

end

function [RHS] = CHRHS(uf,vf,phi,psi,M_change)
global D2adv D2GradX D2GradY
Difx=M_change.*D2GradX(psi);Dify=M_change.*D2GradY(psi);
MLa_Psi=D2GradX(Difx)+D2GradY(Dify);    %M先内乘，再求散度

Fphi=D2adv(uf,vf,uf,vf,phi,phi);%恢复低阶时间格式，但用二阶时间推进
RHS= -Fphi + MLa_Psi;
end

