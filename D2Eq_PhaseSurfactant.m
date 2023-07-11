function [ csn,cpsi] = D2Eq_PhaseSurfactant(uf,vf,cs,phi)
global   dt epsilon Pi Ex cM sigma0
%==============================离散求解CH方程===============================
%参数==========================
Beta = 3*sigma0/(4*epsilon) ; 
%Surfactant化学势===============
cpsi= 4*Beta*Pi*log(cs./(1-cs))... %Zong
     -Beta*((phi.^2-1).^2)...
     +(2*Beta/Ex)*phi.^2;
% cpsi= Pi*log(cs./(1-cs))...  %使用Soligo推出的cpsi,也可得到合理分布的浓度场
%      -0.5*((phi.^2-1).^2)...
%      +(0.5/Ex)*phi.^2;
 
[cpsi]=D2set_BCNeu(cpsi); 
%===============================

K1=CHRHS(uf,vf,cs,cpsi,cM);
K2=CHRHS(uf,vf,cs+0.75*dt*K1,cpsi,cM);

%% 更新相分数=============
csn=cs + 1/3*dt*(K1+2*K2);


end

function [RHS] = CHRHS(uf,vf,phi,psi,cM)
global D2adv D2GradX D2GradY
Difx=cM.*D2GradX(psi);Dify=cM.*D2GradY(psi);
MLa_Psi=D2GradX(Difx)+D2GradY(Dify);    %M先内乘，再求散度

Fphi=D2adv(uf,vf,uf,vf,phi,phi);%恢复低阶时间格式，但用二阶时间推进
RHS= -Fphi + MLa_Psi;
end

