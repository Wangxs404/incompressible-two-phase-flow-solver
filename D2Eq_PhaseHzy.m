function [ phin,phiLn,psi,M_change] = D2Eq_PhaseHzy(uf,vf,phi)
global  M epsilon  dt  How_Many_Phase sigma0
%==============================离散求解CH方程===============================
lamda=(3*sigma0*epsilon) / (2*sqrt(2));%此处sigma|lamda不受外场影响，只计算原始表面张力（在它处修正）
psi=lamda*(1/epsilon^2 * phi .*(phi.^2-1) -D2La_Oper(phi));

% lamda=(3*sigma*epsilon) / (2*sqrt(2));
% psi=lamda.*(1/epsilon^2 * phi .*(phi.^2-1) -D2La_Oper(phi));

% H1=(3*sigma0)/(2*sqrt(2)*epsilon);H2=3*epsilon*sigma0/(2*sqrt(2));%Huang化简结果
% psi= H1*(phi.^3-phi) - H2*D2La_Oper(phi);
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

