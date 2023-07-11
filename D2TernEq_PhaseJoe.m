function [ phin,phiLn,psi] = D2TernEq_PhaseJoe(uf,vf,phi,sigmaTern)
global  M 
%==============================离散求解CH方程===============================
%Tern Psi
[psi]=D2TernPsi(phi,sigmaTern); %形态需稍加核对
[psi.A]=D2set_BCNeu(psi.A); [psi.B]=D2set_BCNeu(psi.B); [psi.C]=D2set_BCNeu(psi.C); 

%Tern M
M_A=M./sigmaTern.A; M_A=M_A.*abs(phi.A.*(phi.A-1));
M_B=M./sigmaTern.B; M_B=M_B.*abs(phi.B.*(phi.B-1));

%CH
phiLn=phi; %截留phiL
phin.A= SolveCH(phi.A,uf,vf,psi.A,M_A);
phin.B= SolveCH(phi.B,uf,vf,psi.A,M_B);
phin.C= 1-phi.A-phi.B;

end

function [phin]=SolveCH(phi,uf,vf,psi,M) % 更新相分数
global dt
K1=CHRHS(uf,vf,phi,psi,M);
K2=CHRHS(uf,vf,phi+0.75*dt*K1,psi,M);
phin=phi + 1/3*dt*(K1+2*K2);
end

function [RHS] = CHRHS(uf,vf,phi,psi,M)
global D2adv D2GradX D2GradY
Difx=M.*D2GradX(psi);Dify=M.*D2GradY(psi);
MLa_Psi=D2GradX(Difx)+D2GradY(Dify);    %M先内乘，再求散度

Fphi=D2adv(uf,vf,uf,vf,phi,phi);%恢复低阶时间格式，但用二阶时间推进
RHS= -Fphi + MLa_Psi;
end

