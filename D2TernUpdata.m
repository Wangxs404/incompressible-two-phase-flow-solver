function [rho,rhoL,rhoLL,mu]=D2TernUpdata(phi,phiL,phiLL)
global  rhoIni muIni
%约束phi/phiL/phiLL在（0,1）
%并由phi更新密度和粘度
range1 = (phi.A > 1);  range2 = (phi.B > 1) ;  range3 = (phi.C > 1);
range4 = (phi.A < 0);  range5 = (phi.B < 0) ;  range6 = (phi.C < 0);
phi.A(range1) = 1 ;   phi.B(range2) = 1 ; phi.C(range3) = 1;
phi.A(range4) = 0 ;   phi.B(range5) = 0 ; phi.C(range6) = 0;

phi.A=D2set_BCNeu(phi.A); phi.B=D2set_BCNeu(phi.B); phi.C=D2set_BCNeu(phi.C); % 默认无润湿性
phiL.A=D2set_BCNeu(phiL.A); phiL.B=D2set_BCNeu(phiL.B); phiL.C=D2set_BCNeu(phiL.C);
phiLL.A=D2set_BCNeu(phiLL.A); phiLL.B=D2set_BCNeu(phiLL.B); phiLL.C=D2set_BCNeu(phiLL.C);

rho=rhoIni.A*phi.A + rhoIni.B*phi.B + rhoIni.C*phi.C;
rhoL=rhoIni.A*phiL.A + rhoIni.B*phiL.B + rhoIni.C*phiL.C;
rhoLL=rhoIni.A*phiLL.A + rhoIni.B*phiLL.B + rhoIni.C*phiLL.C;

mu=muIni.A*phi.A+muIni.B*phi.B+muIni.C*phi.C;   %只需要mu（n+1）
mu=D2set_BCNeu(mu); %rhieChow 需要

end
