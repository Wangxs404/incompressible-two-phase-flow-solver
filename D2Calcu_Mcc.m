function [Mccx,Mccy]=D2Calcu_Mcc(rhoL,u,v,M_change,psi)
global  rho_Light rho_Heavy D2GradX D2GradY
%经测试，Mcc几乎对计算没有影响，这不是相场大密度比的关键
psi = D2set_BCNeu(psi);
Mccx=rhoL.*u-(rho_Light-rho_Heavy)/2*M_change.*D2GradX(psi); %修正质量通量Mcc中的U是否需要先使用一次WENO，待商榷
Mccy=rhoL.*v-(rho_Light-rho_Heavy)/2*M_change.*D2GradY(psi);%由于对流项显式化，所以用rhoL

end

% %DHang
% %     SUx=(rho_Light-rho_Heavy)*MLa_Psi.*u;
% %     SUy=(rho_Light-rho_Heavy)*MLa_Psi.*v;
% %Huang
% SUx=(rho_Light-rho_Heavy)*M_change.*(D2GradX(psi).*D2GradX(u)+ D2GradY(psi).* D2GradY(u));
% SUy=(rho_Light-rho_Heavy)*M_change.*(D2GradX(psi).*D2GradX(v)+ D2GradY(psi).* D2GradY(v));