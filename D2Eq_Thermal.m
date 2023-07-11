function [T,sigma,T_R] = D2Eq_Thermal(T,phi,rho,uf,vf)
global nx ny

T=D2set_BCT(T);  %温度场 |考虑温度场时sigma需为全局矩阵，否则与RhieChow不兼容
[sigma,Tn]=D2Theimal(phi,rho,uf,vf,T);
T_R=norm(Tn-T)/(nx*ny) ;T=Tn;

end