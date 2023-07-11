function [MLa_Psi] = D2DiverGrad(M_change,psi)
global  D2GradX D2GradY
Difx=M_change.*D2GradX(psi);Dify=M_change.*D2GradY(psi);
MLa_Psi=D2GradX(Difx)+D2GradY(Dify);    %M先内乘，再求散度
end