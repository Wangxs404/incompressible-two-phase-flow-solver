function [psi]=D2TernPsi(phi,sigmaTern)
global epsilon 
psi.A= phi.A.*(phi.A-1).*(phi.A-0.5)-(sigmaTern.T)./(sigmaTern.A).*(phi.A.*phi.B.*phi.C)-epsilon^2*D2La_Oper(phi.A);
psi.B= phi.B.*(phi.B-1).*(phi.B-0.5)-(sigmaTern.T)./(sigmaTern.B).*(phi.A.*phi.B.*phi.C)-epsilon^2*D2La_Oper(phi.B);
psi.C= 0-psi.A-psi.B;
end
