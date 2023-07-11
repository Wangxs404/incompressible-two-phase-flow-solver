function [cs,cpsi,sigma,cs_R]=D2Eq_Surfactant(cs,phi,sigma,uf,vf)
global nx ny BetaS
[ csn,cpsi] = D2Eq_PhaseSurfactant(uf,vf,cs,phi);
cs_R=norm(csn-cs)/(nx*ny) ; cs=csn;

%cs调控表面张力------------
CoffWeaken=1+BetaS*log(1-cs);   CoffWeaken(CoffWeaken<0.5)=0.5;  %弱化倍数，最低0.5
sigma= CoffWeaken.*sigma;

end