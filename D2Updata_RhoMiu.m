function[rho,rhoL,rhoLL,mu,phi,phiL,phiLL]=D2Updata_RhoMiu(phi,phiL,phiLL,sigma,Wetting)
global Reconstrct  rho_Heavy rho_Light  mu_Light  mu_Heavy
% 根据CH使用的是[0,1]或[-1,1] 选择更新密度粘度的方式
%约束phi/phiL/phiLL在（-1,1）
%并由phi更新密度和粘度
range1 = (phi > 1);  range2 = (phi < -1) ;   phi(range1) = 1;    phi(range2) = -1;
range3 = (phiL > 1); range4 = (phiL < -1);   phiL(range3) = 1;    phiL(range4) = -1;
range5 = (phiLL > 1); range6 = (phiLL < -1); phiLL(range5) = 1;    phiLL(range6) = -1;

%外插，同步影响以下密度和粘度
[phi]=D2set_BCNeu(phi);  % 默认无润湿性
[phiL]=D2set_BCNeu(phiL);
[phiLL]=D2set_BCNeu(phiLL);
if Wetting=="ON"
    [phi] = D2set_BCWet(phi,phiL,sigma);
    [phiL] = D2set_BCWet(phiL,phiLL,sigma);
    [phiLL] = D2set_BCWet(phiLL,phiLL,sigma);
end


rho=Reconstrct(rho_Light,rho_Heavy,phi);
rhoL=Reconstrct(rho_Light,rho_Heavy,phiL);
rhoLL=Reconstrct(rho_Light,rho_Heavy,phiLL);
mu=Reconstrct(mu_Light,mu_Heavy,phi);  %只需要mu（n+1）

[mu]=D2set_BCNeu(mu); %rhieChow 需要
end

