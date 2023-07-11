function[sigmaTern,epsilonPie]=D2TernConst(sigma_12,  sigma_13, sigma_23,epsilon)
global nx ny
sigma1= sigma_12 + sigma_13 - sigma_23;
sigma2= sigma_12 - sigma_13 + sigma_23;
sigma3=-sigma_12 + sigma_13 + sigma_23;
sigmaT=3/(1/sigma1+1/sigma2+1/sigma3);
epsilonPie=epsilon/(4*sqrt(2)); %或许并不需要
% 将三相界面的sigma转为结构体矩阵，便于今后引入表面张力控制
sigmaTern.A=sigma1*ones(nx+6,ny+6);
sigmaTern.B=sigma2*ones(nx+6,ny+6);
sigmaTern.C=sigma3*ones(nx+6,ny+6);
sigmaTern.T=sigmaT*ones(nx+6,ny+6);
end
