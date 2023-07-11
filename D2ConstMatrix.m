global Au Av au av Ru RuT Rv RvT
global upsilon0 LL LU LAu LAuT LAv LAvT
global Rus RuTs Rvs RvTs LLs LUs L
%常系数矩阵A
[Au,au]= D2Matrix_Au(ones(nx+6,ny+6),upsilon0*ones(nx+6,ny+6));
[Av,av]= D2Matrix_Av(ones(nx+6,ny+6),upsilon0*ones(nx+6,ny+6));
Ru=chol(Au);RuT=Ru';
Rv=chol(Av);RvT=Rv';
% [RuT,Ru]=lu(Au);[RvT,Rv]=lu(Av);
Rus=sparse(Ru);RuTs=sparse(RuT);
Rvs=sparse(Rv);RvTs=sparse(RvT);
LAu = ichol(Au);LAuT=LAu';
LAv = ichol(Av);LAvT=LAv';
%常系数泊松矩阵
[L] = D2Matrix_Laplace(ones(nx+6,ny+6));
[LL,LU]=lu(L);
LLs=sparse(LL);LUs=sparse(LU);