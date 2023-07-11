function [ u,v,uL,vL,p,pL,ppie] = D2Eq_Mome(u,v,uL,vL,uf,vf,ufL,vfL,p,pL,mu,rho,rhoL,rhoLL, Mccx, Mccy,psi,phi,sigma,Fe,Fm)
%%  此函数汇总投影法三步，仅指示输入输出即可，所用全局变量均在子函数内声明
global dt D2GradX D2GradY rho0
%%
%Step1：FVM计算预测步速度场[U_star，V_star]===============================
[u_star,v_star]=D2Eq_Mome_S1Const(u,v,uL,vL,uf,vf,ufL,vfL,p,pL,mu,rho,rhoL,rhoLL, Mccx, Mccy,psi,phi,sigma,Fe,Fm);
uL=u;   vL=v;           %为引入BDF二阶，截断储存一次u，v

%%
%Step2：FVM解Poison方程，计算Ppie=========================================
[u_star,v_star]=D2set_BC(u_star,v_star);[p]=D2set_BCNeu(p); 
[ppie] = D2RhieChow2(u_star,v_star,rho,mu,p);%RC插值，包含三种选择（0，dt，dCap）
[ppie] = D2set_BCNeu(ppie);%补充一次BC，否则外三层全为0，造成下次循环的P无法使用
pL=p;       %截取pL，用于计算Uf
p=p+ppie;   %基于补充压力算法时，由Ppie计算更新uv，同步更新p用于下次循环的Step1

%%
% Step3：计算U_n+1========================================================
[ppie] = D2set_BCNeu(ppie);
% u=u_star-(2*dt/3)./rho.*D2GradX(ppie) ; v=v_star-(2*dt/3)./rho.*D2GradY(ppie);
u=u_star-(2*dt)/(3*rho0)*D2GradX(ppie) ; v=v_star-(2*dt)/(3*rho0)*D2GradY(ppie);
end

