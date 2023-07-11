function [ u,v,uL,vL,p,pL,ppie] = D2Eq_Mome(u,v,uL,vL,uf,vf,ufL,vfL,p,pL,mu,rho,rhoL,rhoLL, Mccx, Mccy,psi,phi,sigma,Fe,Fm)
%%  �˺�������ͶӰ����������ָʾ����������ɣ�����ȫ�ֱ��������Ӻ���������
global dt D2GradX D2GradY rho0
%%
%Step1��FVM����Ԥ�ⲽ�ٶȳ�[U_star��V_star]===============================
[u_star,v_star]=D2Eq_Mome_S1Const(u,v,uL,vL,uf,vf,ufL,vfL,p,pL,mu,rho,rhoL,rhoLL, Mccx, Mccy,psi,phi,sigma,Fe,Fm);
uL=u;   vL=v;           %Ϊ����BDF���ף��ضϴ���һ��u��v

%%
%Step2��FVM��Poison���̣�����Ppie=========================================
[u_star,v_star]=D2set_BC(u_star,v_star);[p]=D2set_BCNeu(p); 
[ppie] = D2RhieChow2(u_star,v_star,rho,mu,p);%RC��ֵ����������ѡ��0��dt��dCap��
[ppie] = D2set_BCNeu(ppie);%����һ��BC������������ȫΪ0������´�ѭ����P�޷�ʹ��
pL=p;       %��ȡpL�����ڼ���Uf
p=p+ppie;   %���ڲ���ѹ���㷨ʱ����Ppie�������uv��ͬ������p�����´�ѭ����Step1

%%
% Step3������U_n+1========================================================
[ppie] = D2set_BCNeu(ppie);
% u=u_star-(2*dt/3)./rho.*D2GradX(ppie) ; v=v_star-(2*dt/3)./rho.*D2GradY(ppie);
u=u_star-(2*dt)/(3*rho0)*D2GradX(ppie) ; v=v_star-(2*dt)/(3*rho0)*D2GradY(ppie);
end

