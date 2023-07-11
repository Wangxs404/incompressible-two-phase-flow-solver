function[b]=D2Calcu_bConst(u,uL,m,g,D2GradFunc,rho,rhoL,rhoLL,mu,uf,vf,ufL,vfL,psi,phi,p,pL,D2GradFunc2,sigma,Dt,Fe,Fm)
global ds  imin imax jmin jmax  nx ny dt D2adv  epsilon FsTension How_Many_Phase upsilon0 rho0
%Ready for Calcu
uR=2*u-uL; pR=2*p-pL;

%=============================计算buConst======================================
BuTime= ds*(2*rhoL.*u-0.5*rhoLL.*uL)./rho;    %未考虑不合理密度商
bu_time=reshape(BuTime(imin:imax,jmin:jmax),nx*ny,1);

BuAdv=-dt*ds*D2adv(uf,vf,ufL,vfL,m,m)./rho; %离散层面-守恒形式 U▽·(ρU）
% BuAdv=-dt*ds*D2adv(uf,vf,ufL,vfL,uR,uR);
bu_adv=reshape(BuAdv(imin:imax,jmin:jmax),nx*ny,1);%基于WENO和面心速度计算对流项，并时空积分

BuVisc1=  dt*ds*D2DiverGrad(mu,uR)./rho;  %黏性项1
bu_visc1=reshape(BuVisc1(imin:imax,jmin:jmax),nx*ny,1);  

BuGradP1= -dt*ds*D2GradFunc(pR)./rho;  %压力梯度项1
bu_p1=reshape(BuGradP1(imin:imax,jmin:jmax),nx*ny,1);

BuVisc2=  -dt*ds*upsilon0*D2La_Oper(uR);  %黏性项2！
bu_visc2=reshape(BuVisc2(imin:imax,jmin:jmax),nx*ny,1); 

BuGradP2=  dt*ds*D2GradFunc(p-pL)/rho0;  %压力梯度项2
bu_p2=reshape(BuGradP2(imin:imax,jmin:jmax),nx*ny,1);

bu_dif= reshape(dt*ds*Dt,nx*ny,1);   %黏性项转置

BuFg= dt*ds*g*ones(nx+6,ny+6);  %重力项
bu_Fg=reshape(BuFg(imin:imax,jmin:jmax),nx*ny,1);

%% 表面张力 (通过 FsTension ~="ON"，关闭表面张力-sigma置零)
if FsTension ~="ON"; sigma= 0*sigma ;end 
    
if How_Many_Phase~=3  
    Fs=(0.75*sqrt(2)*epsilon)*Calcu_Fs(sigma,phi,psi,D2GradFunc,D2GradFunc2);%两相的表面张力计算
    BuFs= dt*ds*Fs./rho;
    bu_Fs=reshape(BuFs(imin:imax,jmin:jmax),nx*ny,1);
elseif How_Many_Phase==3
    Fs.A=(6*sqrt(2)*epsilon)*Calcu_Fs(sigma.A,phi.A,psi.A,D2GradFunc,D2GradFunc2); %三相的表面张力计算
    Fs.B=(6*sqrt(2)*epsilon)*Calcu_Fs(sigma.B,phi.B,psi.B,D2GradFunc,D2GradFunc2);
    Fs.C=(6*sqrt(2)*epsilon)*Calcu_Fs(sigma.C,phi.C,psi.C,D2GradFunc,D2GradFunc2);
    FsTern= Fs.A+Fs.B+Fs.C;
    BuFs= dt*ds*FsTern./rho;
    bu_Fs=reshape(BuFs(imin:imax,jmin:jmax),nx*ny,1);
end

%% 电场力（已在IniM时初始化为0，故不开电场方程时Fe=0）
BuFe= dt*ds*Fe./rho;
bu_Fe=reshape(BuFe(imin:imax,jmin:jmax),nx*ny,1);

%% 磁场力（已在IniM时初始化为0，故不开磁场方程时Fm=0）
BuFm= dt*ds*Fm./rho;
bu_Fm=reshape(BuFm(imin:imax,jmin:jmax),nx*ny,1);

b=  bu_time + bu_adv + bu_p1 + bu_p2 +bu_visc1+ bu_visc2 + bu_dif...
   +bu_Fg + bu_Fs + bu_Fe + bu_Fm;

end

function [Fs]=Calcu_Fs(sigma,phi,psi,D2GradFunc,D2GradFunc2)
global epsilon
Fs1= psi.*D2GradFunc(phi).*(sigma/epsilon^2); % 法向Fs
Fs2= (D2GradFunc(phi).^2 + D2GradFunc2(phi).^2) .* D2GradFunc(sigma);% 切向Fs1
Fs3= -(D2GradFunc(sigma).*D2GradFunc(phi)+D2GradFunc2(sigma).*D2GradFunc2(phi)).*D2GradFunc(phi);%切向Fs2

Fs= Fs1+Fs2+Fs3;
end
