%                               说    明
%     本程序基于同位网格投影算法，采用有限体积法+相场模型求解不可压两相流
clear ; close all ; warning off
%========================= =====Global Variable============================
global imax imin jmax jmin kmin kmax nx ny Coord_x Coord_y Coord_z dx dy dz dxi dyi dzi dt Lx Ly Swith_SM
global Convection_Scheme refere_Pressure  velocity_BCtype Swith_Force rho_Heavy rho_Light dv
global Dim initial_type dat_Freq Dynamic_Draw mat_Freq Res_Freq RhieChow_or_Not How_Many_Phase
global Lz nz 
%==================================划分网格&时间步================================gx gy
Dim=1;
Lx=1*Dim  ;  Ly=1*Dim;  Lz=1*Dim;
nx=50;   ny=50;  nz=50; %nx must equal to ny ???
dt=2.5e-5;  istep_max=ceil(1/dt);  
[imin,imax,jmin,jmax,kmin,kmax,Coord_x ,Coord_y ,Coord_z,dx,dy,dz,dxi,dyi,dzi,dv] = D3Mesh(Lx,Ly,Lz);
%===============================对流格式与界面选择=================
How_Many_Phase=2;
Convection_Scheme=1;                   %对流格式  1:WENO 2:Quick 3:Centre
refere_Pressure=4;                     %压力零点  1/2/3/4/5，左下逆时针四对角-体心
Swith_Force=1;                         %1:G+Fs ; 2: Only G
velocity_BCtype=2;  %Check‘set_BC’   %1:NS slip EW Solid_NoSlip; 2:Solid Bc ; 3:NS slip EW Open_NoSlip ; 4:NS slip EW Period 5：All Period
initial_type=1;     %See‘D3Initial_phi' %1.Bubble_Rising;2.Dam_Break;3.Zalesak's Disk;4.Rayleigh Taylor------
Swith_SM=2;         %See‘Calcu_MccSu' %修正格式 1:Su  2:Mcc-Dh
RhieChow_or_Not=11;  %See‘Eq_Mome_S2_RhieChow' %1.Centre;9. Full RhieChow-LJY;11. Full RhieChow-HZY
property=2;         %See‘D3Property'    %材料参数 (详见‘Property’)
run D3Property                           %===指定材料和界面参数（要求位于IniPhi之前）
%================================
Res_Freq=20  ;  dat_Freq=500;  mat_Freq=1000; Dynamic_Draw='ON'; %===Debug设置===
[phi,phiL]=D3Initial_phi(initial_type); %===相界面初始化=====
istep=0;                                 %===为便于续算，需将初始istep置于“SpecialSet”前
% run DimensionlessData                    %===指定无量纲量并输出，以及计算双HELM_Alpha
run D3HandleFunc                           %===生成全局-算子函数句柄
run D3Def_Ini_Matrix                       %===初始化矩阵大小
run D3SpecialSet                           %===控制续算 计算初始质量以及各算例的个性化设置
run D3Post_InitialPhi                      %===绘制初始界面、输出/保存初始数据
run D3LocationBc   %位置
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  开始闭环运算  $$$$$$$$$$$$$$$$$$$$$$$$$$$
close all %关闭图窗
Residual=zeros(istep_max,3);
R_limit=100;  time_start=clock;A_timeStart=datestr(now,'HH:MM mmm.dd');
while (istep < istep_max)
    %================================设置边条==========================
    [u,v,w,uf,vf,wf]=D3GivenFlow(nx,ny,nz);%给定速度场，检验CH方程
    [u,v,w] = D3set_BC(u,v,w);    [p]=D3set_BCNeu(p);   [phi]=D3set_BCNeu(phi); 
    %================================计算相场方程得[phi(n+1)]================
    phiLL=phiL;                                            %截留phiL，作为phiLL，用于动量方程的计算 
    [phin,phiLn,psi,M_change] = D3Eq_PhaseHzy(uf,vf,wf,phi);
    phi_R=mean(phin,'all')-mean(phi,'all');
    phiL=phiLn; phi=phin;   %计算残差 ；赋值phi & phiL
    [rho,rhoL,rhoLL,mu,phi,phiL,phiLL]=D3Updata_RhoMiu(phi,phiL,phiLL,rho_Light,rho_Heavy,mu_Light,mu_Heavy); %更新ρ、μ，同时约束φ[-1,1]
    %%
%     %==============================计算动量方程得[U V]=========================
   [ un,vn,wn,uLn,vLn,wLn,pn] = D3Eq_Mome(u,v,w,uL,vL,wL,uf,vf,wf,ufL,vfL,wfL,p,mu,rho,rhoL,rhoLL,psi,phi);
    u_R =mean(un,'all')-mean(u,'all');  v_R =mean(vn,'all')-mean(v,'all');  p_R = mean(pn,'all')-mean(p,'all');
    u=un;  v=vn;  w=wn;
    uL=uLn; vL=vLn;  wL=wLn;  p=pn;   
    [uf,vf,wf,ufL,vfL,wfL]= D3Uf_inte(u,v,w,uL,vL,wL) ; %暂用简单线性差分计算面心，该为高阶插值更好
%     [uf,vf,ufL,vfL]=Uf_inte(u,v,uL,vL); 
    %==================================完成单次循环============================
    %% 动态显示残差
    istep=istep+1 ; Residual(istep,1)=u_R ;  Residual(istep,2)=v_R ; Residual(istep,3)=phi_R ;  %计算残差
    if   mod(istep,Res_Freq) ==  0
        fprintf(' 迭代次数  %s / %s\n  phi_R=   %s\n  u_R  =   %s\n  v_R  =   %s  \n\n',...
            num2str(istep),num2str(istep_max),num2str(phi_R),num2str(u_R),num2str(v_R));
    elseif (u_R>R_limit)||(v_R>R_limit)
        fprintf('残差过大不收敛');
        break;
    end
    %% Output and Draw
% %=====输出数据到.dat======
if  mod(istep,dat_Freq) ==  0 || istep==1
    var1=phi(4:imax+1,4:jmax+1,4:kmax+1); %1.指定输出变量及范围
    var2=psi(4:imax+1,4:jmax+1,4:kmax+1);
    var3=u(4:imax+1,4:jmax+1,4:kmax+1);
    var4=v(4:imax+1,4:jmax+1,4:kmax+1);
    var5=w(4:imax+1,4:jmax+1,4:kmax+1);
    var6=p(4:imax+1,4:jmax+1,4:kmax+1);
    varName ='phi psi u v w p\n';  % 2.指定输出变量显示名;注意同一.m文件中后面的varName会被前面的覆盖
    D3export(istep,xx(:),yy(:),zz(:),NGX,NGY,NGZ,varName,var1(:),var2(:),var3(:),var4(:),var5(:),var6(:))%3.输出变量
end
%=====绘制动态图======
switch Dynamic_Draw
    case 'ON'
        if     mod(istep,50) ==  0
            figure(3)  %绘制三维等值面图
            D3show(phi);
            grid minor;
            title('Moving phi face');
            drawnow ;  %强切窗口
            clf

%             figure(4)  %沿z方向切片，绘制等值线图
%             D3slice(phi,[0.3,0.5,0.6])
%             grid minor;
%             title('Moving phi slice');
%             drawnow
        end
end
%=====保存数据到Workspace====
if   mod(istep,mat_Freq) ==  0   %保存数据到Workspace
    filename=strcat('D3Data', num2str(istep));
    save (filename);
end
if   istep==100   
    time_end100=clock;  runtime100=etime(time_end100,time_start)/60;
    PredictTime=runtime100*istep_max/100;
    fprintf(' PredictTime = %s \n',num2str(PredictTime));
end

end
time_end=clock  ; A_timeEnd=datestr(now,'HH:MM mmm.dd');
A_runtime_s=etime(time_end,time_start); A_runtime_min=etime(time_end,time_start)/60;
save D3Data_Multy
% run DataQuantify
EndStep=istep;
fprintf(' EndStep= %s \n ',num2str(istep));
% fprintf('MassLeaf= %s %s\n',num2str(MassLeaf(end)*100),'%');
% run Post_ShearDrop
