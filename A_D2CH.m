%                               说    明
%     本程序采用二阶时间推进，验证CH方程的求解精度
%     Case1为Za圆盘 Case为单涡逆转
%     
clear ; close all ; warning off
%========================= =====Global Variable============================
global imax imin jmax jmin nx ny Coord_x Coord_y dx dy dxi dyi dt Lx Ly 
global Convection_Scheme  Dim initial_type dat_Freq  Dynamic_Draw mat_Freq Res_Freq  How_Many_Phase
%==================================划分网格&时间步================================gx gy
Dim=1;
Lx=1*Dim  ;  Ly=1*Dim;   nx=100;   ny=100; 
[imin,imax,jmin,jmax,Coord_x,Coord_y,dx,dy,dxi,dyi] = D2Mesh(Lx,Ly);

CHCase=1; % Za disk | Reversed Circle
CFL=0.1 ; 
switch CHCase
    case 1
        Pe=800 ; Cn=0.05;  NPeriods=1;% Za时调这里
        M=1e-5 ;  epsilon=Dim*Cn;
        U0=Pe*M/epsilon; sigma=0.01;
        dt=CFL*dx/U0;  istep_max=ceil(NPeriods*(2/(U0*dt)));
    case 2
        Periods=2;  Cn=0.02;      %Reversed Circle时调这里,注释1
        M=(1e-7)*(Cn*Dim/(1/32));
        dt = CFL*Cn ; istep_max=Periods/dt;
end
%===============================对流格式与界面选择=================
How_Many_Phase=2;
Convection_Scheme=1;                   %对流格式  1:WENO 2:Quick 3:Centre
initial_type=3;     %See‘D2Initial_phi' %3.Zalesak's Disk;10.Reversed Circle
%================================
Res_Freq=50  ;  dat_Freq=10;  mat_Freq=1000; Dynamic_Draw='ON'; %===Debug设置===

[phi,phiL]=D2Initial_phi(initial_type) ; %===相界面初始化=====
istep=0;                                 %===为便于续算，需将初始istep置于“SpecialSet”前
run D2HandleFunc                           %算子函数句柄
run D2Def_Ini_Matrix                       %===初始化矩阵大小
run D2Post_InitialPhi                      %===绘制初始界面、输出/保存初始数据
run D2CH_IniUV                             %给定速度场% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  开始闭环运算  $$$$$$$$$$$$$$$$$$$$$$$$$$$
close
Residual=zeros(istep_max,1);
R_limit=100;  time_start=clock;A_timeStart=datestr(now,'HH:MM mmm.dd');
while (istep < istep_max)
switch CHCase
    case 2
    time=istep*dt;
   [u,v] = D2InitialUVRsv(u,v,time,Periods);
   uL=u ; vL=v;
   u=D2set_BCcycle(u) ; v=D2set_BCcycle(v);
   uf=D2fXmean(u);vf=D2fYmean(v);
   ufL=uf ;  vfL=vf; 
end
    %================================计算相场方程得[phi(n+1)]================
    phiLL=phiL;                                            %截留phiL，作为phiLL，用于动量方程的计算 
    [phin,phiLn,psi,M_change] = D2Eq_PhaseHzy(uf,vf,phi);
    phi_R=norm(phin-phi)/(nx*ny) ; 
    phiL=phiLn; phi=phin;   %计算残差 ；赋值phi & phiL
    %==================================完成单次循环============================
    %% 动态显示残差
    istep=istep+1 ;Residual(istep,1)=phi_R ;  %计算残差
    if   mod(istep,Res_Freq) ==  0
        fprintf(' 迭代次数  %s / %s\n  phi_R=   %s\n   \n\n',...
            num2str(istep),num2str(istep_max),num2str(phi_R));
    elseif (phi_R>R_limit)
        fprintf('残差过大不收敛');
        break;
    end
    %% Output and Draw
%=====输出数据到.dat======
if  mod(istep,dat_Freq) ==  0
    var1=phi(4:imax+1,4:jmax+1); %1.指定输出变量及范围
    var2=psi(4:imax+1,4:jmax+1);
    var3=u(4:imax+1,4:jmax+1);
    var4=v(4:imax+1,4:jmax+1);
    varName ='phi psi u v\n';  % 2.指定输出变量显示名;注意同一.m文件中后面的varName会被前面的覆盖
    D2export(istep,xx(:),yy(:),NGX,NGY,varName,var1(:),var2(:),var3(:),var4(:))%3.输出变量
end
%=====绘制动态图======
switch Dynamic_Draw
    case 'ON'
        if     mod(istep,20) ==  0
            contour(Coord_x(4:1+nx),Coord_y(4:1+ny),phi(4:1+nx,4:1+ny)',1) %等高线图
            title('Moving phi');axis equal;
            drawnow
        end
end
%=====保存数据到Workspace====
if   mod(istep,mat_Freq) ==  0   %保存数据到Workspace
    filename=strcat('D2Data', num2str(istep));
    save (filename);
end
if   istep==100   %保存数据到Workspace
    time_end100=clock;  runtime100=etime(time_end100,time_start)/60;
    A_PredictTime=runtime100*istep_max/100;
    fprintf(' PredictTime = %s \n',num2str(A_PredictTime));
end

end
time_end=clock  ; A_timeEnd=datestr(now,'HH:MM mmm.dd');
A_runtime_s=etime(time_end,time_start); A_runtime_min=etime(time_end,time_start)/60;
save D2Data_Multy
fprintf(' EndStep= %s \n ',num2str(istep));

errorphi = sum(abs(phi_initial(:) - phi(:))) / sum(abs(phi_initial(:)));
fprintf('errorphi=%s ',num2str(errorphi))
