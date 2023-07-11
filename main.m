%                               说    明
%  本程序基于同位网格RhieChow插值，常系数投影算法，采用有限体积法+WENO+相场模型求解不可压两相流
%  Copy right by Wangxs

clear ; close all ; warning off
%=================================Global Variable============================
global imax imin jmax jmin nx ny Coord_x Coord_y dx dy dxi dyi dt Lx Ly Swith_Mcc
global Convection_Scheme refere_Pressure  velocity_BCtype  theta FsTension  u_top u_bottom
global Dim initial_type dat_Freq  Dynamic_Draw mat_Freq Res_Freq RhieChow_or_Not How_Many_Phase Datpath

%==================================划分网格&时间步================================
Dim=1;   
Lx=1*Dim    ;  Ly=2*Dim;       nx=64;    ny=128; 
dt=5e-5;    istep_max=ceil(3/dt);  theta=90;  % 默认润湿角
FsTension="ON";    %控制开关Fs
Wetting="Off"; Thermal="Off";  Surfactant="Off";      % 选择多物理场（"ON"）
Electronic="Off"; Magnetic="Off";
Res_Freq=20  ;  dat_Freq=1000;  mat_Freq=10000;  Dynamic_Draw='ON'; % Debug设置
[imin,imax,jmin,jmax,Coord_x,Coord_y,dx,dy,dxi,dyi] = D2Mesh(Lx,Ly);

%===============================对流格式与界面选择=================
How_Many_Phase=2;                      %兼顾单相、两相、三相流动；
Convection_Scheme=1;                   %对流格式  1:WENO 2:Quick 3:Centre
refere_Pressure=4;                     %压力零点  1/2/3/4/5，左下逆时针四对角-体心
velocity_BCtype=1;  % Check“D2set_BC’   %1:NS slip EW Solid_NoSlip; 2:Solid Bc ; 3:NS slip EW Open_NoSlip ; 4:NS slip EW Period 5：All Period
initial_type=1;     % See3"D2Initial_phi'  %1.Bubble_Rising;2.Dam_Break;3.Zalesak's Disk;4.Rayleigh Taylor------
Swith_Mcc="Off";    % See‘D2Calcu_MccSu' %修正通量 1:ON
RhieChow_or_Not=1;  % See'D2Eq_Mome_S2_RhieChow' %1.Centre;9. Full RhieChow-LJY;11. Full RhieChow-HZY
property=2;         % See'D2Property'    %材料参数 (详见‘Property’)
run D2Property      %===指定材料、界面、多物理场参数（要求位于IniPhi之前）

%==================================将数据存储至上级目录=============
name="RB1000-Dense2";        % 在上级目录创建文件夹存储数据
Matpath=strcat("V1-",name,"-Mat");
Datpath=strcat("V1-",name,"-Data");
mkdir('../',strcat("V1-",name,"-Data"))  
mkdir('../',strcat("V1-",name,"-Mat"))

%==============================前处理==============================
[phi,phiL,u_top,u_bottom]=D2Initial_phi(initial_type,0,0) ; %===相界面&平板速度初始化=====
%==============================后处理==============================
ResultDeal="RB10";  %后处理命令包括：{'RB10' 'RB1000' 'laplace' 'MassLeaf' 'ShearDrop' ...
%                                          'LDC100/1000/3200/5000' 'Buoyancy'  'RT'} 

istep=0;                                 %===为便于续算，需将初始istep置于“SpecialSet”前
 %=====边界元位置向量========算子函数句柄=============初始化矩阵大小====
 %====续算/Cases个性化=====绘图、输出/保存初始数据%===指定无量纲量并输出
        run D2LocEleBc ;   run D2HandleFunc  ;  run D2IniMatrix            
        run D2SpecialSet;  run D2PostIni     ;  run D2DimLess ;  
        run D2ConstMatrix


% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  开始闭环运算  $$$$$$$$$$$$$$$$$$$$$$$$$$$
close %关闭前处理图窗

Residual=zeros(istep_max,6);R_limit=100;  time_start=clock ; A_timeStart=datestr(now,'HH:MM mmm.dd');
while (istep < istep_max)
    %=========================Set Boudary Condition==================

    [u,v]=D2set_BC(u,v)    ;    [p]=D2set_BCNeu(p);  
    phi = D2set_BCNeu(phi);      % 默认无润湿性
%     if Wetting=="ON" ; phi = D2set_BCWet(phi,phiL,sigma);end %若开启润湿性，则修改φ_BC
 
    %==============================Phase Field================
     phiLL=phiL; 
     [phin,phiLn,psi,M_change] = D2Eq_PhaseLiuTVD3(uf,vf,phi);
     phi_R=norm(phin-phi)/(nx*ny) ; phiL=phiLn; phi=phin;                %计算残差 ；赋值phi & phiL
    [rho,rhoL,rhoLL,mu,phi,phiL,phiLL]=D2Updata_RhoMiu(phi,phiL,phiLL,sigma,Wetting); %更新ρ、μ，同时约束φ[-1,1]
    [Mccx,Mccy]=D2Calcu_Mcc(rhoL,u,v,M_change,psi);                      %-修正质量通量Mcc

    %==============================MultiPhysics======================
%      if Thermal   =="ON"; [T,sigma,T_R]=D2Eq_Thermal(T,phi,rho,uf,vf) ;end   
%      if Surfactant=="ON"; [cs,cpsi,sigma,cs_R]=D2Eq_Surfactant(cs,phi,sigma,uf,vf) ;end      
%      if Electronic=="ON"; [Fe,PhiE,ChargeQ]= D2Eq_Electronic(phi) ;end 
%      if Magnetic  =="ON"; [Fm,PhiM,H,Emodulus]= D2Eq_Magnetic(phi); end % 完成度60%
        
    %==============================NS Equation=========================
    [un,vn,uLn,vLn,pn,pL,ppien] = D2Eq_Mome(u,v,uL,vL,uf,vf,ufL,vfL,p,pL,mu,rho,rhoL,rhoLL, Mccx, Mccy,psi,phi,sigma,Fe,Fm);
    u_R = norm(un-u)/(nx*ny);    v_R=norm(vn-v)/(nx*ny);  p_R = norm(pn-p)/(nx*ny);
    u=un;   v=vn;  uL=uLn; vL=vLn;  p=pn;   ppie=ppien;
    [uf,vf,ufL,vfL]=D2Uf_inte(u,v,uL,vL);                   %调用插值函数计算Uf .二选一

    %==================================完成单次循环============================
    %% 动态显示残差
    istep=istep+1 ; Residual(istep,1)=u_R ;  Residual(istep,2)=v_R ; Residual(istep,3)=phi_R ;  %计算残差
    if   mod(istep,Res_Freq) ==  0
        fprintf('\n 迭代次数  %s / %s\n  phi_R=   %s\n  u_R  =   %s\n  v_R  =   %s  \n',...
            num2str(istep),num2str(istep_max),num2str(phi_R),num2str(u_R),num2str(v_R));      
        if  Thermal=="ON"
            Residual(istep,4)=T_R ;
            fprintf('  T_R  =   %s\n',num2str(T_R));
        end
        if  Surfactant=="ON"
            Residual(istep,5)=cs_R ;
            fprintf('  cs_R =   %s\n',num2str(cs_R));
        end
    elseif (u_R>R_limit)||(v_R>R_limit)
        fprintf('请再检查一下吧 >o< \n');
        break;
    end
    %% Output and Draw
%=====输出数据到.dat======
if  mod(istep,dat_Freq) ==  0
    var1=phi(imin:imax,jmin:jmax); %1.指定输出变量及范围
    var2=psi(imin:imax,jmin:jmax);
    var3=u(imin:imax,jmin:jmax);
    var4=v(imin:imax,jmin:jmax);
    var5=p(imin:imax,jmin:jmax);
    var6=T(imin:imax,jmin:jmax);
    var7=cs(imin:imax,jmin:jmax);
    var8=cpsi(imin:imax,jmin:jmax);
    var9=sigma(imin:imax,jmin:jmax);
    var10=PhiE(imin:imax,jmin:jmax);
    var11=ChargeQ(imin:imax,jmin:jmax);
    var12=PhiM(imin:imax,jmin:jmax);
    varName ='phi psi u v p T cs cpsi sigma PhiE ChargeQ PhiM\n';  % 2.指定输出变量显示名;注意同一.m文件中后面的varName会被前面的覆盖
    D2export(istep,xx(:),yy(:),NGX,NGY,varName,var1(:),var2(:),var3(:),var4(:),var5(:),var6(:),var7(:),var8(:),var9(:),var10(:),var11(:),var12(:))%3.输出变量
end
%=====绘制动态图======
switch Dynamic_Draw
    case 'ON'
        if     mod(istep,20) ==  0
            contourf(Coord_x(4:1+nx),Coord_y(4:1+ny),phi(4:1+nx,4:1+ny)',3) %相界面 等高线图
%             Vorticity=D2Gradx_Matrix(v)-D2Grady_Matrix(u);
%             contourf(Coord_x(4:1+nx),Coord_y(4:1+ny),Vorticity(4:1+nx,4:1+ny)',20)   %涡量图
            title('Moving phi');axis equal;
            drawnow
        end
end
%=====保存数据到Workspace====
if   mod(istep,mat_Freq) ==  0   %保存数据到Workspace
    filename=strcat('D2Data', num2str(istep));
    MatData=char(strcat('../',Matpath,'/',filename));%将路径和变量名全部拼接成char类型，再save即可
    save (MatData);
end
if   istep==100   %保存数据到Workspace
    time_end100=clock;  runtime100=etime(time_end100,time_start)/60;
    PredictTime=runtime100*istep_max/100;
    fprintf(' PredictTime = %s \n',num2str(PredictTime));
end

end

time_end=clock  ; A_timeEnd=datestr(now,'HH:MM mmm.dd');
A_runtime_s=etime(time_end,time_start); A_runtime_min=etime(time_end,time_start)/60;
save(char(strcat('../',Matpath,'/','D2Data_Multy')))
% run D2DataQuantify
% fprintf(' EndStep= %s \n ',num2str(istep));
% run D2VolLeaf

% if istep==istep_max
% D2ResultDeal(ResultDeal); %后处理命令包括：{'RB10' 'RB1000' 'laplace''MassLeaf'...}
% end



