%                               ˵    ��
%  ���������ͬλ����RhieChow��ֵ����ϵ��ͶӰ�㷨���������������+WENO+�ೡģ����ⲻ��ѹ������
%  Copy right by Wangxs

clear ; close all ; warning off
%=================================Global Variable============================
global imax imin jmax jmin nx ny Coord_x Coord_y dx dy dxi dyi dt Lx Ly Swith_Mcc
global Convection_Scheme refere_Pressure  velocity_BCtype  theta FsTension  u_top u_bottom
global Dim initial_type dat_Freq  Dynamic_Draw mat_Freq Res_Freq RhieChow_or_Not How_Many_Phase Datpath

%==================================��������&ʱ�䲽================================
Dim=1;   
Lx=1*Dim    ;  Ly=2*Dim;       nx=64;    ny=128; 
dt=5e-5;    istep_max=ceil(3/dt);  theta=90;  % Ĭ����ʪ��
FsTension="ON";    %���ƿ���Fs
Wetting="Off"; Thermal="Off";  Surfactant="Off";      % ѡ���������"ON"��
Electronic="Off"; Magnetic="Off";
Res_Freq=20  ;  dat_Freq=1000;  mat_Freq=10000;  Dynamic_Draw='ON'; % Debug����
[imin,imax,jmin,jmax,Coord_x,Coord_y,dx,dy,dxi,dyi] = D2Mesh(Lx,Ly);

%===============================������ʽ�����ѡ��=================
How_Many_Phase=2;                      %��˵��ࡢ���ࡢ����������
Convection_Scheme=1;                   %������ʽ  1:WENO 2:Quick 3:Centre
refere_Pressure=4;                     %ѹ�����  1/2/3/4/5��������ʱ���ĶԽ�-����
velocity_BCtype=1;  % Check��D2set_BC��   %1:NS slip EW Solid_NoSlip; 2:Solid Bc ; 3:NS slip EW Open_NoSlip ; 4:NS slip EW Period 5��All Period
initial_type=1;     % See3"D2Initial_phi'  %1.Bubble_Rising;2.Dam_Break;3.Zalesak's Disk;4.Rayleigh Taylor------
Swith_Mcc="Off";    % See��D2Calcu_MccSu' %����ͨ�� 1:ON
RhieChow_or_Not=1;  % See'D2Eq_Mome_S2_RhieChow' %1.Centre;9. Full RhieChow-LJY;11. Full RhieChow-HZY
property=2;         % See'D2Property'    %���ϲ��� (�����Property��)
run D2Property      %===ָ�����ϡ����桢������������Ҫ��λ��IniPhi֮ǰ��

%==================================�����ݴ洢���ϼ�Ŀ¼=============
name="RB1000-Dense2";        % ���ϼ�Ŀ¼�����ļ��д洢����
Matpath=strcat("V1-",name,"-Mat");
Datpath=strcat("V1-",name,"-Data");
mkdir('../',strcat("V1-",name,"-Data"))  
mkdir('../',strcat("V1-",name,"-Mat"))

%==============================ǰ����==============================
[phi,phiL,u_top,u_bottom]=D2Initial_phi(initial_type,0,0) ; %===�����&ƽ���ٶȳ�ʼ��=====
%==============================����==============================
ResultDeal="RB10";  %�������������{'RB10' 'RB1000' 'laplace' 'MassLeaf' 'ShearDrop' ...
%                                          'LDC100/1000/3200/5000' 'Buoyancy'  'RT'} 

istep=0;                                 %===Ϊ�������㣬�轫��ʼistep���ڡ�SpecialSet��ǰ
 %=====�߽�Ԫλ������========���Ӻ������=============��ʼ�������С====
 %====����/Cases���Ի�=====��ͼ�����/�����ʼ����%===ָ���������������
        run D2LocEleBc ;   run D2HandleFunc  ;  run D2IniMatrix            
        run D2SpecialSet;  run D2PostIni     ;  run D2DimLess ;  
        run D2ConstMatrix


% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  ��ʼ�ջ�����  $$$$$$$$$$$$$$$$$$$$$$$$$$$
close %�ر�ǰ����ͼ��

Residual=zeros(istep_max,6);R_limit=100;  time_start=clock ; A_timeStart=datestr(now,'HH:MM mmm.dd');
while (istep < istep_max)
    %=========================Set Boudary Condition==================

    [u,v]=D2set_BC(u,v)    ;    [p]=D2set_BCNeu(p);  
    phi = D2set_BCNeu(phi);      % Ĭ������ʪ��
%     if Wetting=="ON" ; phi = D2set_BCWet(phi,phiL,sigma);end %��������ʪ�ԣ����޸Ħ�_BC
 
    %==============================Phase Field================
     phiLL=phiL; 
     [phin,phiLn,psi,M_change] = D2Eq_PhaseLiuTVD3(uf,vf,phi);
     phi_R=norm(phin-phi)/(nx*ny) ; phiL=phiLn; phi=phin;                %����в� ����ֵphi & phiL
    [rho,rhoL,rhoLL,mu,phi,phiL,phiLL]=D2Updata_RhoMiu(phi,phiL,phiLL,sigma,Wetting); %���¦ѡ��̣�ͬʱԼ����[-1,1]
    [Mccx,Mccy]=D2Calcu_Mcc(rhoL,u,v,M_change,psi);                      %-��������ͨ��Mcc

    %==============================MultiPhysics======================
%      if Thermal   =="ON"; [T,sigma,T_R]=D2Eq_Thermal(T,phi,rho,uf,vf) ;end   
%      if Surfactant=="ON"; [cs,cpsi,sigma,cs_R]=D2Eq_Surfactant(cs,phi,sigma,uf,vf) ;end      
%      if Electronic=="ON"; [Fe,PhiE,ChargeQ]= D2Eq_Electronic(phi) ;end 
%      if Magnetic  =="ON"; [Fm,PhiM,H,Emodulus]= D2Eq_Magnetic(phi); end % ��ɶ�60%
        
    %==============================NS Equation=========================
    [un,vn,uLn,vLn,pn,pL,ppien] = D2Eq_Mome(u,v,uL,vL,uf,vf,ufL,vfL,p,pL,mu,rho,rhoL,rhoLL, Mccx, Mccy,psi,phi,sigma,Fe,Fm);
    u_R = norm(un-u)/(nx*ny);    v_R=norm(vn-v)/(nx*ny);  p_R = norm(pn-p)/(nx*ny);
    u=un;   v=vn;  uL=uLn; vL=vLn;  p=pn;   ppie=ppien;
    [uf,vf,ufL,vfL]=D2Uf_inte(u,v,uL,vL);                   %���ò�ֵ��������Uf .��ѡһ

    %==================================��ɵ���ѭ��============================
    %% ��̬��ʾ�в�
    istep=istep+1 ; Residual(istep,1)=u_R ;  Residual(istep,2)=v_R ; Residual(istep,3)=phi_R ;  %����в�
    if   mod(istep,Res_Freq) ==  0
        fprintf('\n ��������  %s / %s\n  phi_R=   %s\n  u_R  =   %s\n  v_R  =   %s  \n',...
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
        fprintf('���ټ��һ�°� >o< \n');
        break;
    end
    %% Output and Draw
%=====������ݵ�.dat======
if  mod(istep,dat_Freq) ==  0
    var1=phi(imin:imax,jmin:jmax); %1.ָ�������������Χ
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
    varName ='phi psi u v p T cs cpsi sigma PhiE ChargeQ PhiM\n';  % 2.ָ�����������ʾ��;ע��ͬһ.m�ļ��к����varName�ᱻǰ��ĸ���
    D2export(istep,xx(:),yy(:),NGX,NGY,varName,var1(:),var2(:),var3(:),var4(:),var5(:),var6(:),var7(:),var8(:),var9(:),var10(:),var11(:),var12(:))%3.�������
end
%=====���ƶ�̬ͼ======
switch Dynamic_Draw
    case 'ON'
        if     mod(istep,20) ==  0
            contourf(Coord_x(4:1+nx),Coord_y(4:1+ny),phi(4:1+nx,4:1+ny)',3) %����� �ȸ���ͼ
%             Vorticity=D2Gradx_Matrix(v)-D2Grady_Matrix(u);
%             contourf(Coord_x(4:1+nx),Coord_y(4:1+ny),Vorticity(4:1+nx,4:1+ny)',20)   %����ͼ
            title('Moving phi');axis equal;
            drawnow
        end
end
%=====�������ݵ�Workspace====
if   mod(istep,mat_Freq) ==  0   %�������ݵ�Workspace
    filename=strcat('D2Data', num2str(istep));
    MatData=char(strcat('../',Matpath,'/',filename));%��·���ͱ�����ȫ��ƴ�ӳ�char���ͣ���save����
    save (MatData);
end
if   istep==100   %�������ݵ�Workspace
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
% D2ResultDeal(ResultDeal); %�������������{'RB10' 'RB1000' 'laplace''MassLeaf'...}
% end



