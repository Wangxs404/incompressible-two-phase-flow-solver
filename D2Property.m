global gx gy  M epsilon  rho_Heavy rho_Light  mu_Light  mu_Heavy Dim
global Pi Ex cM sigma0 K Beta BetaS
global EK1 EK2 Eepsilon1 Eepsilon2  PhiE_N  PhiE_S
global miuM1 miuM2 H0  upsilon0 rho0
%改为了sigma0!!
% 1. RisingBubble Case1   1. RisingBubble Case2   3. MassTest

switch property  %常规
    case 1
        mu_Light=1  ;  mu_Heavy=10  ;  rho_Light =100  ;  rho_Heavy=1000; % “-1”represent “the Heavy one”
        gx=0   ;     gy=-0.98;                 % g向上为正
        M=1e-1 ;     epsilon=0.02;      sigma0=24.5;        %迁移率 界面厚度 表面张力系数
    case 2
        mu_Light=0.1  ;  mu_Heavy=10  ;  rho_Light =1  ;  rho_Heavy=1000; % “-1”represent “the Heavy one”
        gx=0   ;     gy=-0.98;                 % g向上为正
        M=1e-2 ;     epsilon=1e-2;      sigma0=1.96;        %迁移率 界面厚度 表面张力系数
    case 3    %用于质量泄漏测试
        mu_Light=1/100  ;  mu_Heavy=1/100   ;  rho_Light =1 ;  rho_Heavy=1; % "-1"represent "the Heavy one"
        gx=0   ;     gy=0;                 % g向上为正
        M=1e-2 ;     epsilon=0.02;       sigma0=1e-12;          %迁移率 界面厚度 表面张力系数
    case 31    %用于磁场测试
        mu_Light=16e-3  ;  mu_Heavy=0.8e-3  ;  rho_Light =1.58e3 ;  rho_Heavy=0.8e3; % “-1”represent “the Heavy one”
        gx=0   ;     gy=0;                 % g向上为正
        M=1e-7 ;     epsilon=0.2e-3;      sigma0=3.07e-3;          %迁移率 界面厚度 表面张力系数

    case 4   %用于静态液滴
        mu_Light=0.1  ;  mu_Heavy=1  ;  rho_Light =1  ;  rho_Heavy=1000; % "-1"represent "the Heavy one"
        gx=0   ;     gy=-9.8;                 % g向上为正
        M=1e-2 ;     epsilon=2e-2;      sigma0=1;        %迁移率 界面厚度 表面张力系数
    case 5 %real
        mu_Light = 1.78e-5 ;  mu_Heavy = 1.002e-3 ;  rho_Light = 1.204 ;  rho_Heavy=998.207; % “-1”represent “the Heavy one”
        gx=0   ;     gy=-9.8;                 % g向上为正
        M=1e-2;     epsilon=0.02*Dim;      sigma0=7.28e-2;        %迁移率 界面厚度 表面张力系数
    case 6 %剪切流液滴变形
        mu_Light=1e-2  ;  mu_Heavy=1e-2  ;  rho_Light =1 ;  rho_Heavy=1; % “-1”represent “the Heavy one”
        gx=0   ;     gy=0;                 % g向上为正
        Ca=0.2;  %毛细数(ratio between viscous and surface tension）
        sigma0=1e-5/Ca; %sigma0=0.01;
        M=1e-1;     epsilon=2e-2;       %迁移率 界面厚度 表面张力系数
   case 7 %双液滴融合
        mu_Light=1  ;  mu_Heavy=10  ;  rho_Light =100  ;  rho_Heavy=1000; % “-1”represent “the Heavy one”
        gx=0   ;     gy=0;                 % g向上为正
        M=1e-5 ;     epsilon=0.02;      sigma0=24.5;        %迁移率 界面厚度 表面张力系数
   case 8 %双液滴融合
        mu_Light=0.1  ;  mu_Heavy=10  ;  rho_Light =1  ;  rho_Heavy=1000; % “-1”represent “the Heavy one”
        gx=0   ;     gy=0;                 % g向上为正
        M=1e-5 ;     epsilon=2e-2;      sigma0=1.96;       %迁移率 界面厚度 表面张力系数
  case 9 %Capillary wave
        mu_Light=0.064720863  ;  mu_Heavy=0.064720863   ;  rho_Light =1  ;  rho_Heavy=1000; % "-1"represent "the Heavy one”
        gx=0   ;     gy=0;                 % g向上为正
        M=1e-1 ;     epsilon=2e-2;      sigma0=2;       %迁移率 界面厚度 表面张力系数
  case 10 %RT
        mu_Light=1/1000  ;  mu_Heavy=1/1000  ;  rho_Light =1  ;  rho_Heavy=3; % "-1"represent "the Heavy one"
        gx=0   ;     gy=-1;                 % g向上为正
        M=1e-2 ;     epsilon=1e-2;      sigma0=0;       %迁移率 界面厚度 表面张力系数
  case 11 %RT2
        mu_Light=1/10  ;  mu_Heavy=1/10  ;  rho_Light =1  ;  rho_Heavy=3; % "-1"represent "the Heavy one"
        gx=0   ;     gy=-0.01;                 % g向上为正
        M=1e-2 ;     epsilon=2e-2;      sigma0=0;       %迁移率 界面厚度 表面张力系数

end

% Surfactant Property=================
Pi=1.35 ; Ex=0.9;  cM=1e-7;   %净表面张力 - 温度常数 - 溶解度 - 迁移率
BetaS=0.9;%弹性数（强度）
K= 3/8*epsilon*sigma0;Beta=0.75*sigma0/epsilon;

% Electronic Property=================
EK1=4e-8   ;  EK2=7e-8;      %导电率   EK(比值决定了变形方向)
Eepsilon1=0.01 ;  Eepsilon2=0.01;   %介电常数 Eepsilon
PhiE_N=10  ;  PhiE_S=0;

% Electronic Property=================
% miuM1=4.02e-6   ;  miuM2=0.5*miuM1;      %磁导率
miuM1=4.02e-6   ;  miuM2=0.5*miuM1;
H0=3.7;

%常系数算法常量
upsilon0=0.5*mu_Heavy/rho_Light; %大于等于
rho0= rho_Light; %小于等于