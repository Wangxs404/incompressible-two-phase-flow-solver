global gx gy  M epsilon  rhoIni muIni
global Pi Ex cM sigma0 K Beta BetaS
global EK1 EK2 Eepsilon1 Eepsilon2  PhiE_N  PhiE_S
%改为了sigma0!!
% 1. RisingBubble Case1   1. RisingBubble Case2   3. MassTest

switch property  %常规
    case 1
        rhoIni.A =100   ;  rhoIni.B=500  ;   rhoIni.C=1000;
        muIni.A  =1    ;  muIni.B =2    ;   muIni.C=3  ;
        sigma_12=1   ;  sigma_13=0.01  ;   sigma_23=0.001;
        gx=0   ;     gy=-0.05;                 % g向上为正
        M=1e-5 ;     epsilon=2e-2;            %迁移率 界面厚度 表面张力系数
end
% 暂且使用固定M
[sigmaTern,epsilonPie]=D2TernConst(sigma_12,sigma_13,sigma_23,epsilon);%参数准备，用于计算3相化学势

%% 多物理场参数
% Surfactant Property=================
Pi=1.35 ; Ex=0.9;  cM=1e-7;   %净表面张力 - 温度常数 - 溶解度 - 迁移率
BetaS=0.9;%弹性数（强度）
K= 3/8*epsilon*sigma0;Beta=0.75*sigma0/epsilon;

% Electronic Property=================
EK1=4e-8   ;  EK2=7e-8;      %导电率   EK(比值决定了变形方向)
Eepsilon1=0.01 ;  Eepsilon2=0.01;   %介电常数 Eepsilon
PhiE_N=10  ;  PhiE_S=0;