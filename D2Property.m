global gx gy  M epsilon  rho_Heavy rho_Light  mu_Light  mu_Heavy Dim
global Pi Ex cM sigma0 K Beta BetaS
global EK1 EK2 Eepsilon1 Eepsilon2  PhiE_N  PhiE_S
global miuM1 miuM2 H0  upsilon0 rho0
%��Ϊ��sigma0!!
% 1. RisingBubble Case1   1. RisingBubble Case2   3. MassTest

switch property  %����
    case 1
        mu_Light=1  ;  mu_Heavy=10  ;  rho_Light =100  ;  rho_Heavy=1000; % ��-1��represent ��the Heavy one��
        gx=0   ;     gy=-0.98;                 % g����Ϊ��
        M=1e-1 ;     epsilon=0.02;      sigma0=24.5;        %Ǩ���� ������ ��������ϵ��
    case 2
        mu_Light=0.1  ;  mu_Heavy=10  ;  rho_Light =1  ;  rho_Heavy=1000; % ��-1��represent ��the Heavy one��
        gx=0   ;     gy=-0.98;                 % g����Ϊ��
        M=1e-2 ;     epsilon=1e-2;      sigma0=1.96;        %Ǩ���� ������ ��������ϵ��
    case 3    %��������й©����
        mu_Light=1/100  ;  mu_Heavy=1/100   ;  rho_Light =1 ;  rho_Heavy=1; % "-1"represent "the Heavy one"
        gx=0   ;     gy=0;                 % g����Ϊ��
        M=1e-2 ;     epsilon=0.02;       sigma0=1e-12;          %Ǩ���� ������ ��������ϵ��
    case 31    %���ڴų�����
        mu_Light=16e-3  ;  mu_Heavy=0.8e-3  ;  rho_Light =1.58e3 ;  rho_Heavy=0.8e3; % ��-1��represent ��the Heavy one��
        gx=0   ;     gy=0;                 % g����Ϊ��
        M=1e-7 ;     epsilon=0.2e-3;      sigma0=3.07e-3;          %Ǩ���� ������ ��������ϵ��

    case 4   %���ھ�̬Һ��
        mu_Light=0.1  ;  mu_Heavy=1  ;  rho_Light =1  ;  rho_Heavy=1000; % "-1"represent "the Heavy one"
        gx=0   ;     gy=-9.8;                 % g����Ϊ��
        M=1e-2 ;     epsilon=2e-2;      sigma0=1;        %Ǩ���� ������ ��������ϵ��
    case 5 %real
        mu_Light = 1.78e-5 ;  mu_Heavy = 1.002e-3 ;  rho_Light = 1.204 ;  rho_Heavy=998.207; % ��-1��represent ��the Heavy one��
        gx=0   ;     gy=-9.8;                 % g����Ϊ��
        M=1e-2;     epsilon=0.02*Dim;      sigma0=7.28e-2;        %Ǩ���� ������ ��������ϵ��
    case 6 %������Һ�α���
        mu_Light=1e-2  ;  mu_Heavy=1e-2  ;  rho_Light =1 ;  rho_Heavy=1; % ��-1��represent ��the Heavy one��
        gx=0   ;     gy=0;                 % g����Ϊ��
        Ca=0.2;  %ëϸ��(ratio between viscous and surface tension��
        sigma0=1e-5/Ca; %sigma0=0.01;
        M=1e-1;     epsilon=2e-2;       %Ǩ���� ������ ��������ϵ��
   case 7 %˫Һ���ں�
        mu_Light=1  ;  mu_Heavy=10  ;  rho_Light =100  ;  rho_Heavy=1000; % ��-1��represent ��the Heavy one��
        gx=0   ;     gy=0;                 % g����Ϊ��
        M=1e-5 ;     epsilon=0.02;      sigma0=24.5;        %Ǩ���� ������ ��������ϵ��
   case 8 %˫Һ���ں�
        mu_Light=0.1  ;  mu_Heavy=10  ;  rho_Light =1  ;  rho_Heavy=1000; % ��-1��represent ��the Heavy one��
        gx=0   ;     gy=0;                 % g����Ϊ��
        M=1e-5 ;     epsilon=2e-2;      sigma0=1.96;       %Ǩ���� ������ ��������ϵ��
  case 9 %Capillary wave
        mu_Light=0.064720863  ;  mu_Heavy=0.064720863   ;  rho_Light =1  ;  rho_Heavy=1000; % "-1"represent "the Heavy one��
        gx=0   ;     gy=0;                 % g����Ϊ��
        M=1e-1 ;     epsilon=2e-2;      sigma0=2;       %Ǩ���� ������ ��������ϵ��
  case 10 %RT
        mu_Light=1/1000  ;  mu_Heavy=1/1000  ;  rho_Light =1  ;  rho_Heavy=3; % "-1"represent "the Heavy one"
        gx=0   ;     gy=-1;                 % g����Ϊ��
        M=1e-2 ;     epsilon=1e-2;      sigma0=0;       %Ǩ���� ������ ��������ϵ��
  case 11 %RT2
        mu_Light=1/10  ;  mu_Heavy=1/10  ;  rho_Light =1  ;  rho_Heavy=3; % "-1"represent "the Heavy one"
        gx=0   ;     gy=-0.01;                 % g����Ϊ��
        M=1e-2 ;     epsilon=2e-2;      sigma0=0;       %Ǩ���� ������ ��������ϵ��

end

% Surfactant Property=================
Pi=1.35 ; Ex=0.9;  cM=1e-7;   %���������� - �¶ȳ��� - �ܽ�� - Ǩ����
BetaS=0.9;%��������ǿ�ȣ�
K= 3/8*epsilon*sigma0;Beta=0.75*sigma0/epsilon;

% Electronic Property=================
EK1=4e-8   ;  EK2=7e-8;      %������   EK(��ֵ�����˱��η���)
Eepsilon1=0.01 ;  Eepsilon2=0.01;   %��糣�� Eepsilon
PhiE_N=10  ;  PhiE_S=0;

% Electronic Property=================
% miuM1=4.02e-6   ;  miuM2=0.5*miuM1;      %�ŵ���
miuM1=4.02e-6   ;  miuM2=0.5*miuM1;
H0=3.7;

%��ϵ���㷨����
upsilon0=0.5*mu_Heavy/rho_Light; %���ڵ���
rho0= rho_Light; %С�ڵ���