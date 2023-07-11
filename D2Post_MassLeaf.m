


clc
clear
run D2DataQuantify

load('D2Data_Initial.mat')
VolumnIniPhi=phi;
VolumnIniPhi(VolumnIniPhi<0)=0;
VolumnIniPhi(VolumnIniPhi>0)=1;
VolumnIni=sum(sum(VolumnIniPhi*dx*dy));

plot((Volumn-VolumnIni)/VolumnIni,'ro');
hold on
plot(zeros(length(Volumn),1),'b')
legend('守恒曲线','质量损失曲线');
ylim([-0.5 0.5])
% global imin imax jmin jmax nx ny sigma0 r
% 
% load('D2Data_Multy');%FVM
% 
% load('D2Data_Initial');
% 
% 
fprintf('MassLeaf= %s',strcat( num2str(((Volumn(end)-VolumnIni)/VolumnIni)*100),"%"));
    