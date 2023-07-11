function[Advection] = D2D_Advc(UF,VF,UFL,VFL,F,FL)  %由面心速度生成离散的对流项
global nx ny dxi dyi imin imax jmin jmax
Advection=zeros(nx+6,ny+6);                    %统一维度，便于计算
UF_ready=2*UF-UFL;     VF_ready=2*VF-VFL;   F_ready=2*F-FL;  %Mcc中F=FL phi中两者不等
[~,Fe,Fw]=D2Matrix_FluxX(F_ready);  %输出界面通量矩阵e、w、n、s,（nx*ny），没有选择别的格式，可以不switch
[~,Fn,Fs]=D2Matrix_FluxY(F_ready);  
UFe=UF_ready(2:end,:);UFw=UF_ready(1:end-1,:);VFn=VF_ready(:,2:end);VFs=VF_ready(:,1:end-1);
Advection(imin:imax,jmin:jmax)= dxi*(UFe .* Fe -  UFw.* Fw) + dyi*(VFn .* Fn -  VFs.* Fs);%对于mcc，计算后需要重构重组形成b-adv

end