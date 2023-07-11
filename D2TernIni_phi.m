function [phi,phiL,phiDisp]=D2TernIni_phi(initial_type)
global epsilon imax imin jmax jmin Coord_x Coord_y Ly Lx nx ny Dim
%1.Bubble_Rising;2.Dam_Break;3.Zalesak's Disk;4.Rayleigh Taylor;
%5.Double drop;6.Drop breakup;7.Ellipse Drop1
% phi=zeros(imax+3,jmax+3);  %和u/v初始化为同样维度，方便定位,[-1 为重]
phi.A=zeros(imax+3,jmax+3) ; phi.B=zeros(imax+3,jmax+3) ;  phi.C=zeros(imax+3,jmax+3); 


for j=jmin: jmax
    for i=imin: imax
        switch initial_type % Light is "1",Heavy is "-1"
            case 1        
                %  Phi.A    
                X0=0.5*Dim;Y0=0.5*Dim;r=0.25*Dim;
                phi.A(i,j)  = 0.5-0.5*tanh((sqrt((Coord_x(i)-X0)^2 + (Coord_y(j)-Y0)^2)-r)/(sqrt(2)*epsilon));%分母是否需要*2                phi(i,j)  = 0.5+0.5*tanh((sqrt((Coord_x(i)-0.5)^2 + (Coord_y(j)-0.5)^2)-0.15)/(2*epsilon));
                %  Phi.B    
                X0=1.5*Dim;Y0=0.5*Dim;r=0.25*Dim;
                phi.B(i,j)  = 0.5-0.5*tanh((sqrt((Coord_x(i)-X0)^2 + (Coord_y(j)-Y0)^2)-r)/(sqrt(2)*epsilon));%                 phi(i,j)  = 0.5+0.5*tanh((sqrt((Coord_x(i)-0.5)^2 + (Coord_y(j)-0.5)^2)-0.15)/(2*epsilon));
                %  Phi.C
                phi.C= 1-phi.A-phi.B; 

%             case 2
%                 X1=3.2*Dim;Y1=2.5*Dim;r1=0.7*Dim;
%                 X2=4.8*Dim;Y2=1.5*Dim;r2=0.7*Dim;
%                 if Coord_x(i) <= 4
%                     phi(i,j)  = -tanh((sqrt((Coord_x(i)-X1)^2 + (Coord_y(j)-Y1)^2)-r1)/(sqrt(2)*epsilon));%                 phi(i,j)  = 0.5+0.5*tanh((sqrt((Coord_x(i)-0.5)^2 + (Coord_y(j)-0.5)^2)-0.15)/(2*epsilon));
%                 else
%                     phi(i,j)  = -tanh((sqrt((Coord_x(i)-X2)^2 + (Coord_y(j)-Y2)^2)-r2)/(sqrt(2)*epsilon));%                 phi(i,j)  = 0.5+0.5*tanh((sqrt((Coord_x(i)-0.5)^2 + (Coord_y(j)-0.5)^2)-0.15)/(2*epsilon));
%                 end
        end
    end
end
%由于未定义外层坐标，所以给一次边界条件 %可兼容润湿性，相当于初始无润湿，此后立刻润湿
phi.A=D2set_BCNeu(phi.A);   phi.B=D2set_BCNeu(phi.B);  phi.B=D2set_BCNeu(phi.B);  
%  初始计算需给定phiL
phiL.A=phi.A  ;phiL.B=phi.B  ;phiL.C=phi.C  ;
%  通过叠1/2/3的后处理方式，使三相界面易于观察。
[phiDisp]=2*phi.A+3*phi.B+phi.C; phiDisp=D2set_BCNeu(phiDisp);
end

