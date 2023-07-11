function [phi,phiL,u_top,u_bottom]=D2Initial_phi(initial_type,InputU_top,InputU_bottom)
global epsilon imax imin jmax jmin Coord_x Coord_y Ly Lx nx ny Dim r
%1.Bubble_Rising;2.Dam_Break;3.Zalesak's Disk;4.Rayleigh Taylor;
%5.Double drop;6.Drop breakup;7.Ellipse Drop1
u_top   = InputU_top;
u_bottom= InputU_bottom;
phi=zeros(imax+3,jmax+3);  %和u/v初始化为同样维度，方便定位,[-1 为重]
rr = 0.4 ;  %圆盘参数
y1 = 0.5 - sqrt(0.46^2 - 0.04^2) ;%圆盘参数
y2 = 0.5 - sqrt(0.4^2 - 0.04^2) ;%圆盘参数

e=epsilon;  %界面宽度
Li=1;     %溃坝宽度
He=1;     %溃坝高度

for j=jmin: jmax
    for i=imin: imax
        switch initial_type % Light is "1",Heavy is "-1"
            case 1          %Bubble_Rising
                %                 phi(i,j)  = -tanh((sqrt((Coord_x(i)-0.5e-2)^2 + (Coord_y(j)-0.5e-2)^2)-0.25e-2)/(sqrt(2)*epsilon));
                X0=0.5*Dim;Y0=0.5*Dim;r=0.25*Dim;
                phi(i,j)  = -tanh((sqrt((Coord_x(i)-X0)^2 + (Coord_y(j)-Y0)^2)-r)/(sqrt(2)*epsilon));%                 phi(i,j)  = 0.5+0.5*tanh((sqrt((Coord_x(i)-0.5)^2 + (Coord_y(j)-0.5)^2)-0.15)/(2*epsilon));

            case 2          %Dam_Break
                if  Coord_y(j) >= He-e    &&    Coord_y(j) <= He+e    && Coord_x(i)>=0   &&  Coord_x(i)<=Li-e  %a
                    phi(i,j)  = tanh((Coord_y(j)-He)/(sqrt(2)*epsilon));
                end

                if  Coord_x(i) >= Li-e && Coord_x(i) <= Li+e          &&Coord_y(j) >= 0    &&    Coord_y(j) <= He-e  %b
                    phi(i,j)  = tanh((Coord_x(i)-Li)/(sqrt(2)*epsilon));
                end

                if       Coord_y(j) <= Ly && Coord_y(j) >= He+e   && Coord_x(i) <= Li+e && Coord_x(i) >= 0 ...
                        ||  Coord_x(i) <=Lx && Coord_x(i) >= Li+e
                    phi(i,j)   = 1;
                end
                if    Coord_y(j) <= He-e && Coord_y(j) >= 0   && Coord_x(i) <= Li-e && Coord_x(i) >=0  %3
                    phi(i,j)   = -1;
                end
                if Coord_y(j) <= He+e && Coord_y(j) >= He-e   && Coord_x(i) <= Li+e && Coord_x(i) >= Li-e %c
                    nd=0;ndx=0;ndy=0;
                    for a=1:ny
                        if Coord_y(a)>= He-e && Coord_y(a)<= He+e
                            nd=nd+1;
                        end
                        if  Coord_y(a)<He-e
                            ndy=ndy+1;
                        end
                    end
                    for a=1:nx
                        if  Coord_x(a)< Li-e
                            ndx=ndx+1;
                        end
                    end
                    phi_vector=phi(imin,ndy+1:ndy+nd);
                    AA=repmat( phi_vector,nd,1);
                    b1=triu(fliplr(rot90(AA)));b2=triu(fliplr(rot90(AA)),1);
                    c=rot90(b1)+flipud(b2);
                    phi(ndx+1:ndx+nd,ndy+1:ndy+nd)=rot90(c);
                end


            case 3          %zalesak
                if((Coord_y(j) <= 0.5 && Coord_y(j) >= y1-epsilon && Coord_x(i)<= 0.46) || (Coord_y(j) <= 0.5 && Coord_y(j) >= y1-epsilon && Coord_x(i)>= 0.54))

                    phi(i,j) = tanh(min(rr - sqrt((Coord_x(i) - 0.5)^2+(Coord_y(j) - 0.5)^2),...
                        max(0.46 - Coord_x(i),Coord_x(i) - 0.54 )) * 2 / epsilon);
                end

                if(Coord_x(i)<=0.54 && Coord_x(i)>=0.46  && Coord_y(j)>= 0.5 )

                    phi(i,j)=tanh(min(rr -sqrt((Coord_x(i)-0.5)^2+(Coord_y(j)-0.5)^2),...
                        Coord_y(j) - 0.5) * 2 / epsilon);
                end

                if(Coord_x(i) <= 0.54  && Coord_x(i) >= 0.46  && Coord_y(j)<= 0.5  && Coord_y(j)>= y1 )

                    phi(i,j) = tanh(max(max(rr - Coord_x(i),Coord_x(i) - 0.54),Coord_y(j) - 0.5) * 2 / epsilon);
                end

                if((Coord_x(i) <= 0.46  && Coord_y(j) >= 0.5 ) || (Coord_x(i) >= 0.54  && Coord_y(j) >= 0.5 ))
                    phi(i,j)=tanh(min(rr -sqrt((Coord_x(i) - 0.5)^2+(Coord_y(j) - 0.5)^2),...
                        min(sqrt((Coord_x(i) - 0.46)^2 + (Coord_y(j) - 0.5)^2),sqrt((Coord_x(i) - 0.54)^2 ...
                        +(Coord_y(j) - 0.5)^2))) * 2 / epsilon);
                end

                if(Coord_y(j) <= y1 - epsilon  )

                    phi(i,j) = tanh((rr - sqrt((Coord_x(i) - 0.5)^2+(Coord_y(j) - 0.5)^2)) * 2 / epsilon);
                end

                if(Coord_y(j) <= y2 - epsilon && Coord_x(i) >= 0.5 )  && Coord_x(i) <= 0.54 && Coord_x(i) >= 0.46  && abs(rr - Coord_x(i)) < 0.04*(rr - Coord_y(j))/sqrt((0.46^2 - 0.04^2))
                    phi(i,j) = tanh(min(sqrt((Coord_x(i) - 0.46)^2+(Coord_y(j) - y2)^2),...
                        sqrt((Coord_x(i) - 0.54)^2+(Coord_y(j) - y2)^2)) * 2 / epsilon);
                end
                % 修剪多余区块，1、抹除底部块；2、镜像翻转使对称
                if Coord_y(j) <= y1
                    phi(i,j)=-1;
                end
                phileft=phi; phiright=flipud(phi); phi=[phileft(1:(imax+3)/2,:) ; phiright((imax+3)/2+1:end,:)];

            case 4     %Rayleigh Taylor
                Coord_x=1/Dim*Coord_x;Coord_y=1/Dim*Coord_y;
                Ly=Ly/Dim;
                eta2= - 0.1 * cos(2 * pi * Coord_x(i));
                phi(i,j) = -tanh((Coord_y(j)-Ly/2 + eta2) / (sqrt(2)*e));
                Coord_x=Dim*Coord_x;Coord_y=Dim*Coord_y;
                Ly=Ly*Dim;
            case 5    % 5.Double drop % 6.Drop breakup

                if 0.5*Coord_x(i)+ Coord_y(j) <= 1
                    phi(i,j)  = -tanh((sqrt((Coord_x(i)-0.84)^2 + (Coord_y(j)-0.35)^2)-0.2)/(sqrt(2)*epsilon));

                else
                    phi(i,j)  = -tanh((sqrt((Coord_x(i)-1.16)^2 + (Coord_y(j)-0.65)^2)-0.2)/(sqrt(2)*epsilon));
                end
            case 6    % 6.Drop breakup
                phi(i,j)  = -tanh((sqrt((Coord_x(i)-1)^2 + (Coord_y(j)-1)^2)-0.5)/(sqrt(2)*epsilon));
            case 7    % 7.ellipse Drop
                x0=2  ;   y0=1  ;   a=1;   b=1;  r=0.5;
                phi(i,j)  = -tanh((sqrt(((Coord_x(i)-x0)/a)^2 + ((Coord_y(j)-y0)/b)^2)-r)/(sqrt(2)*epsilon));
            case 8
                phi(i,j) = -tanh((Coord_y(j)-Ly/2 ) / (sqrt(2)*e));
            case 9
                X1=0.33*Dim;Y1=0.5*Dim;r1=0.15*Dim;
                X2=0.67*Dim;Y2=0.5*Dim;r2=0.15*Dim;
                if Coord_x(i) <= 0.5
                    phi(i,j)  = -tanh((sqrt((Coord_x(i)-X1)^2 + (Coord_y(j)-Y1)^2)-r1)/(sqrt(2)*epsilon));%                 phi(i,j)  = 0.5+0.5*tanh((sqrt((Coord_x(i)-0.5)^2 + (Coord_y(j)-0.5)^2)-0.15)/(2*epsilon));
                else
                    phi(i,j)  = -tanh((sqrt((Coord_x(i)-X2)^2 + (Coord_y(j)-Y2)^2)-r2)/(sqrt(2)*epsilon));%                 phi(i,j)  = 0.5+0.5*tanh((sqrt((Coord_x(i)-0.5)^2 + (Coord_y(j)-0.5)^2)-0.15)/(2*epsilon));
                end
            case 10
                phi(i,j)  = -tanh((sqrt((Coord_x(i)-0.5)^2 + (Coord_y(j)-0.75)^2)-0.15)/(sqrt(2)*epsilon));%                 phi(i,j)  = 0.5+0.5*tanh((sqrt((Coord_x(i)-0.5)^2 + (Coord_y(j)-0.5)^2)-0.15)/(2*epsilon));
            case 11
                phi(i,j)=-1;
                if sqrt((Coord_x(i)-2)^2 + (Coord_y(j)-2.75)^2) <= 0.5 %||(abs(Coord_x(i)-2) <= 0.06  ||  Coord_y(j) <=2.85)
                    phi(i,j)  =  1;end
                if abs(Coord_x(i)-2) <= 0.06 &&  Coord_y(j) <=2.85
                    phi(i,j)  =  -1;end
            case 12
                X1=3.2*Dim;Y1=2.5*Dim;r1=0.7*Dim;
                X2=4.8*Dim;Y2=1.5*Dim;r2=0.7*Dim;
                if Coord_x(i) <= 4
                    phi(i,j)  = -tanh((sqrt((Coord_x(i)-X1)^2 + (Coord_y(j)-Y1)^2)-r1)/(sqrt(2)*epsilon));%                 phi(i,j)  = 0.5+0.5*tanh((sqrt((Coord_x(i)-0.5)^2 + (Coord_y(j)-0.5)^2)-0.15)/(2*epsilon));
                else
                    phi(i,j)  = -tanh((sqrt((Coord_x(i)-X2)^2 + (Coord_y(j)-Y2)^2)-r2)/(sqrt(2)*epsilon));%                 phi(i,j)  = 0.5+0.5*tanh((sqrt((Coord_x(i)-0.5)^2 + (Coord_y(j)-0.5)^2)-0.15)/(2*epsilon));
                end
            case 13 %Capillary wave
                y= pi + 0.01*2*pi*cos(Coord_x(i)-0.5*Lx/nx+pi);
                phi(i,j)  = -tanh((Coord_y(j)-y)/(sqrt(2)*epsilon));%
            case 14 %horizontal shear layer
                delta1=1/30; 
                if Coord_y(j) <= 0.5
                    phi(i,j)=tanh((Coord_y(j)-0.25)/delta1);
                else
                    phi(i,j)=tanh((0.75-Coord_y(j))/delta1);
                end
        end
    end
end

phi=D2set_BCNeu(phi);   %由于未定义外层坐标，所以给一次边界条件 %可兼容润湿性，相当于初始无润湿，此后立刻润湿
phiL=phi ;
end

