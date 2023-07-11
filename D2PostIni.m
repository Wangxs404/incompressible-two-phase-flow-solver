%此子程序功能：绘制并导出初始相界面数据
figure(1)
contourf(Coord_x(4:1+nx),Coord_y(4:1+ny),phi(4:1+nx,4:1+ny)',3)
title('Initial phi ');
axis equal;
grid minor;
hold on;
pause(0.001);
%===================输出到Tecplot
x=linspace (dx/2 ,Lx-dx/2, nx); y=linspace (dy/2 ,Ly-dy/2, ny);%将网格中节点的数据输出到Tecplot
[X1,Y1] = meshgrid(x,y);
xx = X1'                ;  yy = Y1';
NGX = size(x,2)         ;  NGY = size(y,2);
var.I1=phi(imin:imax,jmin:jmax);  %1.指定输出变量及范围
var.I2=cs(imin:imax,jmin:jmax);
var.I3=u(imin:imax,jmin:jmax);
var.I4=v(imin:imax,jmin:jmax);
var.I5=T(imin:imax,jmin:jmax);
varName ='phi_Initial cs_Initial u_Initial v_Initial T-Ini\n';    % 2.指定输出变量显示名;注意同一.m文件中后面的varName会被前面的覆盖
D2export('InitialPhi',xx(:),yy(:),NGX,NGY,varName,var.I1(:),var.I2(:),var.I3(:),var.I4(:),var.I5(:))%3.输出变量

     
        
%================保存初始数据到工作区        
save(char(strcat('../',Matpath,'/','D2Data_Initial')))
