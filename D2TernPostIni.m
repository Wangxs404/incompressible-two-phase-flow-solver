%此子程序功能：绘制并导出初始相界面数据
subplot(221)
contourf(Coord_x(4:1+nx),Coord_y(4:1+ny),phiDisp(4:1+nx,4:1+ny)',5);axis equal;
subplot(222)
contourf(Coord_x(4:1+nx),Coord_y(4:1+ny),phi.A(4:1+nx,4:1+ny)',3);axis equal;
subplot(223)
contourf(Coord_x(4:1+nx),Coord_y(4:1+ny),phi.B(4:1+nx,4:1+ny)',3);axis equal;
subplot(224)
contourf(Coord_x(4:1+nx),Coord_y(4:1+ny),phi.C(4:1+nx,4:1+ny)',3);axis equal;
title('Initial phi ');
grid minor;
hold on;
pause(0.001);

%===================输出到Tecplot
x=linspace (dx/2 ,Lx-dx/2, nx); y=linspace (dy/2 ,Ly-dy/2, ny);%将网格中节点的数据输出到Tecplot
[X1,Y1] = meshgrid(x,y);
xx = X1'                ;  yy = Y1';
NGX = size(x,2)         ;  NGY = size(y,2);
var.I1=phiDisp(imin:imax,jmin:jmax);  %1.指定输出变量及范围
var.I2=phi.A(imin:imax,jmin:jmax);var.I3=phi.B(imin:imax,jmin:jmax);var.I4=phi.C(imin:imax,jmin:jmax);
var.I5=u(imin:imax,jmin:jmax);
var.I6=v(imin:imax,jmin:jmax);
varName ='phiDisp_Initial phiA_Initial phiB_Initial phiC_Initial\n';    % 2.指定输出变量显示名;注意同一.m文件中后面的varName会被前面的覆盖
D2export('InitialPhi',xx(:),yy(:),NGX,NGY,varName,var.I1(:),var.I2(:),var.I3(:),var.I4(:),var.I5(:),var.I6(:))%3.输出变量

     
        
%================保存初始数据到工作区        
save D2Data_Initial