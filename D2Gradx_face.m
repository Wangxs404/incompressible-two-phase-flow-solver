function [Grx] = D2Gradx_face(F)     %矩阵操作构成面心梯度算子，全局操作方便，外层值并不使用
global dxi imin imax jmin jmax                  
Grx=dxi*(circshift(F,-1,1)-F);
Grx=Grx(imin-1:imax,jmin:jmax);   %维度，（nx+1） * ny
end   

