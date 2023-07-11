function [FaceX] = D2FaceMeanX(F) 
%矩阵操作输出变量在4个面的均值矩阵:（nx+1）*ny和nx*(ny+1)
global imin imax jmin jmax

FaceX=0.5*(F + circshift(F,1,1));
FaceX=FaceX(imin:imax+1 , jmin:jmax);          %只提取 内部+1 矩阵

end