function [FaceY] = D2FaceMeanY(F) 
%矩阵操作输出变量在4个面的均值矩阵:（nx+1）*ny和nx*(ny+1)
global imin imax jmin jmax
FaceY=0.5*(F + circshift(F,1,2));
FaceY=FaceY(imin:imax , jmin:jmax+1);%只提取 内部+1 矩阵

end
