function [FaceX,FaceY] = D2Matrix_FaceMeanF(F) 
%矩阵操作输出变量在4个面的均值矩阵:（nx+1）*ny和nx*(ny+1)
global imin imax jmin jmax


FaceX=0.5*(F + circshift(F,1,1));
FaceY=0.5*(F + circshift(F,1,2));

FaceX=FaceX(imin:imax+1 , jmin:jmax);          %只提取 内部+1 矩阵
FaceY=FaceY(imin:imax , jmin:jmax+1);

end





% F=[1,2,3;4,5,6;7,8,9]'
% circshift(F,1,1)
% circshift(F,1,1)