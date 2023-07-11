function [FaceE,FaceW,FaceN,FaceS] = D2Matrix_FaceMeanC(Matrix)  
%矩阵操作输出变量在4个面的均值矩阵：nx*ny
global imin imax jmin jmax
FaceE=0.5*(Matrix + circshift(Matrix,-1,1)); %注意所有矩阵旋转90°故对x均，须上下移位
FaceW=0.5*(Matrix + circshift(Matrix,1,1));  %取全维度统一计算
FaceN=0.5*(Matrix + circshift(Matrix,-1,2));
FaceS=0.5*(Matrix + circshift(Matrix,1,2));

FaceE=FaceE(imin:imax , jmin:jmax);          %只提取内部矩阵
FaceW=FaceW(imin:imax , jmin:jmax);
FaceN=FaceN(imin:imax , jmin:jmax);
FaceS=FaceS(imin:imax , jmin:jmax);

end                                                         
 
