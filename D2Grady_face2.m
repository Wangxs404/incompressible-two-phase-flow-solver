function [Gry] =D2Grady_face2(F)     %矩阵操作构成面心梯度算子，全局操作方便，外层值并不使用
global dyi                          
Gry=dyi*(circshift(F,-1,2)-F); %注意所有矩阵旋转90°
Gry=Gry(:,1:end-1);   %维度，（nx+1） * ny
end   