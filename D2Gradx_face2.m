function [Grx] = D2Gradx_face2(F)     %矩阵操作构成面心梯度算子，全局操作方便，外层值并不使用
global dxi                
Grx=dxi*(circshift(F,-1,1)-F);
Grx=Grx(1:end-1,:);   %维度，（nx+1） * ny
end   

