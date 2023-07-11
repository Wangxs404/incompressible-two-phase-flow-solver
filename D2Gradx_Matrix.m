
function [Grx] = D2Gradx_Matrix(Matrix)     %矩阵操作构成体心梯度算子，全局操作方便，外层值并不使用
global dxi                          
Grx=0.5*dxi*(circshift(Matrix,-1,1)-circshift(Matrix,1,1)); %注意所有矩阵旋转90°
end                                                         %故对x求偏导，须上下移位
                                                            %二位1表示行变换，一位1、-1表示正负移位
                                                
