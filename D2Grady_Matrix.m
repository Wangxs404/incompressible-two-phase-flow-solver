function [Gry] = D2Grady_Matrix(Matrix) %体心梯度
global dyi                          
Gry=0.5*dyi*(circshift(Matrix,-1,2)-circshift(Matrix,1,2));

end


