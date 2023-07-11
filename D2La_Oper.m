function [La_F]=D2La_Oper(F)
global dxi dyi
La_F=-2*(dxi^2+dyi^2)* F ...
               +dxi^2*circshift(F,-1,1)...
               +dxi^2*circshift(F,1,1)...
               +dyi^2*circshift(F,-1,2)...
               +dyi^2*circshift(F,1,2);    
end