function [FLuxMx, Flux_e,Flux_w]=D2Matrix_FluxX(F)  %构造界面通量矩阵▽x·Φ（nx*ny）
global weno5
[Fw2,Fw1,Fpx,Fe1,Fe2,Fe3]= D2weno_eleX(F);  %1、取出元素矩阵,元素维度：(nx+1)*ny

weno_fluxX=0.5*(weno5(Fw2,Fw1,Fpx,Fe1,Fe2) +weno5(Fe3,Fe2,Fe1,Fpx,Fw1 ));
                                          %2、调用WENO矩阵函数，生成界面值矩阵,维度：(nx+1)*ny

Flux_e=weno_fluxX(2:end,:);            %3、不需要移位，删除多余行即可，维度：nx*ny
Flux_w=weno_fluxX(1:end-1,:);

FLuxMx = Flux_e - Flux_w ;

end
