function [FLuxMy,Flux_n,Flux_s]=D2Matrix_FluxY(F)  %构造界面通量矩阵▽y·Φ（nx*ny）
global weno5
[Fs2,Fs1,Fpy,Fn1,Fn2,Fn3]= D2weno_eleY(F);         %1、取出元素矩阵,元素维度：(nx)*(ny+1)

D2weno_fluxY=0.5*(weno5(Fs2,Fs1,Fpy,Fn1,Fn2)+weno5(Fn3,Fn2,Fn1,Fpy,Fs1 ));        %2、调用WENO矩阵函数，生成界面值矩阵 维度：(nx)*(ny+1)


Flux_n=D2weno_fluxY(:,2:end);                   %3、不需要移位，删除多余列即可，维度：nx*ny
Flux_s=D2weno_fluxY(:,1:end-1);

FLuxMy=Flux_n - Flux_s ;

end