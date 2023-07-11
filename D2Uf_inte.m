    function[uf,vf,ufL,vfL]=D2Uf_inte(u,v,uL,vL)
    global nx ny imin imax jmin jmax 
        [u,v]=D2set_BC(u,v); [uL,vL]=D2set_BC(uL,vL);
    [X,Y] = meshgrid(1:ny+6,1:nx+6);% 样本网格点
    Xuf=( imin-0.5 :  imax+0.5  )'  ;   Yuf=(jmin  :  jmax); %  Uf查询点
    Xvf=( imin :  imax )'          ;   Yvf=(jmin-0.5  :  jmax+0.5); %  Vf查询点
%     uf = interp2(X,Y,u,Yuf,Xuf,'spline');  ufL = interp2(X,Y,uL,Yuf,Xuf,'spline');
%     vf = interp2(X,Y,v,Yvf,Xvf,'spline');  vfL = interp2(X,Y,vL,Yvf,Xvf,'spline');
    uf = interp2(X,Y,u,Yuf,Xuf);  ufL = interp2(X,Y,uL,Yuf,Xuf);
    vf = interp2(X,Y,v,Yvf,Xvf);  vfL = interp2(X,Y,vL,Yvf,Xvf);
    end