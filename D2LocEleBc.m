global nx ny D2Loc D2LocR AX AY

D2Loc.s=[ ones(nx,1) ; repmat(zeros(nx,1),ny-1,1)] ;
D2Loc.n=[ repmat(zeros(nx,1),ny-1,1) ;ones(nx,1) ];
D2Loc.w=repmat([1; zeros(nx-1,1)],ny,1);
D2Loc.e=repmat([zeros(nx-1,1) ; 1],ny,1);

%0,1互换的函数句柄-----------------
%去边界-位置向量
D2rver=@(x)-1*x+1; 
D2LocR.s=D2rver(D2Loc.s);
D2LocR.n=D2rver(D2Loc.n);
D2LocR.w=D2rver(D2Loc.w);
D2LocR.e=D2rver(D2Loc.e);

% 对角元素在系数矩阵A中的坐标 
 AXp=(1:nx*ny);AYp=(1:nx*ny);
 AXw=(2:nx*ny);AYw=(1:nx*ny-1);
 AXe=(1:nx*ny-1); AYe=(2:nx*ny);
 AXs=(nx+1:nx*ny); AYs=(1:nx*ny-nx);
 AXn=(1:nx*ny-nx); AYn=(nx+1:nx*ny);
 AX=[AXs , AXw , AXp , AXe , AXn]';
 AY=[AYs , AYw , AYp , AYe , AYn]';


