function [u,v] = D2InitialUVza(U0,u,v)
%…Ë÷√Za's ÀŸ∂»
global imin imax jmin jmax Coord_x Coord_y
for j=jmin: jmax
    for i=imin: imax
        u(i,j)=-U0*pi*(Coord_y(j)-0.5);
        v(i,j)= U0*pi*(Coord_x(i)-0.5);
    end
end
end

