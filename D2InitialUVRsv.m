function [u,v] = D2InitialUVRsv(u,v,time,T)
global imin imax jmin jmax Coord_x Coord_y
for j=jmin: jmax
    for i=imin: imax
        u(i,j)=-2*cos(pi*time/T)*sin(pi*Coord_x(i))^2*sin(pi*Coord_y(j))*cos(pi*Coord_y(j));
        v(i,j)= 2*cos(pi*time/T)*sin(pi*Coord_y(j))^2*sin(pi*Coord_x(i))*cos(pi*Coord_x(i));
    end
end
end

