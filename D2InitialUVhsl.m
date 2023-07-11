function [u,v]  =D2InitialUVhsl(Coord_x,Coord_y,imax,imin,jmax,jmin)
global imin imax jmin jmax
u=zeros(imax+3,jmax+3) ; v=zeros(imax+3,jmax+3);
delta1=1/30; delta2=0.05;
for j=jmin: jmax
    for i=imin: imax
        if Coord_y(j) <= 0.5
            u(i,j)=tanh((Coord_y(j)-0.25)/delta1);
        else
            u(i,j)=tanh((0.75-Coord_y(j))/delta1);
        end

        v(i,j)=delta2*sin(2*pi*Coord_x(i));
    end
end

end