function [u,v] = D2set_BC(u,v)
%ÉèÖÃ±ß½çÌõ¼þ
global imax imin jmax jmin velocity_BCtype
global u_bottom v_bottom u_top  v_top  u_left  v_left  u_right v_right
%1:NS slip EW Solid_NoSlip; 2:Solid Bc ;
%3:NS slip EW Open_NoSlip ; 4:NS slip EW Period
%%
v_top=0;
v_bottom=0;
u_left=0;   v_left=0;
u_right=0;  v_right=0;

%%
switch velocity_BCtype
    case 1
        % NS
        % U £º»¬ÒÆ ¹Ì±Ú£¨k,Const£©=(-1,2*u_bc)
        % V £º»¬ÒÆ ¹Ì±Ú (-1,2*v_bc)
        u(:,jmin-1)=2*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-1)=2*(v_bottom-v(:,jmin))+v(:,jmin);
        u(:,jmin-2)=4*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-2)=4*(v_bottom-v(:,jmin))+v(:,jmin);
        u(:,jmin-3)=6*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-3)=6*(v_bottom-v(:,jmin))+v(:,jmin);

        u(:,jmax+1)=2*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+1)=2*(v_top-v(:,jmax))+v(:,jmax);
        u(:,jmax+2)=4*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+2)=4*(v_top-v(:,jmax))+v(:,jmax);
        u(:,jmax+3)=6*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+3)=6*(v_top-v(:,jmax))+v(:,jmax);
        % EW
        % U £º»¬ÒÆ   ¹Ì±Ú  (-1,2*u_bc)
        % V £ºÎÞ»¬ÒÆ ¹Ì±Ú  (1,0)
        u(imin-1,:)=2*(u_left-u(imin,:))+u(imin,:);
        u(imin-2,:)=4*(u_left-u(imin,:))+u(imin,:);
        u(imin-3,:)=6*(u_left-u(imin,:))+u(imin,:);
        v(imin-1,:)=v(imin,:);
        v(imin-2,:)=v(imin,:);
        v(imin-3,:)=v(imin,:);

        u(imax+1,:)=2*(u_right-u(imax,:))+u(imax,:);
        u(imax+2,:)=4*(u_right-u(imax,:))+u(imax,:);
        u(imax+3,:)=6*(u_right-u(imax,:))+u(imax,:);
        v(imax+1,:)=v(imax,:);
        v(imax+2,:)=v(imax,:);
        v(imax+3,:)=v(imax,:);

    case 2
        % NS
        % U £º»¬ÒÆ ¹Ì±Ú (-1,2*u_bc)
        % V £º»¬ÒÆ ¹Ì±Ú (-1,2*v_bc)
        u(:,jmin-1)=2*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-1)=2*(v_bottom-v(:,jmin))+v(:,jmin);
        u(:,jmin-2)=4*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-2)=4*(v_bottom-v(:,jmin))+v(:,jmin);
        u(:,jmin-3)=6*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-3)=6*(v_bottom-v(:,jmin))+v(:,jmin);

        u(:,jmax+1)=2*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+1)=2*(v_top-v(:,jmax))+v(:,jmax);
        u(:,jmax+2)=4*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+2)=4*(v_top-v(:,jmax))+v(:,jmax);
        u(:,jmax+3)=6*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+3)=6*(v_top-v(:,jmax))+v(:,jmax);
        % EW
        % U £º»¬ÒÆ ¹Ì±Ú (-1,2*u_bc)
        % V £º»¬ÒÆ ¹Ì±Ú (-1,2*v_bc)
        u(imin-1,:)=2*(u_left-u(imin,:))+u(imin,:);
        u(imin-2,:)=4*(u_left-u(imin,:))+u(imin,:);
        u(imin-3,:)=6*(u_left-u(imin,:))+u(imin,:);
        v(imin-1,:)=2*(v_left-v(imin,:))+v(imin,:);
        v(imin-2,:)=4*(v_left-v(imin,:))+v(imin,:);
        v(imin-3,:)=6*(v_left-v(imin,:))+v(imin,:);

        u(imax+1,:)=2*(u_right-u(imax,:))+u(imax,:);
        u(imax+2,:)=4*(u_right-u(imax,:))+u(imax,:);
        u(imax+3,:)=6*(u_right-u(imax,:))+u(imax,:);
        v(imax+1,:)=2*(v_right-v(imax,:))+v(imax,:);
        v(imax+2,:)=4*(v_right-v(imax,:))+v(imax,:);
        v(imax+3,:)=6*(v_right-v(imax,:))+v(imax,:);

    case 3
        % NS
        % U £º»¬ÒÆ ¹Ì±Ú (-1,2*u_bc)
        % V £º»¬ÒÆ ¹Ì±Ú (-1,2*v_bc)
        u(:,jmin-1)=2*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-1)=2*(v_bottom-v(:,jmin))+v(:,jmin);
        u(:,jmin-2)=4*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-2)=4*(v_bottom-v(:,jmin))+v(:,jmin);
        u(:,jmin-3)=6*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-3)=6*(v_bottom-v(:,jmin))+v(:,jmin);

        u(:,jmax+1)=2*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+1)=2*(v_top-v(:,jmax))+v(:,jmax);
        u(:,jmax+2)=4*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+2)=4*(v_top-v(:,jmax))+v(:,jmax);
        u(:,jmax+3)=6*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+3)=6*(v_top-v(:,jmax))+v(:,jmax);
        % EW
        % U £º»¬ÒÆ ¿ª·Å±ß½ç (1,0)
        % V £º»¬ÒÆ ¿ª·Å±ß½ç (1,0)
        u(imin-1,:)=u(imin,:);
        u(imin-2,:)=u(imin,:);
        u(imin-3,:)=u(imin,:);
        v(imin-1,:)=v(imin,:);
        v(imin-2,:)=v(imin,:);
        v(imin-3,:)=v(imin,:);

        u(imax+1,:)=u(imax,:);
        u(imax+2,:)=u(imax,:);
        u(imax+3,:)=u(imax,:);
        v(imax+1,:)=v(imax,:);
        v(imax+2,:)=v(imax,:);
        v(imax+3,:)=v(imax,:);

    case 4
        % NS ¹Ì±Ú
        u(:,jmin-1)=2*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-1)=2*(v_bottom-v(:,jmin))+v(:,jmin);
        u(:,jmin-2)=4*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-2)=4*(v_bottom-v(:,jmin))+v(:,jmin);
        u(:,jmin-3)=6*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-3)=6*(v_bottom-v(:,jmin))+v(:,jmin);

        u(:,jmax+1)=2*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+1)=2*(v_top-v(:,jmax))+v(:,jmax);
        u(:,jmax+2)=4*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+2)=4*(v_top-v(:,jmax))+v(:,jmax);
        u(:,jmax+3)=6*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+3)=6*(v_top-v(:,jmax))+v(:,jmax);
        %EW Period
        u(imax+1,:)=u(imin,:);
        u(imax+2,:)=u(imin+1,:);
        u(imax+3,:)=u(imin+2,:);
        u(imin-1,:)=u(imax,:);
        u(imin-2,:)=u(imax-1,:);
        u(imin-3,:)=u(imax-2,:);
        v(imax+1,:)=v(imin,:);
        v(imax+2,:)=v(imin+1,:);
        v(imax+3,:)=v(imin+2,:);
        v(imin-1,:)=v(imax,:);
        v(imin-2,:)=v(imax-1,:);
        v(imin-3,:)=v(imax-2,:);
    case 5   %È«ÖÜÆÚ±ß½ç
        % NS Period
        u(jmax+1,:)=u(jmin,:);
        u(jmax+2,:)=u(jmin+1,:);
        u(jmax+3,:)=u(jmin+2,:);
        u(jmin-1,:)=u(jmax,:);
        u(jmin-2,:)=u(jmax-1,:);
        u(jmin-3,:)=u(jmax-2,:);
        v(jmax+1,:)=v(jmin,:);
        v(jmax+2,:)=v(jmin+1,:);
        v(jmax+3,:)=v(jmin+2,:);
        v(jmin-1,:)=v(jmax,:);
        v(jmin-2,:)=v(jmax-1,:);
        v(jmin-3,:)=v(jmax-2,:);
        %EW Period
        u(imax+1,:)=u(imin,:);
        u(imax+2,:)=u(imin+1,:);
        u(imax+3,:)=u(imin+2,:);
        u(imin-1,:)=u(imax,:);
        u(imin-2,:)=u(imax-1,:);
        u(imin-3,:)=u(imax-2,:);
        v(imax+1,:)=v(imin,:);
        v(imax+2,:)=v(imin+1,:);
        v(imax+3,:)=v(imin+2,:);
        v(imin-1,:)=v(imax,:);
        v(imin-2,:)=v(imax-1,:);
        v(imin-3,:)=v(imax-2,:);
    case 6   %È«ÖÜÆÚ±ß½ç
        % NS
        % U £º»¬ÒÆ ¹Ì±Ú (-1,2*u_bc)
        % V £º»¬ÒÆ ¹Ì±Ú (-1,2*v_bc)
        u(:,jmin-1)=2*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-1)=2*(v_bottom-v(:,jmin))+v(:,jmin);
        u(:,jmin-2)=4*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-2)=4*(v_bottom-v(:,jmin))+v(:,jmin);
        u(:,jmin-3)=6*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-3)=6*(v_bottom-v(:,jmin))+v(:,jmin);

        u(:,jmax+1)=2*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+1)=2*(v_top-v(:,jmax))+v(:,jmax);
        u(:,jmax+2)=4*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+2)=4*(v_top-v(:,jmax))+v(:,jmax);
        u(:,jmax+3)=6*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+3)=6*(v_top-v(:,jmax))+v(:,jmax);
        %EW Period
        u(imax+1,:)=u(imin,:);
        u(imax+2,:)=u(imin+1,:);
        u(imax+3,:)=u(imin+2,:);
        u(imin-1,:)=u(imax,:);
        u(imin-2,:)=u(imax-1,:);
        u(imin-3,:)=u(imax-2,:);
        v(imax+1,:)=v(imin,:);
        v(imax+2,:)=v(imin+1,:);
        v(imax+3,:)=v(imin+2,:);
        v(imin-1,:)=v(imax,:);
        v(imin-2,:)=v(imax-1,:);
        v(imin-3,:)=v(imax-2,:);
end


end

