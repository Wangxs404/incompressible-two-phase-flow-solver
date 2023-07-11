    function [S,alpha,lamda]=D2Phase_const(epsilon,M,sigma,dt)
        eta=epsilon;
        gamma0=1.5;
        gamma1=M;%qianyilv
        
        lamda=(3*sigma*eta) / (2*sqrt(2));
        S=eta^2 * sqrt( (4*gamma0) / ( lamda*gamma1*dt )  );% >=  
        alpha=(-S/(2*eta^2)) * (  1+   sqrt(  (1-(4*gamma0*eta^4)/(lamda*gamma1*dt*S^2))   )        );

    end