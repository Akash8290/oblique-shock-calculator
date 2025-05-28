function [sol] = oblique_shock(M1,v)
g =1.4;                                                                     % g= gamma
outVar = 0;
i = 1;                                                                      %i =1 for weak shocks, 2 fro strong shocks
% Make sure Mach number is greater than unity
if (M1 <= 1)
    sol = inf;
    fprintf('M1 must be greater than 1\n');
    return;
end

% Solve Turn Angles (Weak)
    theta = v;
    if (theta >= 90)
        sol = inf;
        fprintf('Turning angle is too large\n');
        return;
    end
    
    % Check to make sure turning angle is greater than zero
    if (theta <= 0)
        sol = inf;
        fprintf('Turning angle must be greater than zero\n');
        return;
    end
    
    % Solve for shock angle [deg]
    beta = MDB(g,M1,theta,i);    
                                                                             % Shock angle [deg]
                                                                                 % Turn angle [deg]
% Solve for solution variables
Mn1      = M1*sind(beta);                                                    % Upstream normal Mach number []
Mn2      = M2(g,Mn1);                                                       % Downstream normal Mach number []
Mach2    = Mn2/sind(beta-theta);                                                % Downstream Mach number []
P2P1     = 1+2*g/(g+1)*(Mn1^2-1);                                           % Static pressure ratio []
P02P01   = PP0(g,Mn1)/PP0(g,M2(g,Mn1))*P2P1;                                     % Stagnation pressure ratio []
rho2rho1 = RR0(g,M2(g,Mn1))/RR0(g,Mn1)*P02P01;                              % Static density ratio []
T2T1     = TT0(g,M2(g,Mn1))/TT0(g,Mn1);                                       % Static temperature ratio []

% Set solution variables
    sol.theta  = theta;
    sol.beta   = beta;
    sol.M1     = M1;
    sol.Mn1    = Mn1;
    sol.Mn2    = Mn2;
    sol.M2     = Mach2;
    sol.P2P1   = P2P1;
    sol.P02P01 = P02P01;
    sol.r2r1   = rho2rho1;
    sol.T2T1   = T2T1;
    sol.state  = 'Attached';

% functions

function [M2_Out]  = M2(g,M1)
    M2_Out = sqrt((1+0.5*(g-1)*M1*M1)/(g*M1*M1-0.5*(g-1)));
end

function [pp0_Out] = PP0(g,M)
    pp0_Out = (1+(g-1)/2*(M^2))^(-g/(g-1));
end

function [rr0_Out] = RR0(g,M)
    rr0_Out = (1+(g-1)/2*(M^2))^(-1/(g-1));
end

function [tt0_Out] = TT0(g,M)
    tt0_Out = (1+(g-1)/2*(M^2))^(-1);
end

function [mbd_Out] = MBD(g,M1,b)
    mbd_Out = atand((M1^2*sind(2*b)-2/tand(b))/(2+M1^2*(g+cosd(2*b))));
end

function [mdb_Out] = MDB(g,M1,d,i)
    d = d*(pi/180);
    p = -(M1*M1+2)/M1/M1-g*sin(d)*sin(d);
    q = (2*M1*M1+1)/(M1^4)+((g+1)*(g+1)/4+(g-1)/M1/M1)*sin(d)*sin(d);
    r = -cos(d)*cos(d)/(M1^4);
    
    a = (3*q-p*p)/3;
    b = (2*p*p*p-9*p*q+27*r)/27;
    
    test = b*b/4+a*a*a/27;
    
    if (test > 0)
        mdb_Out = -1;
        return;
    else
        if (test == 0)
            x1 = sqrt(-a/3);
            x2 = x1;
            x3 = 2*x1;
            if (b > 0)
                x1 = -1*x1;
                x2 = -1*x2;
                x3 = -1*x3;
            end
        end

        if (test < 0)
            phi = acos(sqrt(-27*b^2/4/a/a/a));
            x1  = 2*sqrt(-a/3)*cos(phi/3);
            x2  = 2*sqrt(-a/3)*cos(phi/3+pi*2/3);
            x3  = 2*sqrt(-a/3)*cos(phi/3+pi*4/3);
            if (b > 0)
                x1 = -1*x1;
                x2 = -1*x2;
                x3 = -1*x3;
            end
        end
        s1 = x1-p/3;
        s2 = x2-p/3;
        s3 = x3-p/3;

        if (s1 < s2 && s1 < s3)
          t1 = s2;
          t2 = s3;
        elseif (s2 < s1 && s2 < s3)
          t1 = s1;
          t2 = s3;
        else
          t1 = s1;
          t2 = s2;
        end
        
        b1 = asin(sqrt(t1));
        b2 = asin(sqrt(t2));
        
        betas = b1;
        betaw = b2;
        if (b2 > b1)
          betas = b2;
          betaw = b1;
        end
        if (i == 1)
            mdb_Out = betaw*(180/pi);
            return;
        end       
        if (i == 2)
            mdb_Out = betas*(180/pi);
            return;
        end
    end
end

end
