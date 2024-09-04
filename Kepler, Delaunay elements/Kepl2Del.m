function Y = Kepl2Del(mu, X, toKepl)
%FROMKEPL2DEL Converts classical orbital elements to Delaunay and vice
%versa

% X is either [l,g,h,L,G,H] or [a,e,i,omega,RAAN,MA]

if toKepl == true
  
    % to Classical
  
    l = X(1); g = X(2); h = X(3); L = X(4); G = X(5); H = X(6);
    
    % SMA
    Y(1,1) = L^2/mu;
    
    % Ecc
    Y(2,1) = sqrt(1 - (G/L)^2);
    
    % Inc
    Y(3,1) = acos(H/G);
    
    % AOP
    omega = g;
    if omega < 0
        omega = omega + 2*pi;
    end
    if omega > 2*pi
        omega = mod(omega,2*pi);
    end
    
    Y(4,1) = omega;
    
    % RAAN
    RAAN = h;
    if RAAN < 0
        RAAN = RAAN + 2*pi;
    end
    if RAAN > 2*pi
        RAAN = mod(RAAN,2*pi);
    end
    Y(5,1) = RAAN;
    
    % MA
    MA = l;
    if MA < 0
        MA = MA + 2*pi;
    end
    if MA > 2*pi
        MA = mod(MA,2*pi);
    end
    Y(6,1) = MA;
    
else

    % to Delaunay
    
    a = X(1); e = X(2); i = X(3); omega = X(4); RAAN = X(5); MA = X(6);
  
    % L
    L = sqrt(mu*a);
    
    % G
    G = L*sqrt(1 - e^2);
    
    % H
    H = G*cos(i);
    
    % l
    l = MA;
    % g
    g = omega;
    % h
    h = RAAN;
    
    Y = [l,g,h,L,G,H]';
    
end
end

