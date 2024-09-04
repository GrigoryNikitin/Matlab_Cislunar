function [r, v] = Moon_crds_from_Earth(t,M0)
%MOON_CRDS_FROM_EARTH calculates Moon coordinates from Earth centered
%orbital coordinate system

ecc = 0.0549006; % mean eccentricity 0.0549006 %ZERO FOR NOW
a = 384748;
mu_earth = 398600.44;
M = M0 + sqrt(mu_earth/a^3)*t;
EA = Kepler_Eqn_solver(M, ecc, 10^-10);
f = atan2((1-ecc^2)*sin(EA)/(1-ecc*cos(EA)),(cos(EA)-ecc)/(1-ecc*cos(EA)));
if f<0
    f = f + 2*pi;
end
r_M = a*(1-ecc^2)/(1+ecc*cos(f));

%r = [r*cos(f); r*sin(f); 0];
r = [a*(cos(EA) - ecc); a*sqrt(1-ecc^2)*sin(EA); 0];
v = [ -sqrt(mu_earth*a)/r_M*sin(EA); sqrt(mu_earth*a*(1-ecc^2))/r_M*cos(EA); 0];
end

