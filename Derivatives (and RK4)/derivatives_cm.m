function dydt = derivatives_cm(t,y, M0_moon)
%DERIVATIVES_CM calculates dydt by knowing masses and radius vectors of
%Earth and Moon
% y is a state vector in the center of mass reference frame: x,y,z,vx,vy,vz

global Mass_Earth;
global Mass_Moon;
global G;

r = [y(1); y(2); y(3)];
v = [y(4); y(5); y(6)];

[r1, r2] = calc_barycenter(t, M0_moon);

a = G*Mass_Earth*(r1-r)/norm(r1-r)^3 + G*Mass_Moon*(r2-r)/norm(r2-r)^3; %acceleration vector on the spacecraft in barycentric coordinate system

dydt = [v(1); v(2); v(3); a(1); a(2); a(3)];
end

