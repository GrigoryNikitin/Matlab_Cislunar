function dydt = derivatives_moon(t,y)
%DERIVATIVES_CM calculates dydt by knowing masses and radius vectors of
%Earth and Moon
% y is a state vector in the Moon orbital! reference frame: x,y,z,vx,vy,vz

global Mass_Moon;
global G;

r = [y(1); y(2); y(3)];
v = [y(4); y(5); y(6)];

a = -G*Mass_Moon*r/norm(r)^3; %acceleration vector on the spacecraft in Moon centered coordinate system

dydt = [v(1); v(2); v(3); a(1); a(2); a(3)];
end

