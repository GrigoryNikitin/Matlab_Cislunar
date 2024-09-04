function dydt = derivatives_CR3BP_nondim(t,y)
%DERIVATIVES_CM calculates dydt in the rotating frame by knowing masses and radius vectors of
%Earth and Moon
% y is a state vector in the center of mass rotating reference frame: x,y,z,vx,vy,vz


r = [y(1); y(2); y(3)];
v = [y(4); y(5); y(6)];

mu  = 0.0121505856;
omega = 1;

r1 = sqrt((r(1)+mu)^2+r(2)^2+r(3)^2);
r2 = sqrt((r(1)-1+mu)^2+r(2)^2+r(3)^2);

%velocity vector of the spacecraft in barycentric rotating coordinate system
v(1) = v(1);% - r(2)*omega;
v(2) = v(2);% + r(1)*omega;
v(3) = v(3);

%acceleration vector of the spacecraft in barycentric rotating coordinate system
a(1) = 2*v(2)+r(1)-(1-mu)/r1^3*(r(1)+mu)-mu/r2^3*(r(1)-1+mu);
a(2) = -2*v(1)+r(2)-(1-mu)/r1^3*(r(2))-mu/r2^3*(r(2));
a(3) = -(1-mu)/r1^3*(r(3))-mu/r2^3*(r(3));

dydt = [v(1); v(2); v(3); a(1); a(2); a(3)];
end

