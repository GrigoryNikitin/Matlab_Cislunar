function dydt = derivatives_Moon_cntrd(t,y, M0_moon)
%DERIVATIVES_MOON_CNTRD calculates acceleration in central gravity field of
%the Moon and disturbing acceleration of the Earth

% y is a state vector in the Moon orbital (initial CS)! reference frame: x,y,z,vx,vy,vz

global Mass_Earth;
global Mass_Moon;
global G;

r = [y(1); y(2); y(3)];
v = [y(4); y(5); y(6)];

[rM, ~] = Moon_crds_from_Earth(t,M0_moon);
rE = -rM; %vector from Moon to Earth
rE_SC = r - rE; %vector from the Earth to the spacecraft

q = r(1)/norm(rE) * (r(1)/norm(rE) - 2*rE(1)/norm(rE)) + r(2)/norm(rE) * (r(2)/norm(rE) - 2*rE(2)/norm(rE)) + r(3)/norm(rE) * (r(3)/norm(rE) - 2*rE(3)/norm(rE));
f = (3 + 3*q^3) / ((1+q)^1.5 + 1) * q;

a = -G*Mass_Moon*r/norm(r)^3 - G*Mass_Earth/norm(rE_SC)^3 * (r + f*rE); %acceleration vector of the spacecraft in Moon centered coordinate system

dydt = [v(1); v(2); v(3); a(1); a(2); a(3)];
end

