function [r1, r2, inert_sv] = rot2inert(t, mu, rot_sv)
%ROT2INERT transfers barycentric rotating state vector to the inertial frame
%everything is nondimensionalized
%t is nondim too

omega = 1;

rot_angle = -omega*t;
rot = [cos(rot_angle), sin(rot_angle), 0; -sin(rot_angle), cos(rot_angle), 0; 0, 0, 1];

r1 = rot*[-mu,0,0]';
r2 = rot*[1-mu,0,0]';
inert_r = rot*[rot_sv(1); rot_sv(2); rot_sv(3)]; %rotated r of a spacecraft;
inert_v = rot*[rot_sv(4); rot_sv(5); rot_sv(6)]; %rotated v of a spacecraft;

inert_v = inert_v + cross([0; 0; omega], inert_r); %this is right checked

inert_sv = [inert_r', inert_v'];
end

