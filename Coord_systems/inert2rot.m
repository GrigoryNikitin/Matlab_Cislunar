function [rot_r1, rot_r2, rot_sv] = inert2rot(r1, r2, v2, inert_sv)
%INERT2ROT transfers inertial barycentric state vectors of the Moon, Earth
%and a spacecraft to the rotating frame

rot_angle = atan2(r2(2), r2(1));
rot = [cos(rot_angle), sin(rot_angle), 0; -sin(rot_angle), cos(rot_angle), 0; 0, 0, 1];

rot_r1 = rot*r1';
rot_r2 = rot*r2';
rot_r = rot*[inert_sv(1); inert_sv(2); inert_sv(3)]; %rotated r of a spacecraft;
rot_v = rot*[inert_sv(4); inert_sv(5); inert_sv(6)]; %rotated v of a spacecraft;

omega = norm(v2) / norm(r2); % omega of Earth-Moon system
% also maybe omega as omega = sqrt(G*(Mass_Earth+Mass_Moon)/R^3), where R
% is the distance between Earth and Moon
rot_v = rot_v - cross([0; 0; omega], rot_r); %this is right checked

rot_sv = [rot_r', rot_v'];
end

