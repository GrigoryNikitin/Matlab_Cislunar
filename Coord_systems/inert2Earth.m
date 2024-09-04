function [rot_sv] = inert2Earth(r1, v1, inert_sv)
%INERT2Earth transfers inertial barycentric state vectors of a
%spacecraft to the Earth centered rotating frame

%r1 and v1 are r and v of the Moon from the barycenter!

%changing the center of CS to the Earth center
rot_r = [inert_sv(1) - r1(1); inert_sv(2) - r1(2); inert_sv(3) - r1(3)];
rot_v = [inert_sv(4) - v1(1); inert_sv(5) - v1(2); inert_sv(6) - v1(3)];

rot_angle1 = atan2(r1(2), r1(1)) + pi;
rot_z = [cos(rot_angle1), sin(rot_angle1), 0; -sin(rot_angle1), cos(rot_angle1), 0; 0, 0, 1];

rot_angle2 = deg2rad(6.68);
rot_y = [cos(rot_angle2), 0, -sin(rot_angle2); 0, 1, 0; sin(rot_angle2), 0, cos(rot_angle2)];

%first rotation around Z
rot_r = rot_z*rot_r; %rotated r of a spacecraft;
rot_v = rot_z*rot_v; %rotated v of a spacecraft;
%second rotation around Y
rot_r = rot_y*rot_r; %rotated r of a spacecraft;
rot_v = rot_y*rot_v; %rotated v of a spacecraft;

rot_sv = [rot_r', rot_v'];
end