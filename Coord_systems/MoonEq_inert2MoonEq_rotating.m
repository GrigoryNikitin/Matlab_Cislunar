function [rot_sv] = MoonEq_inert2MoonEq_rotating(t, inert_sv, toRotating)
%InertMoonEq2MoonEq_rotating transfers a state vector from Moon equatorial non-rotating inertial
%frame to the Moon equatorial rotating frame (OX points on Earth) and back
% just 1 rotation

global we;

rot_r = [inert_sv(1); inert_sv(2); inert_sv(3)];
rot_v = [inert_sv(4); inert_sv(5); inert_sv(6)];

if toRotating == true %to rotating CS
    
    rot_angle1 = we*t; %nondim w * nondim time
    rot_z = [cos(rot_angle1), sin(rot_angle1), 0; -sin(rot_angle1), cos(rot_angle1), 0; 0, 0, 1];
    
    %rotation around Z
    rot_r = rot_z*rot_r; %rotated r of a spacecraft;
    rot_v = rot_z*rot_v; %rotated v of a spacecraft;

else
    
    rot_angle1 = -we*t; %nondim w * nondim time
    rot_z = [cos(rot_angle1), sin(rot_angle1), 0; -sin(rot_angle1), cos(rot_angle1), 0; 0, 0, 1];
    
    
    %rotation around Z
    rot_r = rot_z*rot_r; %rotated r of a spacecraft;
    rot_v = rot_z*rot_v; %rotated v of a spacecraft;
    
end

rot_sv = [rot_r', rot_v'];
end

