function [rot_sv] = Orbital2MoonEq_inert(inert_sv, toEquatorial)
%Inertial2MoonEq_inert transfers inertial Moon centered state vector (in
%the Moon's orbital plane)
%of a spacecraft to the Moon equatorial non-rotating inertial
%frame and back

% Equatorial inertial frame is opposite of initial inertial frame (OX
% points at the Earth from the Moon's pericenter point) and in equatorial
% plane of the Moon

rot_r = [inert_sv(1); inert_sv(2); inert_sv(3)];
rot_v = [inert_sv(4); inert_sv(5); inert_sv(6)];

if toEquatorial == true %to the Equatorial plane
    
    %the first rotation for X axis to point in the opposite direction
    rot_angle1 = + pi;
    rot_z = [cos(rot_angle1), sin(rot_angle1), 0; -sin(rot_angle1), cos(rot_angle1), 0; 0, 0, 1];
    
    %the second one is to rotate into the equatorial plane
    rot_angle2 = deg2rad(-6.68);
    rot_y = [cos(rot_angle2), 0, -sin(rot_angle2); 0, 1, 0; sin(rot_angle2), 0, cos(rot_angle2)];
    
    %first rotation around Z
    rot_r = rot_z*rot_r; %rotated r of a spacecraft;
    rot_v = rot_z*rot_v; %rotated v of a spacecraft;
    %second rotation around Y
    rot_r = rot_y*rot_r; %rotated r of a spacecraft;
    rot_v = rot_y*rot_v; %rotated v of a spacecraft;
else
    
    rot_angle1 = - pi;
    rot_z = [cos(rot_angle1), sin(rot_angle1), 0; -sin(rot_angle1), cos(rot_angle1), 0; 0, 0, 1];
    
    rot_angle2 = -deg2rad(-6.68);
    rot_y = [cos(rot_angle2), 0, -sin(rot_angle2); 0, 1, 0; sin(rot_angle2), 0, cos(rot_angle2)];
    
    %first rotation around Y
    rot_r = rot_y*rot_r; %rotated r of a spacecraft;
    rot_v = rot_y*rot_v; %rotated v of a spacecraft;
    %third rotation around Z
    rot_r = rot_z*rot_r; %rotated r of a spacecraft;
    rot_v = rot_z*rot_v; %rotated v of a spacecraft;
    
end

rot_sv = [rot_r', rot_v'];
end

