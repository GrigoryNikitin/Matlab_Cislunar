function [sv0] = sv0_earth_circular_orbit(r0, M0_moon)
%SV0_EARTH_CIRCULAR_ORBIT returns initial state vector of a spacecraft on
% a circular orbit around the Earth in barycentric coordinate system
global Mass_Earth;
global G;
%M0_moon is M0 of Moon on its orbit

v0 = sqrt(G*Mass_Earth/r0);
sv0 = [r0; 0; 0; 0; v0; 0]; %Center of mass reference frame
[r_EM, v_EM] = Moon_crds_from_Earth(0, M0_moon); %Moon radius-vector in Earth centered orbital ref frame
[r1, ~] = calc_barycenter(0, M0_moon);
v_cm = norm(v_EM)*norm(r1)/norm(r_EM); %speed of the center of mass with respect to the Earth
v_cmE = -v_EM/norm(v_EM)*v_cm; %Earth velocity from barycenter
sv0 = [sv0(1) + r1(1); sv0(2) + r1(2); sv0(3) + r1(3); sv0(4) + v_cmE(1); sv0(5) + v_cmE(2); sv0(6) + v_cmE(3)]; %new r0 in barycentric ref frame
end