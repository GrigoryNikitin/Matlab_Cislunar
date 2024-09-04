function [sv0] = sv0_moon_circular_orbit(r0, M0_moon)
%SV0_MOON_CIRCULAR_ORBIT returns initial state vector of a spacecraft on
% a circular orbit around the Moon in barycentric coordinate system
global Mass_Moon;
global G;
%M0_moon is M0 of Moon on its orbit

v0 = sqrt(G*Mass_Moon/r0);
sv0 = [r0; 0; 0; 0; v0; 0]; %Center of mass reference frame
[r_EM, v_EM] = Moon_crds_from_Earth(0, M0_moon); %Moon radius-vector in Earth centered orbital ref frame
[r1, r2] = calc_barycenter(0, M0_moon);
v_cm = norm(v_EM)*norm(r1)/norm(r_EM); %speed of the center of mass with respect to the Earth
v_cmM = v_EM/norm(v_EM)*(norm(v_EM) - v_cm); %Moon velocity from barycenter
sv0 = [sv0(1) + r2(1); sv0(2) + r2(2); sv0(3) + r2(3); sv0(4) + v_cmM(1); sv0(5) + v_cmM(2); sv0(6) + v_cmM(3)]; %new r0 in barycentric ref frame
end

