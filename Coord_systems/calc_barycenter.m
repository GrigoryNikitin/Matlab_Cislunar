function [rE, rM] = calc_barycenter(t,M0_moon)
%CALC_CENTER_MASS returns r of the Earth and the Moon in barycentric
%coordinate system

global Mass_Earth;
global Mass_Moon;

r_moon_from_earth = Moon_crds_from_Earth(t, M0_moon); %Moon radius-vector in Earth centered orbital ref frame

Center_Mass = r_moon_from_earth.*Mass_Moon/(Mass_Earth + Mass_Moon);

rE = -Center_Mass; %Earth vector from CM of the system
rM = r_moon_from_earth - Center_Mass; %Moon vector from CM of the system
end

