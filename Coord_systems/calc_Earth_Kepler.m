function [Kepler_Earth] = calc_Earth_Kepler(t, M0_moon)
% Earth's Kepler elements are calculated below
%the vector is in the Inertial Orbital Plane CS (has to be transformed to MoonEq Inertial)
global runit;
global vunit;
global tunit;
global mu_Earth;
global mu_Moon;

[rM, vM] = Moon_crds_from_Earth(t*tunit, M0_moon);
rE = -rM/runit; %nondimensionalized vector from Moon to Earth
vE = -vM/vunit; %nondimensionalized velocity of Earth with respect to Moon

% STOPPED ROTATING FOR NOW AS I NEED EARTH TO BE IN OXY PLANE!
%rotating rE and vE into the Equatorial plane
rot_z = [cos(pi), sin(pi), 0; -sin(pi), cos(pi), 0; 0, 0, 1];
%rot_y = [cosd(-6.68), 0, -sind(-6.68); 0, 1, 0; sind(-6.68), 0, cosd(-6.68)];
rE = rot_z*rE;
%rE = rot_y*rE; %now rE is in Equatorial inertial
vE = rot_z*vE;
%vE = rot_y*vE; % noew vE is in Equatorial inertial
% STOPPED ROTATING FOR NOW AS I NEED EARTH TO BE IN OXY PLANE!

Kepler_Earth = Dec2Kepl(mu_Earth/mu_Moon, [rE, vE]);
end

