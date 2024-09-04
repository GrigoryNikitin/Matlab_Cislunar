function dydt = derivatives_Moon_cntrd_nondim(t,y, M0_moon, Degree, Order, Clm, Slm)
%DERIVATIVES_MOON_CNTRD calculates acceleration in central gravity field of
%the Moon (+harmonics) and disturbing acceleration of the Earth

% y is a state vector in the Moon EQUATORIAL INERTIAL reference frame: x,y,z,vx,vy,vz

global mu_Earth;
global mu_Moon;
global runit;
global tunit;
global we;

r = [y(1); y(2); y(3)];
v = [y(4); y(5); y(6)];

[~, r2] = calc_barycenter(t*tunit,M0_moon); %positions of the Earth and the Moon from the barycenter
% r2 is for atan, doesn't matter if nondim or not

a = GravAcc(t, y, we, Degree, Order, 1, 1, Clm, Slm);

if mu_Earth == 0
    a_Earth = [0;0;0];
else
    % if the Earth is included
    %the vector is in the Inertial Orbital Plane CS (has to be transformed to
    %MoonEq Inertial)
    [rM, ~] = Moon_crds_from_Earth(t*tunit,M0_moon);
    rE = -rM/runit; %nondimensionalized vector from Moon to Earth
    
    % STOPPED ROTATING FOR NOW AS I NEED EARTH TO BE IN OXY PLANE!
    %rotating rE into the Equatorial plane
    rot_z = [cos(pi), sin(pi), 0; -sin(pi), cos(pi), 0; 0, 0, 1];
    %rot_y = [cosd(-6.68), 0, -sind(-6.68); 0, 1, 0; sind(-6.68), 0, cosd(-6.68)];
    %first rotation around Z
    rE = rot_z*rE;
    %second rotation around Y
    %rE = rot_y*rE; %now rE is in Equatorial inertial
    % STOPPED ROTATING FOR NOW AS I NEED EARTH TO BE IN OXY PLANE!
    
    rE_SC = r - rE; %nondimensionalized vector from the Earth to the spacecraft
    
    q = r(1)/norm(rE) * (r(1)/norm(rE) - 2*rE(1)/norm(rE)) + r(2)/norm(rE) * (r(2)/norm(rE) - 2*rE(2)/norm(rE)) + r(3)/norm(rE) * (r(3)/norm(rE) - 2*rE(3)/norm(rE));
    %f = (3 + 3*q^3) / ((1+q)^1.5 + 1) * q; % ПОЧ БЫЛО ТАК????
    f = (3 + 3*q + q^2) / ((1+q)^1.5 + 1) * q;
    a_Earth = -(mu_Earth/mu_Moon)/norm(rE_SC)^3 * (r + f*rE); %acceleration vector of the spacecraft in Moon centered coordinate system
    %a_Earth = (mu_Earth/mu_Moon)*(-rE_SC/norm(rE_SC)^3-rE/norm(rE)^3);
end

a(4:6) = a(4:6) + a_Earth;

dydt = a;
end

