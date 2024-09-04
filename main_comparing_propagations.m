%Main file for comparing 2 propagations. 1st with Earth and/or harmonics,
%2nd without everything

clear;
clc;
global Mass_Earth;
global Mass_Moon;
global G;
global mu_Moon;
global mu_Earth;
global runit;
global tunit;
global vunit;
Mass_Earth = 5.9722*10^24;
Mass_Moon = 7.3477*10^22;
G = 6.674*10^-20;
mu_Moon = 4902.8001224453001;
mu_Earth = G*Mass_Earth;
%mu_Earth = 0;

% Units for nondimensionalization 
runit = 1738.0;
tunit = sqrt(1738.0^3/(mu_Moon));
vunit = runit/tunit;

%%Initial conditions
M0_moon = deg2rad(0); %M0 of Moon on its orbit
[~, r2] = calc_barycenter(0, M0_moon); %positions of the Earth and the Moon from the barycenter
r0 = 1800;
sv0_equatorial = Kepl2Dec(1,[r0/runit,0.03,45*pi/180,0*pi/180,90*pi/180,0*pi/180]);
sv0_nondim = inert2Moonequator_inert(r2, sv0_equatorial, 0);


% Chosen degree and order of Moon's gravitational field
Degree = 2;
Order = 0;

[Clm,Slm] = DenormCS('Moon_GRGM1200A.txt', 70);
Jcoeff = -Clm(2:end,1);

timespan = 2*pi*sqrt(384748^3/398600.44)/10/tunit; %devided!!!
timestep = 20/tunit;
number_of_points = ceil(timespan/timestep);


%% 1. Moon centered integration with the Earth and the Moon's harmonics

t_nondim = zeros(number_of_points,1);
inert_moon_sv_nondim = zeros(number_of_points,6);
inert_moon_sv_nondim(1,:) = sv0_nondim; %nondimensionalized state-vector for propagation
inert_moon_sv = zeros(number_of_points,6); % normal state-vector
inert_moon_sv(1, 1:3) = inert_moon_sv_nondim(1,1:3)*runit;
inert_moon_sv(1, 4:6) = inert_moon_sv_nondim(1,4:6)*vunit;
r1 = zeros(number_of_points,3);
r2 = zeros(number_of_points,3);
rot_r1 = zeros(number_of_points,3);
rot_r2 = zeros(number_of_points,3);
rot_sv = zeros(number_of_points,6);
mooneq_sv = zeros(number_of_points,6);
mooneq_inert_sv = zeros(number_of_points,6);
Kepler = zeros(number_of_points,6);
inert_cm_sv_crds = zeros(number_of_points,3);
inert_cm_sv_velocity = zeros(number_of_points,3);
inert_cm_sv = zeros(number_of_points,6);
radius = zeros(number_of_points,1);
speed = zeros(number_of_points,1);
Hamiltonian_bary = zeros(number_of_points,1);
Hamiltonian_Moon = zeros(number_of_points,1);
Hamiltonian_rot = zeros(number_of_points,1);
omega = zeros(number_of_points,1);

%nondimensionalized propagation
i = 1;
while t_nondim(i) < timespan
    
    [t_nondim(i+1), inert_moon_sv_nondim(i+1,:)] = RK4_onestep_Moon_cntrd_nondim(t_nondim(i), inert_moon_sv_nondim(i,:), timestep, M0_moon, Degree, Order, Clm, Slm);
    inert_moon_sv(i+1, 1:3) = inert_moon_sv_nondim(i+1,1:3)*runit;
    inert_moon_sv(i+1, 4:6) = inert_moon_sv_nondim(i+1,4:6)*vunit;
    i = i + 1;

end

% Deleting last elements, because those arrays are bigger than others
t_nondim(end) = [];
inert_moon_sv_nondim(end,:) = [];
inert_moon_sv(end,:) = [];

t = t_nondim * tunit; %normal time

for i = 1:number_of_points
    
    [r1(i,:), r2(i,:)] = calc_barycenter(t(i), M0_moon); %positions of the Earth and the Moon from the barycenter

    %transfering to the barycentric coordinate system
    [r_EM, v_EM] = Moon_crds_from_Earth(t(i), M0_moon); %Moon velocity in the Earth centered orbital ref frame
    v_cm = norm(v_EM)*norm(r1(i,:))/norm(r_EM); %speed of the center of mass with respect to the Earth
    v2 = v_EM/norm(v_EM)*(norm(v_EM) - v_cm); %Moon velocity from barycenter
    inert_cm_sv_crds(i,:) = [r2(i,1) + inert_moon_sv(i,1); r2(i,2) + inert_moon_sv(i,2); r2(i,3) + inert_moon_sv(i,3)];
    inert_cm_sv_velocity(i,:) = [v2(1) + inert_moon_sv(i,4); v2(2) + inert_moon_sv(i,5); v2(3) + inert_moon_sv(i,6)];
    inert_cm_sv(i,:) = [inert_cm_sv_crds(i,:), inert_cm_sv_velocity(i,:)];
    %transfering to rotating ref frame
    [rot_r1(i,:), rot_r2(i,:), rot_sv(i,:)] = inert2rot(r1(i,:), r2(i,:), v2, inert_cm_sv(i,:));
    
    %transfering to the Moon equator coordinate system
    %rotating
    [mooneq_sv(i,:)] = inert2Moonequator(r2(i,:), inert_moon_sv(i,:), 1);
    %non-rotating (inertial)
    [mooneq_inert_sv(i,:)] = inert2Moonequator_inert(r2(i,:), inert_moon_sv(i,:), 1);
    
    %calculating Kepler elements
    %Kepler(i,:) = Dec2Kepl(mu_Moon, mooneq_sv(i,:));
    Kepler(i,:) = Dec2Kepl(mu_Moon, mooneq_inert_sv(i,:)); %FROM MOON EQUATORIAL NON ROTATIONAL FRAME
    
    %abs of radius and Velocity  with respect to the Moon
    radius(i) = norm([inert_moon_sv(i,1), inert_moon_sv(i,2), inert_moon_sv(i,3)]);
    speed(i) = norm([inert_moon_sv(i,4), inert_moon_sv(i,5), inert_moon_sv(i,6)]);
    
    %calculating Hamiltonian in barycentric CS
    inert_r = [inert_cm_sv(i,1), inert_cm_sv(i,2), inert_cm_sv(i,3)];
    Hamiltonian_bary(i) = (inert_cm_sv(i,4)^2 + inert_cm_sv(i,5)^2 + inert_cm_sv(i,6)^2)/2 - mu_Earth/norm(inert_r-r1(i,:)) - mu_Moon/norm(inert_r-r2(i,:));
    
    %calculating Hamiltonian in Moon inertial CS
    inert_r = [inert_moon_sv(i,1), inert_moon_sv(i,2), inert_moon_sv(i,3)];
    Hamiltonian_Moon(i) = (inert_moon_sv(i,4)^2 + inert_moon_sv(i,5)^2 + inert_moon_sv(i,6)^2)/2 - mu_Earth/norm(inert_r+r_EM) - mu_Moon/norm(inert_r);

    %calculating Hamiltonian in barycentric rotating CS
    inert_r = [rot_sv(i,1), rot_sv(i,2), rot_sv(i,3)];
    omega(i) = norm(v2) / norm(r2(i,:)); % omega of Earth-Moon system
    Hamiltonian_rot(i) = (rot_sv(i,4)^2 + rot_sv(i,5)^2 + rot_sv(i,6)^2)/2 - mu_Earth/norm(inert_r-rot_r1(i,:)) - mu_Moon/norm(inert_r-rot_r2(i,:)) - omega(i)^2*(rot_sv(i,1)^2 + rot_sv(i,2)^2)/2;

end

%% 2. Moon centered integration w/o the Earth and the Moon's harmonics

%new variables are called ..._initial
inert_moon_sv_initial = zeros(number_of_points,6); % normal state-vector
inert_moon_sv_initial(1, 1:3) = inert_moon_sv_nondim(1,1:3)*runit;
inert_moon_sv_initial(1, 4:6) = inert_moon_sv_nondim(1,4:6)*vunit;
mooneq_inert_sv_initial = zeros(number_of_points,6);
Kepler_initial = zeros(number_of_points,6);
radius_initial = zeros(number_of_points,1);
speed_initial = zeros(number_of_points,1);

mu_Earth = 0;
% Chosen degree and order of Moon's gravitational field
Degree = 0;
Order = 0;

%nondimensionalized propagation
i = 1;
while t_nondim(i) < timespan
    
    [t_nondim(i+1), inert_moon_sv_nondim(i+1,:)] = RK4_onestep_Moon_cntrd_nondim(t_nondim(i), inert_moon_sv_nondim(i,:), timestep, M0_moon, Degree, Order, Clm, Slm);
    inert_moon_sv_initial(i+1, 1:3) = inert_moon_sv_nondim(i+1,1:3)*runit;
    inert_moon_sv_initial(i+1, 4:6) = inert_moon_sv_nondim(i+1,4:6)*vunit;
    i = i + 1;

end

% Deleting last elements, because those arrays are bigger than others
t_nondim(end) = [];
inert_moon_sv_nondim(end,:) = [];
inert_moon_sv_initial(end,:) = [];

t = t_nondim * tunit; %normal time

for i = 1:number_of_points
    
    [r1(i,:), r2(i,:)] = calc_barycenter(t(i), M0_moon); %positions of the Earth and the Moon from the barycenter

    %transfering to the barycentric coordinate system
    [r_EM, v_EM] = Moon_crds_from_Earth(t(i), M0_moon); %Moon velocity in the Earth centered orbital ref frame
    v_cm = norm(v_EM)*norm(r1(i,:))/norm(r_EM); %speed of the center of mass with respect to the Earth
    v2 = v_EM/norm(v_EM)*(norm(v_EM) - v_cm); %Moon velocity from barycenter
    inert_cm_sv_crds(i,:) = [r2(i,1) + inert_moon_sv_initial(i,1); r2(i,2) + inert_moon_sv_initial(i,2); r2(i,3) + inert_moon_sv_initial(i,3)];
    inert_cm_sv_velocity(i,:) = [v2(1) + inert_moon_sv_initial(i,4); v2(2) + inert_moon_sv_initial(i,5); v2(3) + inert_moon_sv_initial(i,6)];
    inert_cm_sv(i,:) = [inert_cm_sv_crds(i,:), inert_cm_sv_velocity(i,:)];
    %transfering to rotating ref frame
    [rot_r1(i,:), rot_r2(i,:), rot_sv(i,:)] = inert2rot(r1(i,:), r2(i,:), v2, inert_cm_sv(i,:));
    
    %transfering to the Moon equator coordinate system
    %rotating
    [mooneq_sv(i,:)] = inert2Moonequator(r2(i,:), inert_moon_sv_initial(i,:), 1);
    %non-rotating (inertial)
    [mooneq_inert_sv_initial(i,:)] = inert2Moonequator_inert(r2(i,:), inert_moon_sv_initial(i,:), 1);
    
    %calculating Kepler elements
    %Kepler(i,:) = Dec2Kepl(mu_Moon, mooneq_sv(i,:));
    Kepler_initial(i,:) = Dec2Kepl(mu_Moon, mooneq_inert_sv_initial(i,:)); %FROM MOON EQUATORIAL NON ROTATIONAL FRAME
    
    %abs of radius and Velocity  with respect to the Moon
    radius_initial(i) = norm([inert_moon_sv_initial(i,1), inert_moon_sv_initial(i,2), inert_moon_sv_initial(i,3)]);
    speed_initial(i) = norm([inert_moon_sv_initial(i,4), inert_moon_sv_initial(i,5), inert_moon_sv_initial(i,6)]);
end

%% Plots of differences

figure;
tiledlayout(2,3)

nexttile
plot(t/86400,Kepler(:,1) - Kepler_initial(:,1));
title('delta a(t)');

nexttile
plot(t/86400,Kepler(:,2) - Kepler_initial(:,2));
title('delta e(t)');

nexttile
plot(t/86400,rad2deg(Kepler(:,3) - Kepler_initial(:,3)));
title('delta i(t)');

nexttile
plot(t/86400,rad2deg(Kepler(:,4) - Kepler_initial(:,4)));
title('delta omega(t)');

nexttile
plot(t/86400,rad2deg(Kepler(:,5) - Kepler_initial(:,5)));
title('delta RAAN(t)');

nexttile
plot(t/86400,rad2deg(Kepler(:,6) - Kepler_initial(:,6)));
title('delta M(t)');

figure;
tiledlayout(2,1)
nexttile
plot(t/86400, radius - radius_initial);
title('|r|(t)');

nexttile
plot(t/86400, speed - speed_initial);
title('|V|(t)');

% figure;
% plot(t/86400, Hamiltonian_bary);
% title('H(t) barycentric');
% 
% figure;
% plot(t/86400, Hamiltonian_Moon);
% title('H(t) Moon CS');
% 
% figure;
% plot(t/86400, Hamiltonian_rot);
% title('H(t) rotating CS');