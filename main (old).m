%Main file for propagation and plots

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

% Units for nondimensionalization 
runit = 1738.0;
tunit = sqrt(1738.0^3/(mu_Moon));
vunit = runit/tunit;

%%Initial conditions
M0_moon = deg2rad(0); %M0 of Moon on its orbit
%sv0 = [320000; 5000; 0; 0; 1; 0]; %Center of mass reference frame
r0 = 10000;
sv0 = sv0_moon_circular_orbit(r0, M0_moon);
%sv0(5) = sv0(5) + 0.02;
%sv0 = sv0_earth_circular_orbit(r0, M0_moon);

%checking if the initial conditions are correct
[r1(:), r2(:)] = calc_barycenter(0, M0_moon); %positions of the Earth and the Moon from the barycenter   
%transfering to the Moon equator coordinate system
[r_EM, v_EM] = Moon_crds_from_Earth(0, M0_moon); %Moon velocity in the Earth centered orbital ref frame
v_cm = norm(v_EM)*norm(r1(:))/norm(r_EM); %speed of the center of mass with respect to the Earth
v2 = v_EM/norm(v_EM)*(norm(v_EM) - v_cm); %Moon velocity from barycenter
[mooneq_sv(:)] = inert2Moonequator(r2(:), sv0(:), 1);

%calculating Kepler elements
Kepler(:) = Dec2Kepl(Mass_Moon*G, mooneq_sv(:));

%% Propagation (3 bodies)

number_of_points = 1000;
timespan = linspace(0,2*pi*sqrt(384748^3/398600.44),number_of_points);
[t, inert_sv] = ode113(@(t, inert_sv) derivatives_cm(t,inert_sv, M0_moon), timespan, sv0);

r1 = zeros(number_of_points,3);
r2 = zeros(number_of_points,3);
rot_r1 = zeros(number_of_points,3);
rot_r2 = zeros(number_of_points,3);
rot_sv = zeros(number_of_points,6);
mooneq_sv = zeros(number_of_points,6);
Kepler = zeros(number_of_points,6);
radius = zeros(number_of_points,1);
speed = zeros(number_of_points,1);
Hamiltonian = zeros(number_of_points,1);

for i = 1:number_of_points
    [r1(i,:), r2(i,:)] = calc_barycenter(t(i), M0_moon); %positions of the Earth and the Moon from the barycenter   

    %transfering to rotating ref frame
    [rot_r1(i,:), rot_r2(i,:), rot_sv(i,:)] = inert2rot(r1(i,:), r2(i,:), inert_sv(i,:));
    %transfering to the Moon equator coordinate system
    [r_EM, v_EM] = Moon_crds_from_Earth(t(i), M0_moon); %Moon velocity in the Earth centered orbital ref frame
    v_cm = norm(v_EM)*norm(r1(i,:))/norm(r_EM); %speed of the center of mass with respect to the Earth
    v2 = v_EM/norm(v_EM)*(norm(v_EM) - v_cm); %Moon velocity from barycenter
    [mooneq_sv(i,:)] = inert2Moonequator(r2(i,:), v2, inert_sv(i,:), 1);
    
    %calculating Kepler elements
    Kepler(i,:) = Dec2Kepl(Mass_Moon*G, mooneq_sv(i,:));
    
    %abs of radius and Velocity with respect to the Moon
    radius(i) = norm([mooneq_sv(i,1), mooneq_sv(i,2), mooneq_sv(i,3)]);
    speed(i) = norm([mooneq_sv(i,4), mooneq_sv(i,5), mooneq_sv(i,6)]);
    %calculating Hamiltonian
    inert_r = [inert_sv(i,1), inert_sv(i,2), inert_sv(i,3)];
    Hamiltonian(i) = (inert_sv(i,4)^2 + mooneq_sv(i,5)^2 + mooneq_sv(i,6)^2)/2 - Mass_Earth*G/(norm(inert_r-r1(i,:))) - Mass_Moon*G/(norm(inert_r-r2(i,:)));

end

% Plots
figure(1)
tiledlayout(2,1)

nexttile
plot(inert_sv(:,1), inert_sv(:,2));
hold on;
plot(r1(:,1),r1(:,2));
plot(r2(:,1),r2(:,2));
axis equal;
%plot([r1(1,1) r2(1,1)], [r1(1,2) r2(1,2)],'r--');
%plot([r1(250,1) r2(250,1)], [r1(250,2) r2(250,2)],'r--');
%plot([r1(750,1) r2(750,1)], [r1(750,2) r2(750,2)],'r--');
hold off;

nexttile
plot(rot_sv(:,1), rot_sv(:,2));
hold on;
plot(rot_r1(:,1), 0, 'o');
plot(rot_r2(:,1), 0, 'o');
axis equal;
hold off;

figure(2)
[X_sph, Y_sph, Z_sph] = sphere(50);
X_sph = 1737.4*X_sph;
Y_sph = 1737.4*Y_sph;
Z_sph = 1737.4*Z_sph;
plot3(X_sph, Y_sph, Z_sph);
hold on;
plot3(mooneq_sv(:,1), mooneq_sv(:,2), mooneq_sv(:,3));
axis equal;

figure(3)
tiledlayout(2,3)

nexttile
plot(t,Kepler(:,1));
title('a(t)');

nexttile
plot(t,Kepler(:,2));
title('e(t)');

nexttile
plot(t,rad2deg(Kepler(:,3)));
title('i(t)');

nexttile
plot(t,rad2deg(Kepler(:,4)));
title('omega(t)');

nexttile
plot(t,rad2deg(Kepler(:,5)));
title('RAAN(t)');

nexttile
plot(t,rad2deg(Kepler(:,6)));
title('M(t)');

figure(4)
tiledlayout(2,1)
nexttile
plot(t, radius);
title('|r|(t)');

nexttile
plot(t, speed);
title('|V|(t)');

figure(5)
plot(t, Hamiltonian);
title('H(t)');

%% Propagation RK4 (3 bodies)

timespan = 2*pi*sqrt(384748^3/398600.44);
timestep = 120;
number_of_points = ceil(timespan/timestep);

t = zeros(number_of_points,1);
inert_sv = zeros(number_of_points,6);
inert_sv(1,:) = sv0;
r1 = zeros(number_of_points,3);
r2 = zeros(number_of_points,3);
rot_r1 = zeros(number_of_points,3);
rot_r2 = zeros(number_of_points,3);
rot_sv = zeros(number_of_points,6);
mooneq_sv = zeros(number_of_points,6);
Kepler = zeros(number_of_points,6);
radius = zeros(number_of_points,1);
speed = zeros(number_of_points,1);
speed_barycenter = zeros(number_of_points,1);
Hamiltonian = zeros(number_of_points,1);

i = 1;

while t(i) < timespan
    
    
    [r1(i,:), r2(i,:)] = calc_barycenter(t(i), M0_moon); %positions of the Earth and the Moon from the barycenter

    [r_EM, v_EM] = Moon_crds_from_Earth(t(i), M0_moon); %Moon velocity in the Earth centered orbital ref frame
    v_cm = norm(v_EM)*norm(r1(i,:))/norm(r_EM); %speed of the center of mass with respect to the Earth
    v2 = v_EM/norm(v_EM)*(norm(v_EM) - v_cm); %Moon velocity from barycenter
    %v1 = -v_EM/norm(v_EM)*v_cm; %Earth velocity from barycenter
    %transfering to rotating ref frame
    [rot_r1(i,:), rot_r2(i,:), rot_sv(i,:)] = inert2rot(r1(i,:), r2(i,:), inert_sv(i,:));
    %transfering to the Moon equator coordinate system
    [mooneq_sv(i,:)] = inert2Moonequator(r2(i,:), v2, inert_sv(i,:), 1);
    %[mooneq_sv(i,:)] = inert2Earth(r1(i,:), v1, inert_sv(i,:));
    
    %calculating Kepler elements
    Kepler(i,:) = Dec2Kepl(Mass_Moon*G, mooneq_sv(i,:));
    %Kepler(i,:) = Dec2Kepl(Mass_Earth*G, mooneq_sv(i,:));
    
    %abs of radius and Velocity  with respect to the Moon
    radius(i) = norm([mooneq_sv(i,1), mooneq_sv(i,2), mooneq_sv(i,3)]);
    speed(i) = norm([mooneq_sv(i,4), mooneq_sv(i,5), mooneq_sv(i,6)]);
    speed_barycenter(i) = norm([inert_sv(i,4), inert_sv(i,5), inert_sv(i,6)]);
    
    %calculating Hamiltonian
    inert_r = [inert_sv(i,1), inert_sv(i,2), inert_sv(i,3)];
    Hamiltonian(i) = (inert_sv(i,4)^2 + inert_sv(i,5)^2 + inert_sv(i,6)^2)/2 - Mass_Earth*G/norm(inert_r-r1(i,:)) - Mass_Moon*G/norm(inert_r-r2(i,:));

    [t(i+1), inert_sv(i+1,:)] = RK4_onestep(t(i), inert_sv(i,:), timestep, M0_moon, 1);
    i = i + 1;
    
end

% Deleting last elements, because those arrays are bigger than others
t(end) = [];
inert_sv(end) = [];

% Plots
figure(1)
tiledlayout(2,1)

nexttile
plot(inert_sv(:,1), inert_sv(:,2));
hold on;
plot(r1(:,1),r1(:,2));
plot(r2(:,1),r2(:,2));
axis equal;
%plot([r1(1,1) r2(1,1)], [r1(1,2) r2(1,2)],'r--');
%plot([r1(250,1) r2(250,1)], [r1(250,2) r2(250,2)],'r--');
%plot([r1(750,1) r2(750,1)], [r1(750,2) r2(750,2)],'r--');
hold off;

nexttile
plot(rot_sv(:,1), rot_sv(:,2));
hold on;
plot(rot_r1(:,1), 0, 'o');
plot(rot_r2(:,1), 0, 'o');
axis equal;
hold off;

figure(2)
[X_sph, Y_sph, Z_sph] = sphere(50);
X_sph = 1737.4*X_sph;
Y_sph = 1737.4*Y_sph;
Z_sph = 1737.4*Z_sph;
plot3(X_sph, Y_sph, Z_sph);
hold on;
plot3(mooneq_sv(:,1), mooneq_sv(:,2),mooneq_sv(:,3));
axis equal;

figure(3)
tiledlayout(2,3)

nexttile
plot(t,Kepler(:,1));
title('a(t)');

nexttile
plot(t,Kepler(:,2));
title('e(t)');

nexttile
plot(t,rad2deg(Kepler(:,3)));
title('i(t)');

nexttile
plot(t,rad2deg(Kepler(:,4)));
title('omega(t)');

nexttile
plot(t,rad2deg(Kepler(:,5)));
title('RAAN(t)');

nexttile
plot(t,rad2deg(Kepler(:,6)));
title('M(t)');

figure(4)
tiledlayout(2,1)
nexttile
plot(t, radius);
title('|r|(t)');

nexttile
plot(t, speed);
title('|V|(t)');

figure(5)
plot(t, Hamiltonian);
title('H(t)');

figure(6)
plot(t, speed_barycenter);
title('|V|_bary(t)');


%% Propagation (only the Moon)

number_of_points = 2000;
timespan = linspace(0,2*pi*sqrt(384748^3/398600.44),number_of_points);
sv0 = [r0; 0; 0; 0; sqrt(G*Mass_Moon/r0); 0];
[t, inert_moon_sv] = ode113(@derivatives_moon, timespan, sv0);

r1 = zeros(number_of_points,3);
r2 = zeros(number_of_points,3);
rot_r1 = zeros(number_of_points,3);
rot_r2 = zeros(number_of_points,3);
rot_sv = zeros(number_of_points,6);
mooneq_sv = zeros(number_of_points,6);
Kepler = zeros(number_of_points,6);
inert_cm_sv_crds = zeros(number_of_points,3);
inert_cm_sv_velocity = zeros(number_of_points,3);
inert_cm_sv = zeros(number_of_points,6);
radius = zeros(number_of_points,1);
speed = zeros(number_of_points,1);

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
    [rot_r1(i,:), rot_r2(i,:), rot_sv(i,:)] = inert2rot(r1(i,:), r2(i,:), inert_cm_sv(i,:));
    %transfering to the Moon equator coordinate system
    [mooneq_sv(i,:)] = inert2Moonequator(r2(i,:), v2, inert_cm_sv(i,:), 1);
    
    %calculating Kepler elements
    Kepler(i,:) = fromDec2Kepl(Mass_Moon*G, mooneq_sv(i,:));
    
    %abs of radius and Velocity  with respect to the Moon
    radius(i) = norm([inert_moon_sv(i,1), inert_moon_sv(i,2), inert_moon_sv(i,3)]);
    speed(i) = norm([inert_moon_sv(i,4), inert_moon_sv(i,5), inert_moon_sv(i,6)]);
end

% Plots
figure(1)
tiledlayout(2,1)

nexttile
plot(inert_cm_sv(:,1), inert_cm_sv(:,2));
hold on;
plot(r1(:,1),r1(:,2));
plot(r2(:,1),r2(:,2));
axis equal;
hold off;

nexttile
plot(rot_sv(:,1), rot_sv(:,2));
hold on;
plot(rot_r1(:,1), 0, 'o');
plot(rot_r2(:,1), 0, 'o');
axis equal;
hold off;

figure(2)
[X_sph, Y_sph, Z_sph] = sphere(50);
X_sph = 1737.4*X_sph;
Y_sph = 1737.4*Y_sph;
Z_sph = 1737.4*Z_sph;
plot3(X_sph, Y_sph, Z_sph);
hold on;
plot3(mooneq_sv(:,1), mooneq_sv(:,2),mooneq_sv(:,3));
axis equal;

figure(3)
tiledlayout(2,3)

nexttile
plot(t,Kepler(:,1));
title('a(t)');

nexttile
plot(t,Kepler(:,2));
title('e(t)');

nexttile
plot(t,rad2deg(Kepler(:,3)));
title('i(t)');

nexttile
plot(t,rad2deg(Kepler(:,4)));
title('omega(t)');

nexttile
plot(t,rad2deg(Kepler(:,5)));
title('RAAN(t)');

nexttile
plot(t,rad2deg(Kepler(:,6)));
title('M(t)');

figure(4)
tiledlayout(2,1)
nexttile
plot(t, radius);
title('|r|(t)');

nexttile
plot(t, speed);
title('|V|(t)');

%% Propagation RK4 (only the Moon)

timespan = 2*pi*sqrt(384748^3/398600.44);
timestep = 120;
number_of_points = ceil(timespan/timestep);
sv0 = [r0; 0; 0; 0; sqrt(G*Mass_Moon/r0); 0];

t = zeros(number_of_points,1);
inert_moon_sv = zeros(number_of_points,6);
inert_moon_sv(1,:) = sv0;
r1 = zeros(number_of_points,3);
r2 = zeros(number_of_points,3);
rot_r1 = zeros(number_of_points,3);
rot_r2 = zeros(number_of_points,3);
rot_sv = zeros(number_of_points,6);
mooneq_sv = zeros(number_of_points,6);
Kepler = zeros(number_of_points,6);
inert_cm_sv_crds = zeros(number_of_points,3);
inert_cm_sv_velocity = zeros(number_of_points,3);
inert_cm_sv = zeros(number_of_points,6);
radius = zeros(number_of_points,1);
speed = zeros(number_of_points,1);

i = 1;

while t(i) < timespan
    
    [r1(i,:), r2(i,:)] = calc_barycenter(t(i), M0_moon); %positions of the Earth and the Moon from the barycenter

    %transfering to the barycentric coordinate system
    [r_EM, v_EM] = Moon_crds_from_Earth(t(i), M0_moon); %Moon velocity in the Earth centered orbital ref frame
    v_cm = norm(v_EM)*norm(r1(i,:))/norm(r_EM); %speed of the center of mass with respect to the Earth
    v2 = v_EM/norm(v_EM)*(norm(v_EM) - v_cm); %Moon velocity from barycenter
    inert_cm_sv_crds(i,:) = [r2(i,1) + inert_moon_sv(i,1); r2(i,2) + inert_moon_sv(i,2); r2(i,3) + inert_moon_sv(i,3)];
    inert_cm_sv_velocity(i,:) = [v2(1) + inert_moon_sv(i,4); v2(2) + inert_moon_sv(i,5); v2(3) + inert_moon_sv(i,6)];
    inert_cm_sv(i,:) = [inert_cm_sv_crds(i,:), inert_cm_sv_velocity(i,:)];
    %transfering to rotating ref frame
    [rot_r1(i,:), rot_r2(i,:), rot_sv(i,:)] = inert2rot(r1(i,:), r2(i,:), inert_cm_sv(i,:));
    %transfering to the Moon equator coordinate system
    [mooneq_sv(i,:)] = inert2Moonequator(r2(i,:), v2, inert_cm_sv(i,:), 1);
    
    %calculating Kepler elements
    Kepler(i,:) = Dec2Kepl(Mass_Moon*G, mooneq_sv(i,:));
    
    %abs of radius and Velocity  with respect to the Moon
    radius(i) = norm([inert_moon_sv(i,1), inert_moon_sv(i,2), inert_moon_sv(i,3)]);
    speed(i) = norm([inert_moon_sv(i,4), inert_moon_sv(i,5), inert_moon_sv(i,6)]);
    
    [t(i+1), inert_moon_sv(i+1,:)] = RK4_onestep(t(i), inert_moon_sv(i,:), timestep, M0_moon, 0);
    i = i + 1;

end

% Deleting last elements, because those arrays are bigger than others
t(end) = [];
inert_moon_sv(end) = [];

% Plots
figure(1)
tiledlayout(2,1)

nexttile
plot(inert_cm_sv(:,1), inert_cm_sv(:,2));
hold on;
plot(r1(:,1),r1(:,2));
plot(r2(:,1),r2(:,2));
axis equal;
hold off;

nexttile
plot(rot_sv(:,1), rot_sv(:,2));
hold on;
plot(rot_r1(:,1), 0, 'o');
plot(rot_r2(:,1), 0, 'o');
axis equal;
hold off;

figure(2)
[X_sph, Y_sph, Z_sph] = sphere(50);
X_sph = 1737.4*X_sph;
Y_sph = 1737.4*Y_sph;
Z_sph = 1737.4*Z_sph;
plot3(X_sph, Y_sph, Z_sph);
hold on;
plot3(mooneq_sv(:,1), mooneq_sv(:,2),mooneq_sv(:,3));
axis equal;

figure(3)
tiledlayout(2,3)

nexttile
plot(t,Kepler(:,1));
title('a(t)');

nexttile
plot(t,Kepler(:,2));
title('e(t)');

nexttile
plot(t,rad2deg(Kepler(:,3)));
title('i(t)');

nexttile
plot(t,rad2deg(Kepler(:,4)));
title('omega(t)');

nexttile
plot(t,rad2deg(Kepler(:,5)));
title('RAAN(t)');

nexttile
plot(t,rad2deg(Kepler(:,6)));
title('M(t)');

figure(4)
tiledlayout(2,1)
nexttile
plot(t, radius);
title('|r|(t)');

nexttile
plot(t, speed);
title('|V|(t)');


%% Comparing math models

number_of_points = 1000;
timespan = linspace(0,2*pi*sqrt(384748^3/398600.44),number_of_points);
delta_r = [];
initial_r0 = [];
i = 1;

for r0 = 2000:200:20000
    %propagate with the Earth in barycentric cs
    sv0 = sv0_moon_circular_orbit(r0, M0_moon);
    [~, with_Earth_sv] = ode113(@(t, with_Earth_sv) derivatives_cm(t,with_Earth_sv, M0_moon), timespan, sv0);
    
    %propagate without the Earth in selenocentric cs
    sv0 = [r0; 0; 0; 0; sqrt(G*Mass_Moon/r0); 0];
    [t, without_Earth_sv] = ode113(@derivatives_moon, timespan, sv0);
    %transfer to barycenter
    [~, r_Moon_end] = calc_barycenter(t(end), M0_moon);
    
    %difference in r between the different approaches
    delta_r(i) = sqrt((with_Earth_sv(end,1) - (without_Earth_sv(end,1)+r_Moon_end(1)))^2 + (with_Earth_sv(end,2) - (without_Earth_sv(end,2)+r_Moon_end(2)))^2 + (with_Earth_sv(end,3) - (without_Earth_sv(end,3)+r_Moon_end(3)))^2);
    initial_r0(i) = r0;
    i = i+1;
end
figure(3)
plot(initial_r0, delta_r);

%% Moon centered integration with the Earth (nondimensionalized)

timespan = 2*pi*sqrt(384748^3/398600.44)/tunit;
timestep = 120/tunit;
number_of_points = ceil(timespan/timestep);
sv0_nondim = [r0/runit; 0; 0; 0; sqrt(mu_Moon/r0)/vunit; 0];

t_nondim = zeros(number_of_points,1);
t = zeros(number_of_points,1);
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
Kepler = zeros(number_of_points,6);
inert_cm_sv_crds = zeros(number_of_points,3);
inert_cm_sv_velocity = zeros(number_of_points,3);
inert_cm_sv = zeros(number_of_points,6);
radius = zeros(number_of_points,1);
speed = zeros(number_of_points,1);
Hamiltonian_bary = zeros(number_of_points,1);
Hamiltonian_Moon = zeros(number_of_points,1);
Hamiltonian_rot = zeros(number_of_points,1);

%nondimensionalized propagation
i = 1;
while t_nondim(i) < timespan
    
    [t_nondim(i+1), inert_moon_sv_nondim(i+1,:)] = RK4_onestep_Moon_cntrd_nondim(t_nondim(i), inert_moon_sv_nondim(i,:), timestep, M0_moon);
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
    [rot_r1(i,:), rot_r2(i,:), rot_sv(i,:)] = inert2rot(r1(i,:), r2(i,:), inert_cm_sv(i,:));
    %transfering to the Moon equator coordinate system
    [mooneq_sv(i,:)] = inert2Moonequator(r2(i,:), v2, inert_cm_sv(i,:), 1);
    
    %calculating Kepler elements
    Kepler(i,:) = Dec2Kepl(mu_Moon, mooneq_sv(i,:));
    
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
    Hamiltonian_rot(i) = (rot_sv(i,4)^2 + rot_sv(i,5)^2 + rot_sv(i,6)^2)/2 - mu_Earth/norm(inert_r-rot_r1(i,:)) - mu_Moon/norm(inert_r-rot_r2(i,:));

end

% Plots
figure(1)
tiledlayout(2,1)

nexttile
plot(inert_cm_sv(:,1), inert_cm_sv(:,2));
hold on;
plot(r1(:,1),r1(:,2));
plot(r2(:,1),r2(:,2));
axis equal;
hold off;

nexttile
plot(rot_sv(:,1), rot_sv(:,2));
hold on;
plot(rot_r1(:,1), 0, 'o');
plot(rot_r2(:,1), 0, 'o');
axis equal;
hold off;

figure(2)
[X_sph, Y_sph, Z_sph] = sphere(50);
X_sph = 1737.4*X_sph;
Y_sph = 1737.4*Y_sph;
Z_sph = 1737.4*Z_sph;
plot3(X_sph, Y_sph, Z_sph);
hold on;
plot3(mooneq_sv(:,1), mooneq_sv(:,2),mooneq_sv(:,3));
axis equal;

figure(3)
tiledlayout(2,3)

nexttile
plot(t/86400,Kepler(:,1));
title('a(t)');

nexttile
plot(t/86400,Kepler(:,2));
title('e(t)');

nexttile
plot(t/86400,rad2deg(Kepler(:,3)));
title('i(t)');

nexttile
plot(t/86400,rad2deg(Kepler(:,4)));
title('omega(t)');

nexttile
plot(t/86400,rad2deg(Kepler(:,5)));
title('RAAN(t)');

nexttile
plot(t/86400,rad2deg(Kepler(:,6)));
title('M(t)');

figure(4)
tiledlayout(2,1)
nexttile
plot(t/86400, radius);
title('|r|(t)');

nexttile
plot(t/86400, speed);
title('|V|(t)');

figure(5)
plot(t/86400, Hamiltonian_bary);
title('H(t) barycentric');

figure(6)
plot(t/86400, Hamiltonian_Moon);
title('H(t) Moon CS');

figure;
plot(t/86400, Hamiltonian_rot);
title('H(t) rotating CS');