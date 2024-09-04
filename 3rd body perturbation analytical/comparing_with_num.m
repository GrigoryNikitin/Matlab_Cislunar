% Here I am comparing numerical integration with the 3rd body and the
% analytical stuff I got

%% Initial

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
global we;
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

%normilized r_Moon and mu_Moon
re = 1;
mu = 1;

% Chosen degree and order of Moon's gravitational field
Degree = 0;
Order = 0;

% ODE tolerances
tol = 1e-13;
quadtol = 1e-6;
optn = odeset('RelTol',tol,'AbsTol',tol*1e-3);

[Clm,Slm] = DenormCS('Moon_GRGM1200A.txt', 70);
Jcoeff = -Clm(2:end,1);

%%Initial conditions
% Moon rotation rate
we = 2*pi/(27.3220*86400/tunit); %rad/sec nondimensionalized
M0_moon = deg2rad(0); %M0 of Moon on its orbit
[~, r2] = calc_barycenter(0, M0_moon); %positions of the Earth and the Moon from the barycenter
r0 = 10000;
sv0_equatorial_inert = Kepl2Dec(mu,[r0/runit,0.05,30*pi/180,30*pi/180,30*pi/180,60*pi/180]); % state vector in the MoonEq_inertial (not pointing at the Earth!) coordinate system

%% Numerical

finaltime = 2*pi*sqrt(384748^3/398600.44)/4/tunit; %nondimensionalized
number_of_points = 2000;
timespan = linspace(0,finaltime,number_of_points);

% arrays for integration in Moon Equatorial inertial CS
inert_moon_sv = zeros(number_of_points,6); % normal state-vector
%Orbital Inertial CS
orbital_moon_sv = zeros(number_of_points,6);

Kepler = zeros(number_of_points,6);
Delaunay = zeros(number_of_points,6);
Delaunay_nondim = zeros(number_of_points,6);

%ODE integration
optn = odeset('RelTol',tol,'AbsTol',tol*1e-10);
[tspan,inert_moon_sv_nondim] = ode113(@(tspan, inert_moon_sv_nondim) derivatives_Moon_cntrd_nondim(tspan, inert_moon_sv_nondim, M0_moon, Degree, Order, Clm, Slm),timespan,sv0_equatorial_inert,optn);
%ODE integration ^^^

t = tspan * tunit; %normal time

for i = 1:number_of_points
    
    inert_moon_sv(i, 1:3) = inert_moon_sv_nondim(i,1:3)*runit;
    inert_moon_sv(i, 4:6) = inert_moon_sv_nondim(i,4:6)*vunit;
    
    %calculating Kepler elements
    Kepler(i,:) = Dec2Kepl(mu_Moon, inert_moon_sv(i,:)); %FROM MOON EQUATORIAL INERTIAL FRAME
    Delaunay(i,:) = Kepl2Del(mu_Moon, Kepler(i,:), false);
    Delaunay_nondim(i,:) = Kepl2Del(1, [Kepler(i,1)/runit, Kepler(i,2:6)], false);
end

%Deleting leaps in Kepler elements
    for i = 1:number_of_points
        for j = 3:5
            if (Kepler(i,j) > 10*pi/180) && (abs(Kepler(i,j) - 2*pi) < 10*pi/180)
                Kepler(i,j) = Kepler(i,j) - 2*pi;
            end
        end
    end
    
%% Analytical

Delaunay_an0 = Kepl2Del(1, [Kepler(1,1)/runit, Kepler(1,2:6)], false);

% Earth's Kepler elements are calculated below
%the vector is in the Inertial Orbital Plane CS (has to be transformed to
%MoonEq Inertial)
[rM, vM] = Moon_crds_from_Earth(0,M0_moon);
rE = -rM/runit; %nondimensionalized vector from Moon to Earth
vE = -vM/vunit; %nondimensionalized velocity of Earth with respect to Moon

% STOPPED ROTATING FOR NOW AS I NEED EARTH TO BE IN OXY PLANE!
% %rotating rE and vE into the Equatorial plane
rot_z = [cos(pi), sin(pi), 0; -sin(pi), cos(pi), 0; 0, 0, 1];
% rot_y = [cosd(-6.68), 0, -sind(-6.68); 0, 1, 0; sind(-6.68), 0, cosd(-6.68)];
rE = rot_z*rE;
% rE = rot_y*rE; %now rE is in Equatorial inertial
vE = rot_z*vE;
% vE = rot_y*vE; % now vE is in Equatorial inertial
% STOPPED ROTATING FOR NOW AS I NEED EARTH TO BE IN OXY PLANE!

Kepler_Earth = Dec2Kepl(mu_Earth/mu_Moon, [rE, vE]); %Should use the same mu??????

Delaunay_mean0 = Osc2Mean_thirdbody(Delaunay_an0, Kepler_Earth, 1);
Delaunay_mean00 = Osc2Mean_thirdbody_new_one(Delaunay_an0, Kepler_Earth, 1);

%Checking if Osc2Mean and Mean2Osc gives the same
% DO SUCCESSIVE APPROX
aa = Mean2Osc_thirdbody(Delaunay_mean0, Kepler_Earth);
aa00 = Mean2Osc_thirdbody_new_one(Delaunay_mean00, Kepler_Earth);
disp(Delaunay_an0 - aa);
disp(Delaunay_an0 - aa00);
% if aa == Delaunay_an0 then we are good


% Propagate Mean States
[tspan, Delaunay_mean] = ode113(@Secular_rates_third_body,timespan,Delaunay_mean00,optn);
 
Kepler_mean = zeros(length(tspan),6);
for i = 1:length(tspan)
    Delaunay_mean(i,1) = rem(Delaunay_mean(i,1), 2*pi);
    Kepler_mean(i,:) = Kepl2Del( mu, Delaunay_mean(i,1:6)', true )'; %mean Kepler elements
end

% Mean to Osculating Transformation
Delaunay_an = zeros(length(tspan),6);
Earth_sv = zeros(length(tspan),6);
Kepler_Earth = zeros(length(tspan),6);
Delaunay_mean_num = zeros(length(tspan),6);

for i = 1:length(tspan)
    % Earth's Kepler elements are calculated below
    %the vector is in the Inertial Orbital Plane CS (has to be transformed to
    %MoonEq Inertial)
    [rM, vM] = Moon_crds_from_Earth(tspan(i)*tunit, M0_moon);
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
    
    Earth_sv(i,:) = [rE', vE'];
    %disp(rE(1) * vE(1) + rE(2) * vE(2));
    Kepler_Earth(i,:) = Dec2Kepl(mu_Earth/mu_Moon, Earth_sv(i,:));
    Delaunay_an(i,:) = Mean2Osc_thirdbody_new_one(Delaunay_mean(i,:)', Kepler_Earth(i,:)');
    Delaunay_an(i,1) = rem(Delaunay_an(i,1), 2*pi);
    Delaunay_mean_num(i,:) = Osc2Mean_thirdbody_new_one(Delaunay_nondim(i,:)', Kepler_Earth(i,:), 0);
end

Kepler_an = zeros(length(tspan),6);
inert_moon_sv_nondim_an = zeros(length(tspan),6);
inert_moon_sv_an = zeros(length(tspan),6);
for i = 1:length(tspan)
    Kepler_an(i,:) = Kepl2Del(mu, Delaunay_an(i,:)', true )';
    inert_moon_sv_nondim_an(i,:) = Kepl2Dec(mu, Kepler_an(i,:)); %this is equatorial inertial cs
    %giving back the dimensions
    inert_moon_sv_an(i,1:3) = inert_moon_sv_nondim_an(i,1:3).*runit;
    inert_moon_sv_an(i,4:6) = inert_moon_sv_nondim_an(i,4:6).*vunit;
end

%Deleting leaps in Kepler elements
for i = 1:length(tspan)
    for j = 3:5
        if (Kepler_an(i,j) > 10*pi/180) && (abs(Kepler_an(i,j) - 2*pi) < 10*pi/180)
            Kepler_an(i,j) = Kepler_an(i,j) - 2*pi;
        end
        if (Kepler(i,j) > 10*pi/180) && (abs(Kepler(i,j) - 2*pi) < 10*pi/180)
            Kepler(i,j) = Kepler(i,j) - 2*pi;
        end 
    end
end

%% Plots
t_sample = downsample(t, 1);

figure;
tiledlayout(2,3)

nexttile
plot(t_sample/86400, Kepler(:,1));
hold on
plot(t_sample/86400, Kepler_an(:,1)*runit);
legend('numerical', 'analytical');
title('a(t), Km');

nexttile
plot(t_sample/86400, Kepler(:,2));
hold on
plot(t_sample/86400, Kepler_an(:,2));
legend('numerical', 'analytical');
title('e(t)');

nexttile
plot(t_sample/86400, rad2deg(Kepler(:,3)));
hold on
plot(t_sample/86400, rad2deg(Kepler_an(:,3)));
legend('numerical', 'analytical');
title('i(t), deg');

nexttile
plot(t_sample/86400, rad2deg(Kepler(:,4)));
hold on
plot(t_sample/86400, rad2deg(Kepler_an(:,4)));
legend('numerical', 'analytical');
title('omega(t), deg');

nexttile
plot(t_sample/86400, rad2deg(Kepler(:,5)));
hold on
plot(t_sample/86400, rad2deg(Kepler_an(:,5)));
legend('numerical', 'analytical');
title('RAAN(t), deg');

nexttile
plot(t_sample/86400, rad2deg(Kepler(:,6)));
hold on
plot(t_sample/86400, rad2deg(Kepler_an(:,6)));
legend('numerical', 'analytical');
title('M(t), deg');

Error_r = zeros(length(tspan),1);
Error_v = zeros(length(tspan),1);

for i = 1:length(tspan)
    Error_r(i) = norm(inert_moon_sv_an(i,1:3) - inert_moon_sv(i,1:3));
    Error_v(i) = norm(inert_moon_sv_an(i,4:6) - inert_moon_sv(i,4:6));
end

%decreasing amount of data for plots
Error_r = downsample(Error_r, 1);
Error_v = downsample(Error_v, 1);

figure;
tiledlayout(2,1)

nexttile
plot(t_sample/86400, Error_r);
title('\Deltar, km');
xlabel('Time [days]','FontWeight','bold');

nexttile
plot(t_sample/86400, Error_v);
title('\Deltav, km/s');
xlabel('Time [days]','FontWeight','bold');

Error_Kepler = zeros(length(tspan),6);
Error_omega_MA = zeros(length(tspan),1); %error omega+MA
for i = 1:length(tspan)

    Error_Kepler(i,1) = Kepler_an(i,1)*runit - Kepler(i,1);
    Error_Kepler(i,2:6) = Kepler_an(i,2:6) - Kepler(i,2:6);
    Error_omega_MA(i,1) = Kepler_an(i,4) + Kepler_an(i,6) - (Kepler(i,4) + Kepler(i,6));
end

%Deleting leaps in Kepler element's errors
for i = 1:length(tspan)
    if (Error_Kepler(i,6) > 0*pi/180) && (abs(Error_Kepler(i,6) - 2*pi) < 30*pi/180)
        Error_Kepler(i,6) = Error_Kepler(i,6) - 2*pi;
    end
    if (Error_Kepler(i,6) < 0*pi/180) && (abs(Error_Kepler(i,6) + 2*pi) < 180*pi/180)
        Error_Kepler(i,6) = Error_Kepler(i,6) + 2*pi;
    end
    if (Error_omega_MA(i,1) > 0*pi/180) && (abs(Error_omega_MA(i,1) - 2*pi) < 20*pi/180)
        Error_omega_MA(i,1) = Error_omega_MA(i,1) - 2*pi;
    end
    if (Error_omega_MA(i,1) < 0*pi/180) && (abs(Error_omega_MA(i,1) + 2*pi) < 40*pi/180)
        Error_omega_MA(i,1) = Error_omega_MA(i,1) + 2*pi;
    end
end

%decreasing amount of data for plots
Error_Kepler = downsample(Error_Kepler, 1);
Error_omega_MA = downsample(Error_omega_MA, 1);

figure;
tiledlayout(2,3)
nexttile
plot(t_sample/86400, Error_Kepler(:,1));
title('\Deltaa(t), Km');
nexttile
plot(t_sample/86400, Error_Kepler(:,2));
title('\Deltae(t)');
nexttile
plot(t_sample/86400, rad2deg(Error_Kepler(:,3)));
title('\Deltai(t), deg');
nexttile
plot(t_sample/86400, rad2deg(Error_Kepler(:,4)));
title('\Deltaomega(t), deg');
nexttile
plot(t_sample/86400, rad2deg(Error_Kepler(:,5)));
title('\DeltaRAAN(t), deg');
nexttile
plot(t_sample/86400, rad2deg(Error_Kepler(:,6)));
title('\DeltaM(t), deg');

% figure;
% plot(t_sample/86400, rad2deg(Error_omega_MA(:,1)));
% title('\Deltaomega+MA, deg');

figure;
tiledlayout(2,3)

nexttile
plot(t_sample/86400, Delaunay_mean(:,1));
hold on
plot(t_sample/86400, Delaunay_mean_num(:,1));
title('l_{mean}');
legend('mean analytical', 'mean numerical');

nexttile
plot(t_sample/86400, Delaunay_mean(:,2));
hold on
plot(t_sample/86400, Delaunay_mean_num(:,2));
title('g_{mean}');
legend('mean analytical', 'mean numerical');

nexttile
plot(t_sample/86400,Delaunay_mean(:,3));
hold on
plot(t_sample/86400, Delaunay_mean_num(:,3));
title('h_{mean}');
legend('mean analytical', 'mean numerical');

nexttile
plot(t_sample/86400, Delaunay_mean(:,4));
hold on
plot(t_sample/86400, Delaunay_mean_num(:,4));
title('L_{mean}');
legend('mean analytical', 'mean numerical');

nexttile
plot(t_sample/86400, Delaunay_mean(:,5));
hold on
plot(t_sample/86400, Delaunay_mean_num(:,5));
title('G_{mean}');
legend('mean analytical', 'mean numerical');

nexttile
plot(t_sample/86400, Delaunay_mean(:,6));
hold on
plot(t_sample/86400, Delaunay_mean_num(:,6));
title('H_{mean}');
legend('mean analytical', 'mean numerical');

figure;
tiledlayout(2,3)

nexttile
plot(t_sample/86400, rad2deg(Delaunay_nondim(:,1)));
hold on
plot(t_sample/86400, rad2deg(Delaunay_an(:,1)));
legend('numerical', 'analytical');
title('l, deg');
xlabel('Time [days]','FontWeight','bold');

nexttile
plot(t_sample/86400, rad2deg(Delaunay_nondim(:,2)));
hold on
plot(t_sample/86400, rad2deg(Delaunay_an(:,2)));
legend('numerical', 'analytical');
title('g, deg');
xlabel('Time [days]','FontWeight','bold');

nexttile
plot(t_sample/86400, rad2deg(Delaunay_nondim(:,3)));
hold on
plot(t_sample/86400, rad2deg(Delaunay_an(:,3)));
legend('numerical', 'analytical');
title('h, deg');
xlabel('Time [days]','FontWeight','bold');

nexttile
plot(t_sample/86400, Delaunay_nondim(:,4));
hold on
plot(t_sample/86400, Delaunay_an(:,4));
legend('numerical', 'analytical');
title('L');
xlabel('Time [days]','FontWeight','bold');

nexttile
plot(t_sample/86400, Delaunay_nondim(:,5));
hold on
plot(t_sample/86400, Delaunay_an(:,5));
legend('numerical', 'analytical');
title('G');
xlabel('Time [days]','FontWeight','bold');

nexttile
plot(t_sample/86400, Delaunay_nondim(:,6));
hold on
plot(t_sample/86400, Delaunay_an(:,6));
legend('numerical', 'analytical');
title('H');
xlabel('Time [days]','FontWeight','bold');

Error_Delaunay_l = Delaunay_an(:,1)-Delaunay_nondim(:,1);
for i = 1:length(tspan)
    if (Error_Delaunay_l(i,1) < 0) && (abs(Error_Delaunay_l(i,1) + 2*pi) < 40*pi/180)
        Error_Delaunay_l(i,1) = Error_Delaunay_l(i,1) + 2*pi;
    end
end

figure;
tiledlayout(2,3)

nexttile
plot(t_sample/86400, rad2deg(Error_Delaunay_l));
title('\Deltal, deg');
xlabel('Time [days]','FontWeight','bold');

nexttile
plot(t_sample/86400, rad2deg(Delaunay_an(:,2)-Delaunay_nondim(:,2)));
title('\Deltag, deg');
xlabel('Time [days]','FontWeight','bold');

nexttile
plot(t_sample/86400, rad2deg(Delaunay_an(:,3)-Delaunay_nondim(:,3)));
title('\Deltah, deg');
xlabel('Time [days]','FontWeight','bold');

nexttile
plot(t_sample/86400, Delaunay_an(:,4)-Delaunay_nondim(:,4));
title('\DeltaL');
xlabel('Time [days]','FontWeight','bold');

nexttile
plot(t_sample/86400, Delaunay_an(:,5)-Delaunay_nondim(:,5));
title('\DeltaG');
xlabel('Time [days]','FontWeight','bold');

nexttile
plot(t_sample/86400, Delaunay_an(:,6)-Delaunay_nondim(:,6));
title('\DeltaH');
xlabel('Time [days]','FontWeight','bold');

figure;
[X_sph, Y_sph, Z_sph] = sphere(1000);
X_sph = 1737.4*X_sph;
Y_sph = 1737.4*Y_sph;
Z_sph = 1737.4*Z_sph;
plot3(X_sph, Y_sph, Z_sph);
hold on;
plot3(inert_moon_sv_an(:,1), inert_moon_sv_an(:,2),inert_moon_sv_an(:,3)); %MOON EQUATORIAL NON ROTATIONAL FRAME
plot3(inert_moon_sv(:,1), inert_moon_sv(:,2),inert_moon_sv(:,3)); %MOON EQUATORIAL NON ROTATIONAL FRAME
axis equal;