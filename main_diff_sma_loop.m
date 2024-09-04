% here are big loops for different r and inc of orbits

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
thirdbody_analytic_bool = true;

% Units for nondimensionalization 
runit = 1738.0;
tunit = sqrt(1738.0^3/(mu_Moon));
vunit = runit/tunit;

%normilized r_Moon and mu_Moon
re = 1;
mu = 1;

% Chosen degree and order of Moon's gravitational field
Degree = 6;
Order = 6;

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


%% The beginning of the FOR loop

% for r0 = 1800:10000, ecc = 0.03, inc = ~ deg

%For r0 = 2000, step = 10-15 sec
%For r0 = 4000, step = 20-30 sec
%For r0 = 8000, step = 40 sec
%For r0 = 10000, step = 60 sec

iteration = 0;
%arrays for making plots
array_of_r = zeros(33,1);
LP_of_r = zeros(33,6);
Error_r_of_r = zeros(33,1);
Error_v_of_r = zeros(33,1);
Error_Kepler_of_r = zeros(33,6);

for r0 = 1800:200:8200
    
    iteration = iteration + 1 %#ok<NOPTS>
    
    array_of_r(iteration) = r0;
    
    sv0_equatorial_inert = Kepl2Dec(mu,[r0/runit,0.03,1*pi/180,0*pi/180,0*pi/180,0*pi/180]);
    
    % Moon centered integration with the Earth (nondimensionalized)
    
    finaltime = 2*pi*sqrt(384748^3/398600.44)/4/tunit; %nondimensionalized
    number_of_points = 5000;
    timespan = linspace(0,finaltime,number_of_points);
    
    % arrays for integration in Moon Equatorial inertial CS
    inert_moon_sv = zeros(number_of_points,6); % normal state-vector
    %Orbital Inertial CS
    orbital_moon_sv = zeros(number_of_points,6);
    % Barycentric CS for plots
    r1 = zeros(number_of_points,3);
    r2 = zeros(number_of_points,3);
    rot_r1 = zeros(number_of_points,3);
    rot_r2 = zeros(number_of_points,3);
    rot_sv = zeros(number_of_points,6);
    
    mooneq_sv_rotating = zeros(number_of_points,6); %rotating equatorial CS
    
    Kepler = zeros(number_of_points,6);
    Delaunay = zeros(number_of_points,6);
    inert_cm_sv_crds = zeros(number_of_points,3);
    inert_cm_sv_velocity = zeros(number_of_points,3);
    inert_cm_sv = zeros(number_of_points,6);
    radius = zeros(number_of_points,1);
    speed = zeros(number_of_points,1);
    Hamiltonian_bary = zeros(number_of_points,1);
    Hamiltonian_Moon = zeros(number_of_points,1);
    Hamiltonian_rot = zeros(number_of_points,1);
    omega = zeros(number_of_points,1);
    
    %ODE integration
    optn = odeset('RelTol',tol,'AbsTol',tol*1e-10);
    [tspan,inert_moon_sv_nondim] = ode113(@(tspan, inert_moon_sv_nondim) derivatives_Moon_cntrd_nondim(tspan, inert_moon_sv_nondim, M0_moon, Degree, Order, Clm, Slm),timespan,sv0_equatorial_inert,optn);
%ODE integration ^^^
    
    t = tspan * tunit; %normal time
    
    for i = 1:number_of_points
        
        inert_moon_sv(i, 1:3) = inert_moon_sv_nondim(i,1:3)*runit;
        inert_moon_sv(i, 4:6) = inert_moon_sv_nondim(i,4:6)*vunit;
        
        [r1(i,:), r2(i,:)] = calc_barycenter(t(i), M0_moon); %positions of the Earth and the Moon from the barycenter
        orbital_moon_sv(i,:) = Orbital2MoonEq_inert(inert_moon_sv(i,:), 0);
        %transfering to the barycentric coordinate system
        [r_EM, v_EM] = Moon_crds_from_Earth(t(i), M0_moon); %Moon velocity in the Earth centered orbital ref frame
        v_cm = norm(v_EM)*norm(r1(i,:))/norm(r_EM); %speed of the center of mass with respect to the Earth
        v2 = v_EM/norm(v_EM)*(norm(v_EM) - v_cm); %Moon velocity from barycenter
        inert_cm_sv_crds(i,:) = [r2(i,1) + orbital_moon_sv(i,1); r2(i,2) + orbital_moon_sv(i,2); r2(i,3) + orbital_moon_sv(i,3)];
        inert_cm_sv_velocity(i,:) = [v2(1) + orbital_moon_sv(i,4); v2(2) + orbital_moon_sv(i,5); v2(3) + orbital_moon_sv(i,6)];
        inert_cm_sv(i,:) = [inert_cm_sv_crds(i,:), inert_cm_sv_velocity(i,:)];
        %transfering to rotating ref frame
        [rot_r1(i,:), rot_r2(i,:), rot_sv(i,:)] = inert2rot(r1(i,:), r2(i,:), v2, inert_cm_sv(i,:));
        
        %transfering to the Moon equator coordinate system
        %rotating
        [mooneq_sv_rotating(i,:)] = MoonEq_inert2MoonEq_rotating(t(i), inert_moon_sv(i,:), 1);
        
        %calculating Kepler elements
        %Kepler(i,:) = Dec2Kepl(mu_Moon, mooneq_sv(i,:));
        Kepler(i,:) = Dec2Kepl(mu_Moon, inert_moon_sv(i,:)); %FROM MOON EQUATORIAL INERTIAL FRAME
        Delaunay(i,:) = Kepl2Del(mu_Moon, Kepler(i,:), false);
        
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
    
    %Deleting leaps in Kepler elements
    for i = 1:number_of_points
        for j = 3:5
            if (Kepler(i,j) > 10*pi/180) && (abs(Kepler(i,j) - 2*pi) < 10*pi/180)
                Kepler(i,j) = Kepler(i,j) - 2*pi;
            end
        end
    end
    
    % Moon centered ANALYTICAL integration w/o Earth (nondimensionalized)
    
    %Compute Absolute Initial Mean States
    
    % analytical Kepler and Delaunay elements at time = 0
    %GOING BACK TO NONDIMINSIONALIZED
    Kepler_an0 = [Kepler(1,1)/runit, Kepler(1,2:6)]';
    Delaunay_an0 = Kepl2Del(mu, Kepler_an0, false);
    
    Kepler_Earth = calc_Earth_Kepler(0, M0_moon);
    
    %approximation for mean elements
    Delaunay_mean00 = Osc2Mean(Delaunay_an0, thirdbody_analytic_bool, Kepler_Earth, 0, 0, 0, Degree, Order, Clm, Slm, mu, re, we, tol, quadtol);
    Delaunay_mean0 = Osc2Mean_approx(Delaunay_an0, thirdbody_analytic_bool, Kepler_Earth, 0, 0, 0, Degree, Order, Clm, Slm, mu, re, we, tol, quadtol);
    
    % Propagate Mean States
    %a new tspan just of first, last and n elements before the last from the original tspan!
    n = 700; %the number of point that is checked before the last one for the MAX of error
    tspan_new = zeros(2 + n, 1); %ONLY 2 + n ELEMENTS!
    r2_new = zeros(2 + n, 3);
    tspan_new(1) = tspan(1); 
    tspan_new(end) = tspan(end);
    r2_new(1,:) = r2(1,:);
    r2_new(end,:) = r2(end,:);
    for i = 1:n
        tspan_new(end-i) = tspan(end-i);
        r2_new(end-i,:) = r2(end-i,:);
    end
    
    Delaunay_mean = MeanProp(tspan_new, Delaunay_mean0, thirdbody_analytic_bool, mu, re, Jcoeff);
    
    Kepler_mean = zeros(length(tspan_new),6);
    for i = 1:length(tspan_new)
        Kepler_mean(i,:) = Kepl2Del( mu, Delaunay_mean(i,1:6)', true )'; %mean Kepler elements
    end
    
    % Mean to Osculating Transformation
    Delaunay_an = zeros(length(tspan_new),6);
    Kepler_Earth = zeros(length(tspan_new),6);
    
    for i = 1:length(tspan_new)
        Kepler_Earth(i,:) = calc_Earth_Kepler(tspan_new(i), M0_moon);
        Delaunay_an(i,:) = Mean2Osc(Delaunay_mean(i,:)', thirdbody_analytic_bool, Kepler_Earth(i,:), 0, tspan_new(i), 0, Degree, Order, Clm, Slm, mu, re, we, tol, quadtol)';
    end
    
    sv_analytical = zeros(length(tspan_new),6);
    Kepler_an = zeros(length(tspan_new),6);
    inert_moon_sv_nondim_an = zeros(length(tspan_new),6);
    inert_moon_sv_an = zeros(length(tspan_new),6);
    for i = 1:length(tspan_new)
        Kepler_an(i,:) = Kepl2Del( mu, Delaunay_an(i,:)', true )';
        inert_moon_sv_nondim_an(i,:) = Kepl2Dec(mu, Kepler_an(i,:)); %this is equatorial inertial cs
        %giving back the dimensions
        inert_moon_sv_an(i,1:3) = inert_moon_sv_nondim_an(i,1:3).*runit;
        inert_moon_sv_an(i,4:6) = inert_moon_sv_nondim_an(i,4:6).*vunit;
    end
    
    %Deleting leaps in Kepler elements
    for i = 1:length(tspan_new)
        for j = 3:5
            if (Kepler_an(i,j) > 10*pi/180) && (abs(Kepler_an(i,j) - 2*pi) < 10*pi/180)
                Kepler_an(i,j) = Kepler_an(i,j) - 2*pi;
            end
        end
    end
    
    % Long Period variations
    LP_of_r(iteration,:) = LPVarJ26(Delaunay_mean(end,:)', mu, re, Jcoeff);
    
    % Errors between two methods
    
    %Choosing MAX error between the last n elements in r,v and Kepler
    Error_r = zeros(length(tspan_new),1);
    Error_v = zeros(length(tspan_new),1);
    Error_Kepler = zeros(length(tspan_new),1);
    for i = 0:n
        %PROBLEM was HERE!!!
        Error_r(i+1) = norm(inert_moon_sv_an(end-i,1:3) - inert_moon_sv(end-i,1:3));
        Error_v(i+1) = norm(inert_moon_sv_an(end-i,4:6) - inert_moon_sv(end-i,4:6));
        Error_Kepler(i+1,1) = Kepler_an(end-i,1)*runit - Kepler(end-i,1);
        Error_Kepler(i+1,2:6) = Kepler_an(end-i,2:6) - Kepler(end-i,2:6);
        %Deleting leaps in MA error
        if (Error_Kepler(i+1,6) > 20*pi/180) && (abs(Error_Kepler(i+1,6) - 2*pi) < 20*pi/180)
            Error_Kepler(i+1,6) = Error_Kepler(i+1,6) - 2*pi;
        end
        %Deleting leaps in RAAN error
        if (Error_Kepler(i+1,5) > 20*pi/180) && (abs(Error_Kepler(i+1,5) - 2*pi) < 20*pi/180)
            Error_Kepler(i+1,5) = Error_Kepler(i+1,5) - 2*pi;
        end
    end

    if max(Error_r) > abs(min(Error_r))
        Error_r_of_r(iteration) = max(Error_r);
    else
        Error_r_of_r(iteration) = min(Error_r);
    end
    if max(Error_v) > abs(min(Error_v))
        Error_v_of_r(iteration) = max(Error_v);
    else
        Error_v_of_r(iteration) = min(Error_v);
    end
    if max(Error_Kepler(:,1)) > abs(min(Error_Kepler(:,1)))
        [~, index] = max(Error_Kepler(:,1));
        Error_Kepler_of_r(iteration,:) = Error_Kepler(index,:);
        
    else
        [~, index] = min(Error_Kepler(:,1));
        Error_Kepler_of_r(iteration,:) = Error_Kepler(index,:);
    end
    
end

save('errors_of_r_variables.mat','array_of_r', 'LP_of_r','Error_r_of_r', 'Error_v_of_r', 'Error_Kepler_of_r');
%%
% figure;
% 
% plot(tspan_new*tunit/86400, Error_r);
% title('\Deltar, km');
% xlabel('Time [days]','FontWeight','bold');


%% PLOTS
figure;
tiledlayout(3,2);
nexttile;
plot(array_of_r, LP_of_r(:,1)*180/pi, '-');
xlabel('a0, Km','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('l long period(r0)','FontWeight','bold');
grid on;

nexttile;
plot(array_of_r, LP_of_r(:,2)*180/pi, '-');
xlabel('a0, Km','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('g long period(r0)','FontWeight','bold');
grid on;

nexttile;
plot(array_of_r, LP_of_r(:,3)*180/pi, '-');
xlabel('a0, Km','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('h long period(r0)','FontWeight','bold');
grid on;

nexttile;
plot(array_of_r, LP_of_r(:,4), '-');
xlabel('a0, Km','FontWeight','bold');
ylabel('[-]','FontWeight','bold');
title('L long period (nondim)(r0)','FontWeight','bold');
grid on;

nexttile;
plot(array_of_r, LP_of_r(:,5), '-');
xlabel('a0, Km','FontWeight','bold');
ylabel('[-]','FontWeight','bold');
title('G long period (nondim)(r0)','FontWeight','bold');
grid on;

nexttile;
plot(array_of_r, LP_of_r(:,6), '-');
xlabel('a0, Km','FontWeight','bold');
ylabel('[-]','FontWeight','bold');
title('H long period (nondim)(r0)','FontWeight','bold');
grid on;

figure;
tiledlayout(2,1)

nexttile
plot(array_of_r, Error_r_of_r);
title('\Deltar(r0), Km');

nexttile
plot(array_of_r, Error_v_of_r);
title('\Deltav(r0), Km/s');

figure;
tiledlayout(2,3)

nexttile
plot(array_of_r, Error_Kepler_of_r(:,1));
title('\Deltaa(r0), Km');

nexttile
plot(array_of_r, Error_Kepler_of_r(:,2));
title('\Deltae(r0)');

nexttile
plot(array_of_r, rad2deg(Error_Kepler_of_r(:,3)));
title('\Deltai(r0), deg');

nexttile
plot(array_of_r, rad2deg(Error_Kepler_of_r(:,4)));
title('\Deltaomega(r0), deg');

nexttile
plot(array_of_r, rad2deg(Error_Kepler_of_r(:,5)));
title('\DeltaRAAN(r0), deg');

nexttile
plot(array_of_r, rad2deg(Error_Kepler_of_r(:,6)));
title('\DeltaM(r0), deg');

%% Local Functions

function Em = MeanProp(tspan, Em0, thirdbody_included, mu, Re, Jcoeff)

m = length(tspan);

Emdot = SecularRates(Em0, mu, Re,  Jcoeff);

if thirdbody_included %if we want to include third body effects, then a new set of secular rates is calculated
    Emdot_tb = Secular_rates_third_body(0, Em0);
    Emdot = Emdot + Emdot_tb;
end

Em = zeros(m,6);

for ctr = 1:6
    Em(:,ctr) = Em0(ctr) + tspan*Emdot(ctr);
end

end

function Eo = Mean2Osc(Em, thirdbody_included, Kepler_Earth, GMST0, t, t0, AnICDegree, AnICOrder, Clm, Slm, mu, Re, we, tol, quadtol)

Jcoeff = -Clm(2:end,1);

% Zonal Long-period variations
DelZLP = LPVarJ26(Em, mu, Re, Jcoeff);

% Add long-period variations
Elp = Em + DelZLP;

% Zonal short-period variations
DelZSP = SPVarJ26(Elp, mu, Re, Jcoeff);    

% Add short-period variations
Ezsp = Elp + DelZSP;

% tesseral Short-period variations
GMST = GMST0 + we*(t - t0);
DelTsp = SPMDTess(Ezsp,GMST,we,mu,Re,Clm,Slm,AnICDegree,AnICOrder,tol,quadtol, false);

% Add short-period variations due to tesserals
Eo = Ezsp + DelTsp;

%adding short and long period terms from a 3rd body effects if needed
if thirdbody_included
    Eo = Mean2Osc_thirdbody_new_one(Eo, Kepler_Earth);
end

end

function Em = Osc2Mean(Eo, thirdbody_included, Kepler_Earth, GMST0, t, t0, AnICDegree, AnICOrder, Clm, Slm, mu, Re, we, tol, quadtol)

%removing short and long period terms from a 3rd body effects if needed
if thirdbody_included
    Eo = Osc2Mean_thirdbody_new_one(Eo, Kepler_Earth, 0);
end

Jcoeff = -Clm(2:end,1);

% tesseral Short-period variations
GMST = GMST0 + we*(t - t0);
DelTsp = SPMDTess(Eo,GMST,we,mu,Re,Clm,Slm,AnICDegree,AnICOrder,tol,quadtol, false);

% subtract short-period variations due to tesserals
Ezsp = Eo - DelTsp;

% Zonal short-period variations
DelZSP = SPVarJ26(Ezsp, mu, Re, Jcoeff); 

% subtract short-period variations
Ezlp = Ezsp - DelZSP;

% Zonal Long-period variations
DelZLP = LPVarJ26(Ezlp, mu, Re, Jcoeff);

% Add long-period variations
Em = Ezlp - DelZLP;

end

%Successive approximation to find new Mean elements
%Here I do 3 approximations in a row
function Em = Osc2Mean_approx(Eo, thirdbody_included, Kepler_Earth, GMST0, t, t0, AnICDegree, AnICOrder, Clm, Slm, mu, Re, we, tol, quadtol)

%removing short and long period terms from a 3rd body effects if needed
if thirdbody_included
    Eo = Osc2Mean_thirdbody_new_one(Eo, Kepler_Earth, 1);
end

Jcoeff = -Clm(2:end,1);
GMST = GMST0 + we*(t - t0);

if AnICOrder > 0 % if there are Tesseral harmonics
    error = 1000;
    i = 0;
    % Substracting tesseral short period
    X_sp_i = Eo - SPMDTess(Eo,GMST,we,mu,Re,Clm,Slm,AnICDegree,AnICOrder,tol,quadtol, false);
    while error > 10^-12
        i = i + 1;
        
        X_sp_i_ = Eo - SPMDTess(X_sp_i,GMST,we,mu,Re,Clm,Slm,AnICDegree,AnICOrder,tol,quadtol, false);
        
        error = norm(X_sp_i_ - X_sp_i); %calculating error
        X_sp_i = X_sp_i_; %ith iteration equals (i+1)th
    end
    
    error = 1000;
    k = 0;
    % Substracting zonal short-period
    X_sp_k = X_sp_i_ - SPVarJ26(X_sp_i_, mu, Re, Jcoeff);
    while error > 10^-12
        k = k + 1;
        
        X_sp_k_ = X_sp_i_ - SPVarJ26(X_sp_k, mu, Re, Jcoeff);
        
        error = norm(X_sp_k_ - X_sp_k); %calculating error
        X_sp_k = X_sp_k_; %kth iteration equals (k+1)th
    end
else % if no Tesserals
    error = 1000;
    k = 0;
    % Substracting zonal short-period
    X_sp_k = Eo - SPVarJ26(Eo, mu, Re, Jcoeff);
    while error > 10^-12
        k = k + 1;
        
        % ONLY IF TESSERALS ARE 0
        X_sp_k_ = Eo - SPVarJ26(X_sp_k, mu, Re, Jcoeff);
        
        error = norm(X_sp_k_ - X_sp_k); %calculating error
        X_sp_k = X_sp_k_; %kth iteration equals (k+1)th
    end
end

error = 1000;
j = 0;
% Substracting zonal Long-period
X_lp_k = X_sp_k_ - LPVarJ26(X_sp_k_, mu, Re, Jcoeff);
while error > 10^-12
    j = j + 1;
    
    X_lp_k_ = X_sp_k_ - LPVarJ26(X_lp_k, mu, Re, Jcoeff);
    
    error = norm(X_lp_k_ - X_lp_k); %calculating error
    X_lp_k = X_lp_k_; %kth iteration equals (k+1)th
end

Em = X_lp_k_;
end

function DelXm = SPMDTess(Xm,GMST,we,mu,Re,Clm,Slm,MaxDegree,MaxOrder,tol,quadtol, QuadTesseralsOn) 

    % current Greenwich mean sidereal time
    theta = GMST;
    
    % Specified order must be smaller than degree
    MaxOrder = min([MaxDegree, MaxOrder]);
    
    % Delaunay elements
    l = Xm(1);
    g = mod(Xm(2),2*pi);
    h = mod(Xm(3),2*pi);
    L = Xm(4);
    G = Xm(5);
    H = Xm(6);
 
    a = L^2/mu;
    eta = G/L; 
    e = sqrt(1 - eta^2);
    i = acos(H/G);
    
    [~, f] = KeplerEqSolver(l, e, tol);
    
    f0 = 0;
    
    DelXm = zeros(6,1);
    
    % corrections for each m-Daily, sectorial and tesseral
    for n = 2:1:MaxDegree
        for m = 1:1:min(MaxOrder,n)
            
            % C and S coefficients
            C = Clm(n + 1,m + 1);
            S = Slm(n + 1,m + 1);
            C20 = Clm(3,1);
    
            % short-period and m-daily variations for each tesseral
            %coe = [a,e,i,h,g,f]';
            hr = h - theta;
            coer = [a,e,i,hr,g,f]';
            
            % GF Partials w.r.t. classical elements in the ECF frame
            [~, Wx, ~, ~,~,~,~] = TesseralQGFr(coer,f0,n,m,C20,C,S,mu,we,Re,quadtol,false,QuadTesseralsOn);
            Wx = Wx*C20^2/2;

            % l
            Delnm(1,1) = 2*L/mu*Wx(1) + eta^2/(e*L)*Wx(2);
            % g
            Delnm(2,1) = -eta^2/(e*G)*Wx(2) + cos(i)/(G*sin(i))*Wx(3);
            % h
            Delnm(3,1) = -1/(G*sin(i))*Wx(3);
            % L, G, H
            Delnm(4,1) = -Wx(6);
            Delnm(5,1) = -Wx(5);
            Delnm(6,1) = -Wx(4);

            DelXm = DelXm + Delnm;

        end
    end
end


function Xdot = SPAvgDiffEq(t, X, mu, Re, Jcoeff)

L = X(4);
G = X(5);
H = X(6);
g = X(2);

R__e = Re;

J(2) = Jcoeff(2);
J(3) = Jcoeff(3);
J(4) = Jcoeff(4);
J(5) = Jcoeff(5);
J(6) = Jcoeff(6);

ldot = (0.3675e4 / 0.2048e4 * (G + H) * (G ^ 4 - 18 * H ^ 2 * G ^ 2 + 33 * H ^ 4) * ((G ^ 4) - 0.20e2 / 0.7e1 * (G ^ 2) * (L ^ 2) + 0.9e1 / 0.7e1 * (L ^ 4)) * mu ^ 8 * (G - H) / (L ^ 8) / (G ^ 17) * cos((2 * g)) - 0.2205e4 / 0.4096e4 * (G + L) * ((G + H) ^ 2) * (G ^ 2 - 11 * H ^ 2) * mu ^ 8 * ((G - H) ^ 2) * ((G ^ 2) - 0.3e1 / 0.7e1 * (L ^ 2)) * (G - L) / (L ^ 8) / (G ^ 17) * cos((4 * g)) - 0.2625e4 / 0.2048e4 * ((G ^ 4) - 0.10e2 / 0.3e1 * (G ^ 2) * (L ^ 2) + 0.9e1 / 0.5e1 * (L ^ 4)) * mu ^ 8 * ((G ^ 6) - (21 * G ^ 4 * H ^ 2) + (63 * G ^ 2 * H ^ 4) - 0.231e3 / 0.5e1 * (H ^ 6)) / (L ^ 8) / (G ^ 17)) * J(6) * R__e ^ 6 + (-0.135e3 / 0.64e2 * sqrt((G ^ 2 - H ^ 2)) * (G ^ 4 - 14 * H ^ 2 * G ^ 2 + 21 * H ^ 4) * mu ^ 7 * ((G ^ 4) - 0.43e2 / 0.18e2 * (G ^ 2) * (L ^ 2) + 0.7e1 / 0.6e1 * (L ^ 4)) * ((-G ^ 2 + L ^ 2) ^ (-0.1e1 / 0.2e1)) / (G ^ 14) / (L ^ 7) * sin(g) + 0.105e3 / 0.128e3 * sqrt((G ^ 2 - H ^ 2)) * (G + L) * (G + H) * (G - 3 * H) * ((G ^ 2) - (L ^ 2) / 0.2e1) * mu ^ 7 * (G - H) * ((-G ^ 2 + L ^ 2) ^ (-0.1e1 / 0.2e1)) * (G + 3 * H) * (G - L) / (G ^ 14) / (L ^ 7) * sin((3 * g))) * J(5) * R__e ^ 5 + (-0.75e2 / 0.64e2 * (G + H) * mu ^ 6 * (G - H) * ((G ^ 2) - 0.3e1 / 0.5e1 * (L ^ 2)) * (G ^ 2 - 7 * H ^ 2) / (G ^ 11) / (L ^ 6) * cos((2 * g)) + 0.135e3 / 0.128e3 * (G + L) * ((G ^ 4) - (10 * H ^ 2 * G ^ 2) + 0.35e2 / 0.3e1 * (H ^ 4)) * mu ^ 6 * (G - L) / (G ^ 11) / (L ^ 6)) * J(4) * R__e ^ 4 + 0.3e1 / 0.2e1 * J(3) * ((-G ^ 2 + L ^ 2) ^ (-0.1e1 / 0.2e1)) * mu ^ 5 * ((G ^ 2) - 0.3e1 / 0.4e1 * (L ^ 2)) * sqrt((G ^ 2 - H ^ 2)) * (G ^ 2 - 5 * H ^ 2) / (G ^ 8) / (L ^ 5) * sin(g) * R__e ^ 3 - 0.3e1 / 0.4e1 * J(2) * (G ^ 2 - 3 * H ^ 2) / (L ^ 4) / (G ^ 5) * R__e ^ 2 * mu ^ 4;

gdot = (0.525e3 / 0.2048e4 * mu ^ 8 * (7 * G ^ 10 - 171 * G ^ 8 * H ^ 2 - 36 * G ^ 8 * L ^ 2 + 561 * G ^ 6 * H ^ 4 + 836 * G ^ 6 * H ^ 2 * L ^ 2 + 33 * G ^ 6 * L ^ 4 - 429 * G ^ 4 * H ^ 6 - 2652 * G ^ 4 * H ^ 4 * L ^ 2 - 741 * G ^ 4 * H ^ 2 * L ^ 4 + 1980 * G ^ 2 * H ^ 6 * L ^ 2 + 2295 * G ^ 2 * H ^ 4 * L ^ 4 - 1683 * H ^ 6 * L ^ 4) / (L ^ 7) / (G ^ 18) * cos((2 * g)) - 0.315e3 / 0.4096e4 * mu ^ 8 * (G - L) * (G + L) * (G - H) * (G + H) * (7 * G ^ 6 - 110 * G ^ 4 * H ^ 2 - 11 * G ^ 4 * L ^ 2 + 143 * G ^ 2 * H ^ 4 + 158 * G ^ 2 * H ^ 2 * L ^ 2 - 187 * H ^ 4 * L ^ 2) / (L ^ 7) / (G ^ 18) * cos((4 * g)) - 0.2625e4 / 0.2048e4 * ((G ^ 10) + ((-27 * H ^ 2 - 6 * L ^ 2) * G ^ 8) + ((99 * H ^ 4) + (154 * H ^ 2 * L ^ 2) + 0.33e2 / 0.5e1 * (L ^ 4)) * (G ^ 6) + (-0.429e3 / 0.5e1 * (H ^ 6) - (546 * H ^ 4 * L ^ 2) - 0.819e3 / 0.5e1 * (H ^ 2) * (L ^ 4)) * (G ^ 4) + ((462 * H ^ 6 * L ^ 2 + 567 * H ^ 4 * L ^ 4) * G ^ 2) - 0.11781e5 / 0.25e2 * (H ^ 6) * (L ^ 4)) * mu ^ 8 / (L ^ 7) / (G ^ 18)) * J(6) * R__e ^ 6 + (-0.15e2 / 0.128e3 / (G ^ 15) / (L ^ 6) * mu ^ 7 * (18 * G ^ 10 - 357 * G ^ 8 * H ^ 2 - 77 * G ^ 8 * L ^ 2 + 1008 * G ^ 6 * H ^ 4 + 1445 * G ^ 6 * H ^ 2 * L ^ 2 + 63 * G ^ 6 * L ^ 4 - 693 * G ^ 4 * H ^ 6 - 3955 * G ^ 4 * H ^ 4 * L ^ 2 - 1148 * G ^ 4 * H ^ 2 * L ^ 4 + 2667 * G ^ 2 * H ^ 6 * L ^ 2 + 3087 * G ^ 2 * H ^ 4 * L ^ 4 - 2058 * H ^ 6 * L ^ 4) * ((G ^ 2 - H ^ 2) ^ (-0.1e1 / 0.2e1)) * ((-G ^ 2 + L ^ 2) ^ (-0.1e1 / 0.2e1)) * sin(g) + 0.105e3 / 0.256e3 * ((G ^ 2 - H ^ 2) ^ (-0.1e1 / 0.2e1)) * (G + L) * (2 * G ^ 6 - 27 * G ^ 4 * H ^ 2 - 3 * G ^ 4 * L ^ 2 + 33 * G ^ 2 * H ^ 4 + 37 * G ^ 2 * H ^ 2 * L ^ 2 - 42 * H ^ 4 * L ^ 2) * (G + H) * mu ^ 7 * (G - H) * (G - L) * ((-G ^ 2 + L ^ 2) ^ (-0.1e1 / 0.2e1)) / (G ^ 15) / (L ^ 6) * sin((3 * g))) * J(5) * R__e ^ 5 + (-0.15e2 / 0.64e2 / (G ^ 12) / (L ^ 5) * mu ^ 6 * (5 * G ^ 6 - 56 * G ^ 4 * H ^ 2 - 7 * G ^ 4 * L ^ 2 + 63 * G ^ 2 * H ^ 4 + 72 * G ^ 2 * H ^ 2 * L ^ 2 - 77 * H ^ 4 * L ^ 2) * cos((2 * g)) + 0.15e2 / 0.128e3 / (G ^ 12) / (L ^ 5) * mu ^ 6 * (9 * G ^ 6 - 126 * G ^ 4 * H ^ 2 - 21 * G ^ 4 * L ^ 2 + 189 * G ^ 2 * H ^ 4 + 270 * G ^ 2 * H ^ 2 * L ^ 2 - 385 * H ^ 4 * L ^ 2)) * J(4) * R__e ^ 4 + 0.3e1 / 0.8e1 / (G ^ 9) / (L ^ 4) * mu ^ 5 * (4 * G ^ 6 - 35 * G ^ 4 * H ^ 2 - 5 * G ^ 4 * L ^ 2 + 35 * G ^ 2 * H ^ 4 + 41 * G ^ 2 * H ^ 2 * L ^ 2 - 40 * H ^ 4 * L ^ 2) * ((G ^ 2 - H ^ 2) ^ (-0.1e1 / 0.2e1)) * ((-G ^ 2 + L ^ 2) ^ (-0.1e1 / 0.2e1)) * sin(g) * J(3) * R__e ^ 3 - 0.3e1 / 0.4e1 / (G ^ 6) / (L ^ 3) * (G ^ 2 - 5 * H ^ 2) * mu ^ 4 * J(2) * R__e ^ 2;

hdot = (0.9975e4 / 0.1024e4 * (G ^ 2 - 3 * L ^ 2) * (G + L) * H * mu ^ 8 * ((G ^ 4) - 0.102e3 / 0.19e2 * (H ^ 2) * (G ^ 2) + 0.99e2 / 0.19e2 * (H ^ 4)) * (G - L) / (L ^ 7) / (G ^ 17) * cos((2 * g)) - 0.4095e4 / 0.2048e4 * ((G + L) ^ 2) * (G + H) * H * ((G ^ 2) - 0.33e2 / 0.13e2 * (H ^ 2)) * mu ^ 8 * (G - H) * ((G - L) ^ 2) / (L ^ 7) / (G ^ 17) * cos((4 * g)) - 0.105e3 / 0.1024e4 * (15 * G ^ 4 - 70 * G ^ 2 * L ^ 2 + 63 * L ^ 4) * H * mu ^ 8 * (5 * G ^ 4 - 30 * H ^ 2 * G ^ 2 + 33 * H ^ 4) / (L ^ 7) / (G ^ 17)) * J(6) * R__e ^ 6 + (0.15e2 / 0.128e3 * ((G ^ 2 - H ^ 2) ^ (-0.1e1 / 0.2e1)) * (3 * G ^ 2 - 7 * L ^ 2) * (29 * G ^ 4 - 126 * H ^ 2 * G ^ 2 + 105 * H ^ 4) * H * mu ^ 7 * sqrt((-G ^ 2 + L ^ 2)) / (G ^ 14) / (L ^ 6) * sin(g) - 0.105e3 / 0.256e3 * ((G ^ 2 - H ^ 2) ^ (-0.1e1 / 0.2e1)) * (G + L) * (G + H) * H * mu ^ 7 * (G - H) * sqrt((-G ^ 2 + L ^ 2)) * (G - L) * (7 * G ^ 2 - 15 * H ^ 2) / (G ^ 14) / (L ^ 6) * sin((3 * g))) * J(5) * R__e ^ 5 + (-0.15e2 / 0.4e1 * (G + L) * H * mu ^ 6 * ((G ^ 2) - 0.7e1 / 0.4e1 * (H ^ 2)) * (G - L) / (G ^ 11) / (L ^ 5) * cos((2 * g)) + 0.15e2 / 0.32e2 * (3 * G ^ 2 - 5 * L ^ 2) * H * mu ^ 6 * (3 * G ^ 2 - 7 * H ^ 2) / (L ^ 5) / (G ^ 11)) * J(4) * R__e ^ 4 - 0.33e2 / 0.8e1 * ((G ^ 2 - H ^ 2) ^ (-0.1e1 / 0.2e1)) * H * mu ^ 5 * ((G ^ 2) - 0.15e2 / 0.11e2 * (H ^ 2)) * sqrt((-G ^ 2 + L ^ 2)) / (G ^ 8) / (L ^ 4) * sin(g) * J(3) * R__e ^ 3 - 0.3e1 / 0.2e1 * J(2) * H / (L ^ 3) / (G ^ 5) * R__e ^ 2 * mu ^ 4;

Ldot = 0; 

Gdot = (-0.525e3 / 0.1024e4 * mu ^ 8 * (G - L) * (G + L) * (G ^ 2 - 3 * L ^ 2) * (G - H) * (G + H) * (G ^ 4 - 18 * H ^ 2 * G ^ 2 + 33 * H ^ 4) / (L ^ 7) / (G ^ 17) * sin((2 * g)) + 0.315e3 / 0.1024e4 * mu ^ 8 * ((G - L) ^ 2) * ((G + L) ^ 2) * (G ^ 2 - 11 * H ^ 2) * ((G - H) ^ 2) * ((G + H) ^ 2) / (L ^ 7) / (G ^ 17) * sin((4 * g))) * J(6) * R__e ^ 6 + (0.45e2 / 0.128e3 * sqrt((G ^ 2 - H ^ 2)) * ((G ^ 2) - 0.7e1 / 0.3e1 * (L ^ 2)) * (G ^ 4 - 14 * H ^ 2 * G ^ 2 + 21 * H ^ 4) * mu ^ 7 * sqrt((-G ^ 2 + L ^ 2)) / (G ^ 14) / (L ^ 6) * cos(g) - 0.105e3 / 0.256e3 * sqrt((-G ^ 2 + L ^ 2)) * sqrt((G ^ 2 - H ^ 2)) * mu ^ 7 * (G - L) * (G + L) * (G - H) * (G - 3 * H) * (G + 3 * H) * (G + H) / (G ^ 14) / (L ^ 6) * cos((3 * g))) * J(5) * R__e ^ 5 + 0.15e2 / 0.32e2 * mu ^ 6 * (G - L) * (G + L) * (G - H) * (G + H) * (G ^ 2 - 7 * H ^ 2) / (G ^ 11) / (L ^ 5) * sin((2 * g)) * J(4) * R__e ^ 4 - 0.3e1 / 0.8e1 / (G ^ 8) / (L ^ 4) * sqrt((-G ^ 2 + L ^ 2)) * sqrt((G ^ 2 - H ^ 2)) * mu ^ 5 * (G ^ 2 - 5 * H ^ 2) * cos(g) * J(3) * R__e ^ 3;

Hdot = 0;

Xdot = [mu^2/L^3 + ldot, gdot, hdot, Ldot, Gdot, Hdot]';

end