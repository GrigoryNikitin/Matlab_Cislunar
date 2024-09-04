%Here are the long-period terms being calculated for diff semi-major axis
%and different combinations of J coefficients
%% Initial Conditions and stuff
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

thirdbody_analytic_bool = false;

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

number_of_points = 3000;

%%Initial conditions
% Moon rotation rate
we = 2*pi/(27.3220*86400/tunit); %rad/sec nondimensionalized
M0_moon = deg2rad(0); %M0 of Moon on its orbit
%[~, r2] = calc_barycenter(0, M0_moon); %positions of the Earth and the Moon from the barycenter

%% Setting J_coeff to zeros
% Set whatever I need each time
%Jcoeff(3:5,1) = 0;
Jcoeff(7:end,1) = 0;

%% Moon centered ANALYTICAL integration w/ or w/o Earth (nondimensionalized)
num_of_a = 2;
Zonal_LP = zeros(number_of_points,6,num_of_a);
a0 = linspace(2000, 3000, num_of_a);
for a = 1:length(a0)
    
    disp(a);
    
    Kepler_an0 = [a0(a)/runit,0.03,30*pi/180,30*pi/180,30*pi/180,0*pi/180]';
    finaltime = 2*pi*sqrt(384748^3/398600.44)*130/tunit; %many years
    timespan = linspace(0,finaltime,number_of_points);
    tspan = timespan;
    t = tspan * tunit; %normal time
    
    Delaunay_an0 = Kepl2Del(mu, Kepler_an0, false);
    
    Kepler_Earth = calc_Earth_Kepler(0, M0_moon);
    
    Delaunay_mean0 = Osc2Mean_approx(Delaunay_an0, thirdbody_analytic_bool, Kepler_Earth, 0, 0, 0, Degree, Order, Clm, Slm, mu, re, we, tol, quadtol);
    
    
    %Delaunay_an0_check = Mean2Osc(Delaunay_mean0, 0, 0, 0, Degree, Order, Clm, Slm, mu, re, we, tol, quadtol);
    
    % Propagate Mean States
    Delaunay_mean = MeanProp(tspan, Delaunay_mean0, thirdbody_analytic_bool, mu, re, Jcoeff);
    
    % LP variations
    %calculating long-period variations of Delaunay elements
    for i = 1:length(tspan)
        % Zonal Long-period variations
        aaa = LPVarJ26(Delaunay_mean(i,:)', mu, re, Jcoeff);
        Zonal_LP(i,:,a) = LPVarJ26(Delaunay_mean(i,:)', mu, re, Jcoeff);
    end
end


%% Plots

%plots of LP variations
figure;
tiledlayout(2,2);

% Initialize a cell array to store plot handles for the legend
plotHandles = cell(9,1);

nexttile;
plot(t/86400/365, Zonal_LP(:,1,1)*180/pi, '-');
hold on;
for k = 2:num_of_a
    plot(t/86400/365, Zonal_LP(:,1,k)*180/pi, '-');
end
hold off;
xlabel('Time [years]','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('l long period','FontWeight','bold');
xlim([0 10]);
grid on;

nexttile;
plot(t/86400/365, Zonal_LP(:,2,1)*180/pi, '-');
hold on;
for k = 2:num_of_a
    plot(t/86400/365, Zonal_LP(:,2,k)*180/pi, '-');
end
hold off;
xlabel('Time [years]','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('g long period','FontWeight','bold');
xlim([0 10]);
grid on;

nexttile;
plot(t/86400/365, Zonal_LP(:,3,1)*180/pi, '-');
hold on;
for k = 2:num_of_a
    plot(t/86400/365, Zonal_LP(:,3,k)*180/pi, '-');
end
hold off;
xlabel('Time [years]','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('h long period','FontWeight','bold');
xlim([0 10]);
grid on;

% nexttile;
% plot(t/86400/365, Zonal_LP(:,4,1), '-');
% hold on;
% for k = 2:num_of_a
%     plot(t/86400/365, Zonal_LP(:,4,k)*180/pi, '-');
% end
% hold off;
% xlabel('Time [years]','FontWeight','bold');
% ylabel('[-]','FontWeight','bold');
% title('L long period (nondim)','FontWeight','bold');
% xlim([0 10]);
% grid on;

nexttile;
plot(t/86400/365, Zonal_LP(:,5,1), '-');
hold on;
for k = 2:num_of_a
    plot(t/86400/365, Zonal_LP(:,5,k)*180/pi, '-');
end
hold off;
xlabel('Time [years]','FontWeight','bold');
ylabel('[-]','FontWeight','bold');
title('G long period (nondim)','FontWeight','bold');
xlim([0 10]);
grid on;

% nexttile;
% plot(t/86400/365, Zonal_LP(:,6,1), '-');
% hold on;
% for k = 2:num_of_a
%     plot(t/86400/365, Zonal_LP(:,6,k)*180/pi, '-');
% end
% hold off;
% xlabel('Time [years]','FontWeight','bold');
% ylabel('[-]','FontWeight','bold');
% title('H long period (nondim)','FontWeight','bold');
% xlim([0 10]);
% grid on;

fig = gcf; % Get current figure handle
ax = findall(fig, 'Type', 'axes'); % Find all axes in the figure

for i = 1:length(ax)
    %ytickformat(ax(i), '%.1f'); % Format the y-ticks with 2 decimal places
end

legendLabels = arrayfun(@num2str, a0, 'UniformOutput', false);
lgd = legend([plotHandles{:}], legendLabels);

lgd.Orientation = 'vertical';
lgd.Location = 'eastoutside';
lgd.Position = [0.91, 0.5, 0.085, 0.4]; % [left, bottom, width, height]

%% PLots of the differences
clear;
clc;

number_of_points = 3000;
num_of_a = 2;
a0 = linspace(2000, 3000, num_of_a);
mu_Moon = 4902.8001224453001;
tunit = sqrt(1738.0^3/(mu_Moon));
finaltime = 2*pi*sqrt(384748^3/398600.44)*130/tunit; %many years
timespan = linspace(0,finaltime,number_of_points);
tspan = timespan;
t = tspan * tunit; %normal time

Zonal_LP1 = Zonal_LP;
Zonal_LP2 = Zonal_LP;

%plots of LP variations
figure;
tiledlayout(2,2);

% Initialize a cell array to store plot handles for the legend
plotHandles = cell(9,1);

nexttile;
plot(t/86400/365, Zonal_LP2(:,1,1)*180/pi-Zonal_LP1(:,1,1)*180/pi, '-');
hold on;
for k = 2:num_of_a
    plot(t/86400/365, Zonal_LP2(:,1,k)*180/pi-Zonal_LP1(:,1,k)*180/pi, '-');
end
hold off;
xlabel('Time [years]','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('\deltal long period','FontWeight','bold');
xlim([0 10]);
grid on;

nexttile;
plot(t/86400/365, Zonal_LP2(:,2,1)*180/pi-Zonal_LP1(:,2,1)*180/pi, '-');
hold on;
for k = 2:num_of_a
    plot(t/86400/365, Zonal_LP2(:,2,k)*180/pi-Zonal_LP1(:,2,k)*180/pi, '-');
end
hold off;
xlabel('Time [years]','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('\deltag long period','FontWeight','bold');
xlim([0 10]);
grid on;

nexttile;
plot(t/86400/365, Zonal_LP2(:,3,1)*180/pi-Zonal_LP1(:,3,1)*180/pi, '-');
hold on;
for k = 2:num_of_a
    plot(t/86400/365, Zonal_LP2(:,3,k)*180/pi-Zonal_LP1(:,3,k)*180/pi, '-');
end
hold off;
xlabel('Time [years]','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('\deltah long period','FontWeight','bold');
xlim([0 10]);
grid on;

nexttile;
plot(t/86400/365, Zonal_LP2(:,5,1)-Zonal_LP1(:,5,1), '-');
hold on;
for k = 2:num_of_a
    plot(t/86400/365, Zonal_LP2(:,5,k)*180/pi-Zonal_LP1(:,5,k)*180/pi, '-');
end
hold off;
xlabel('Time [years]','FontWeight','bold');
ylabel('[-]','FontWeight','bold');
title('\deltaG long period (nondim)','FontWeight','bold');
xlim([0 10]);
grid on;


fig = gcf; % Get current figure handle
ax = findall(fig, 'Type', 'axes'); % Find all axes in the figure

for i = 1:length(ax)
    %ytickformat(ax(i), '%.1f'); % Format the y-ticks with 2 decimal places
end

legendLabels = arrayfun(@num2str, a0, 'UniformOutput', false);
lgd = legend([plotHandles{:}], legendLabels);

lgd.Orientation = 'vertical';
lgd.Location = 'eastoutside';
lgd.Position = [0.91, 0.5, 0.085, 0.4]; % [left, bottom, width, height]

%% Local functions
%Successive approximation to find new Mean elements
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