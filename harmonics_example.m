% trying to make Moon harmonics work here

%% Initial conditions
clear;
clc;
Mass_Earth = 5.9722*10^24;
Mass_Moon = 7.3477*10^22;
G = 6.674*10^-20;

mu_Moon = 4902.8001224453001;

% ODE tolerances
tol = 1e-13;
quadtol = 1e-6;
optn = odeset('RelTol',tol,'AbsTol',tol*1e-3);

% Units for nondimensionalization 
lunit = 1738.0;
tunit = sqrt(1738.0^3/(mu_Moon));
vunit = lunit/tunit;

%re = 1738.0;
re = 1;
%mu = mu_Moon;
mu = 1;

r0 = 2000;
%sv0 = [r0; 0; 0; 0; 0.5; sqrt(G*Mass_Moon/r0)];
%sv0 = Kepl2Dec(mu, [r0,0.09,50*pi/180,50*pi/180,50*pi/180,50*pi/180]);
%sv0 = Kepl2Dec(mu, [r0/lunit,0.03,50*pi/180,50*pi/180,50*pi/180,0*pi/180]);
sv0 = Kepl2Dec(mu, [r0/lunit,0.03,50*pi/180,50*pi/180,50*pi/180,0*pi/180]);
sv0_norm = [r0/lunit; 0; 0; 0; 0.5/vunit; sqrt(G*Mass_Moon/r0)/vunit];

orb_period = 2*pi*sqrt((r0/lunit)^3/mu);

% Chosen degree and order of Moon's gravitational field
Degree = 6;
Order = 0;

% Moon rotation rate
we = (360/27.3220)/180*pi/86400*tunit;


[Clm,Slm] = DenormCS('Moon_GRGM1200A.txt', 70);
Jcoeff = -Clm(2:end,1);

%% Numerical propagation

%tspan = linspace(0,2*pi*sqrt(384748^3/398600.44)/tunit,30000);
%for n periods only
tspan = linspace(0,10*orb_period,3000);

tspan_norm = linspace(0,2*pi*sqrt(384748^3/398600.44)/tunit,30000);

% non-normalized (normalized too)
[T,sv] = ode113(@(T, sv) GravAcc_original(T, sv, 0, 0,we, Degree, Order, mu,re,Clm,Slm),tspan,sv0,optn);
% normalized
%[T_norm,ORV_norm] = ode113(@(T_norm, ORV_norm) GravAcc(T_norm, ORV_norm, 0, 0,we, Degree, Order, 1,1,Clm,Slm),tspan_norm,sv0_norm, optn);

% numerical Kepler and Delaunay elements
Kepler_num = zeros(length(T), 6);
Delaunay_num = zeros(length(T), 6);
for i = 1:length(T)
    Kepler_num(i,:) = Dec2Kepl(mu, sv(i,:));
    Delaunay_num(i,:) = Kepl2Del(mu, Kepler_num(i,:), false);
end

%% Plots

% figure;
% [X_sph, Y_sph, Z_sph] = sphere(50);
% X_sph = 1737.4*X_sph;
% Y_sph = 1737.4*Y_sph;
% Z_sph = 1737.4*Z_sph;
% plot3(X_sph, Y_sph, Z_sph);
% hold on;
% plot3(ORV_norm(:,1)*lunit, ORV_norm(:,2)*lunit, ORV_norm(:,3)*lunit);
% axis equal;
% hold off;

% ПОКА ЗАКОММЕНТИЛ ЭТО ТОЖЕ
% figure;
% [X_sph, Y_sph, Z_sph] = sphere(50);
% X_sph = lunit*X_sph;
% Y_sph = lunit*Y_sph;
% Z_sph = lunit*Z_sph;
% plot3(X_sph, Y_sph, Z_sph);
% hold on;
% plot3(sv(:,1)*lunit, sv(:,2)*lunit, sv(:,3)*lunit);
% axis equal;
% hold off;

% difference = zeros(size(T));
% for i = 1:length(T)
%     difference(i) = norm(sv(i,1:3)) - lunit*norm(ORV_norm(i,1:3));
% end
% figure;
% plot(T/86400, difference, '-');
% title('r - r_n');
% ylabel('delta, km');
% xlabel('time, days');


%% Compute Absolute Initial Mean States

% analytical Kepler and Delaunay elements at time = 0
Kepler_an0 = Kepler_num(1,:)';
Delaunay_an0 = Delaunay_num(1,:)';

Delaunay_mean00 = Osc2Mean(Delaunay_an0, 0, 0, 0, Degree, Order, Clm, Slm, mu, re, we, tol, quadtol);
%approximation for mean elements
Delaunay_mean0 = Osc2Mean_approx(Delaunay_an0, 0, 0, 0, Degree, Order, Clm, Slm, mu, re, we, tol, quadtol);

%Delaunay_an0_check = Mean2Osc(Delaunay_mean0, 0, 0, 0, Degree, Order, Clm, Slm, mu, re, we, tol, quadtol);

%disp(['Mean IC Computation Done: ' num2str(toc)]);

%% Propagate Mean States

Delaunay_mean = MeanProp(tspan, Delaunay_mean0, mu, re, Jcoeff, tol);

Kepler_mean = zeros(length(tspan),6);
for i = 1:length(tspan)
    Kepler_mean(i,:) = Kepl2Del( mu, Delaunay_mean(i,1:6)', true )'; %mean Kepler elements
end


%disp(['Mean Propagation Done: ' num2str(toc)]);

%% Mean to Osculating Transformation

Delaunay_an = zeros(length(tspan),6);

for i = 1:length(tspan)

    Delaunay_an(i,:) = Mean2Osc(Delaunay_mean(i,:)', 0, tspan(i), 0, Degree, Order, Clm, Slm, mu, re, we, tol, quadtol)';

end

%examples!!!!!!!
a1 = Osc2Mean_approx_2(Delaunay_an0, 0, 0, 0, Degree, Order, Clm, Slm, mu, re, we, tol, quadtol)';
b1 = Mean2Osc(a1', 0, 0, 0, Degree, Order, Clm, Slm, mu, re, we, tol, quadtol);
%

sv_analytical = zeros(length(tspan),6);
Kepler_an = zeros(length(tspan),6);
for i = 1:length(tspan)
    Kepler_an(i,:) = Kepl2Del( mu, Delaunay_an(i,:)', true )';
    sv_analytical(i,:) = Kepl2Dec(mu, Kepler_an(i,:));
end

% calculating difference in position and velocity after Osc2Mean and back
% for t = 0 between num and analytical
fprintf('e = %5.3f \n', Kepler_num(1,2));
fprintf('delta ecc = %7.5f \n', Kepler_num(1,2) - Kepler_an(1,2));
delta_dist = (norm(sv(1,1:3)) - norm(sv_analytical(1,1:3)))*lunit;
fprintf('delta R = %5.3f Km \n', delta_dist);
delta_veloc = (norm(sv(1,4:6)) - norm(sv_analytical(1,4:6)))*vunit;
fprintf('delta V = %6.4f Km/s \n', delta_veloc);
%the same difference but for Delaunay elements
Delaunay_errors = Delaunay_num(1,:) - Delaunay_an(1,:);
fprintf('delta l = %6.4f Deg \n', Delaunay_errors(1)*180/pi);
fprintf('delta g = %6.4f Deg \n', Delaunay_errors(2)*180/pi);
fprintf('delta h = %6.4f Deg \n', Delaunay_errors(3)*180/pi);
fprintf('delta l+g = %6.4f Deg \n', (Delaunay_errors(1)+Delaunay_errors(2))*180/pi);
fprintf('delta L = %7.5f \n', Delaunay_errors(4));
fprintf('delta G = %7.5f \n', Delaunay_errors(5));
fprintf('delta H = %7.5f \n', Delaunay_errors(6));

%% Plot Results

%close all;

adiff = @(x,y) atan((tan(x)-tan(y))./(1+tan(x).*tan(y)));

% 433 Eros
figure;
[x,y,z] = ellipsoid(0,0,0, 1738.0, 1738.0, 1738.0,1000);
surf(x,y,z, 'EdgeColor','none','FaceColor',[0.6,0.6,0.6]);
axis equal;
hold on;

% orbit
plot3(sv_analytical(:,1)*lunit,sv_analytical(:,2)*lunit,sv_analytical(:,3)*lunit);
xlabel('x [km]','FontWeight','bold');
ylabel('y [km]','FontWeight','bold');
zlabel('z [km]','FontWeight','bold');

% saveas(gcf, 'C:\Users\bmahajan\data\Research\Journal Papers\2019\Non-Earth AST\plots\eros200_3d','fig');
% saveas(gcf, 'C:\Users\bmahajan\data\Research\Journal Papers\2019\Non-Earth AST\plots\eros200_3d','png');

% classical elements plots

figure;
%tplot = T/86400;
tplot = T/86400*tunit;
subplot(3,1,1);
plot(tplot, Kepler_an(:,1),'b','LineWidth',1);
hold on;
plot(tplot, Kepler_num(:,1),'g:','LineWidth',1.2);
plot(tplot, Kepler_mean(:,1),'r--','LineWidth',2);
hold off;
% grid on;
% xlabel('Time [days]','FontWeight','bold');
ylabel('SMA [km]','FontWeight','bold');
legend('osc','num','mean');
subplot(3,1,2);
plot(tplot, Kepler_an(:,2),'b','LineWidth',1);
hold on;
plot(tplot, Kepler_num(:,2),'g:','LineWidth',1.2);
plot(tplot, Kepler_mean(:,2),'r--','LineWidth',2);
hold off;
% grid on;
% xlabel('Time [days]','FontWeight','bold');
ylabel('ECC','FontWeight','bold');
subplot(3,1,3); 
plot(tplot, Kepler_an(:,3)*180/pi,'b','LineWidth',1);
hold on;
plot(tplot, Kepler_num(:,3)*180/pi,'g:','LineWidth',1.2);
plot(tplot, Kepler_mean(:,3)*180/pi,'r--','LineWidth',2);
hold off;
% grid on;
xlabel('Time [days]','FontWeight','bold');
ylabel('INC [deg]','FontWeight','bold');

% saveas(gcf, 'C:\Users\bmahajan\data\Research\Journal Papers\2019\Non-Earth AST\plots\eros200_oe1','fig');
% saveas(gcf, 'C:\Users\bmahajan\data\Research\Journal Papers\2019\Non-Earth AST\plots\eros200_oe1','png');


figure;
subplot(3,1,1);
plot(tplot, Kepler_an(:,4)*180/pi,'b','LineWidth',1);
hold on;
plot(tplot, Kepler_num(:,4)*180/pi,'g:','LineWidth',1.2);
plot(tplot, Kepler_mean(:,4)*180/pi,'r--','LineWidth',2);
legend('osc','num','mean');
hold off;
% grid on;
% xlabel('Time [days]','FontWeight','bold');
ylabel('AOP [deg]','FontWeight','bold');
subplot(3,1,2);
plot(tplot, Kepler_an(:,5)*180/pi,'b','LineWidth',1);
hold on;
plot(tplot, Kepler_num(:,5)*180/pi,'g:','LineWidth',1.2);
plot(tplot, Kepler_mean(:,5)*180/pi,'r--','LineWidth',2);
hold off;
% grid on;
% xlabel('Time [days]','FontWeight','bold');
ylabel('RAAN [deg]','FontWeight','bold');
subplot(3,1,3);
plot(tplot, Kepler_an(:,6)*180/pi,'b','LineWidth',1);
hold on;
plot(tplot, Kepler_num(:,6)*180/pi,'g:','LineWidth',1.2);
plot(tplot, Kepler_mean(:,6)*180/pi,'r--','LineWidth',2);
hold off;
% grid on;
xlabel('Time [days]','FontWeight','bold');
ylabel('MA [deg]','FontWeight','bold');


% errors in classical elements plots

figure;
tplot = T/86400*tunit;
subplot(3,1,1);
plot(tplot, (Kepler_an(:,1)-Kepler_num(:,1)),'b','LineWidth',1.2);
hold on;
plot(tplot, (Kepler_mean(:,1)-Kepler_num(:,1)),'r:','LineWidth',1.2);
hold off;
% grid on;
% xlabel('Time [days]','FontWeight','bold');
ylabel('SMA [km]','FontWeight','bold');
legend('osc','mean');
subplot(3,1,2);
plot(tplot, (Kepler_an(:,2)-Kepler_num(:,2)),'b','LineWidth',1.2);
hold on;
plot(tplot, (Kepler_mean(:,2)-Kepler_num(:,2)),'r:','LineWidth',1.2);
hold off;
% grid on;
% xlabel('Time [days]','FontWeight','bold');
ylabel('ECC','FontWeight','bold');
subplot(3,1,3); 
plot(tplot, (Kepler_an(:,3)-Kepler_num(:,3))*180/pi,'b','LineWidth',1.2);
hold on;
plot(tplot, (Kepler_mean(:,3)-Kepler_num(:,3))*180/pi,'r:','LineWidth',1.2);
hold off;
% grid on;
xlabel('Time [days]','FontWeight','bold');
ylabel('INC [deg]','FontWeight','bold');


% saveas(gcf, 'C:\Users\bmahajan\data\Research\Journal Papers\2019\Non-Earth AST\plots\eros200_oeerr1','fig');
% saveas(gcf, 'C:\Users\bmahajan\data\Research\Journal Papers\2019\Non-Earth AST\plots\eros200_oeerr1','png');


figure;
subplot(3,1,1);
plot(tplot, ((Kepler_an(:,4)-Kepler_num(:,4)))*180/pi,'b','LineWidth',1.2);
hold on;
plot(tplot, ((Kepler_mean(:,4)-Kepler_num(:,4)))*180/pi-360,'r:','LineWidth',1.2);
legend('osc','mean');
hold off;
% grid on;
% xlabel('Time [days]','FontWeight','bold');
ylabel('AOP [deg]','FontWeight','bold');
subplot(3,1,2);
plot(tplot, (Kepler_an(:,5)-Kepler_num(:,5))*180/pi,'b','LineWidth',1.2);
hold on;
plot(tplot, (Kepler_mean(:,5)-Kepler_num(:,5))*180/pi,'r:','LineWidth',1.2);
hold off;
% grid on;
% xlabel('Time [days]','FontWeight','bold');
ylabel('RAAN [deg]','FontWeight','bold');
subplot(3,1,3);
plot(tplot, (Kepler_an(:,6)-Kepler_num(:,6))*180/pi-360,'b','LineWidth',1.2);
hold on;
plot(tplot, (Kepler_mean(:,6)-Kepler_num(:,6))*180/pi-360,'r--','LineWidth',1.2);
hold off;
% grid on;
xlabel('Time [days]','FontWeight','bold');
ylabel('MA [deg]','FontWeight','bold');

% saveas(gcf, 'C:\Users\bmahajan\data\Research\Journal Papers\2019\Non-Earth AST\plots\eros200_oeerr2','fig');
% saveas(gcf, 'C:\Users\bmahajan\data\Research\Journal Papers\2019\Non-Earth AST\plots\eros200_oeerr2','png');


% pose vel error
figure;
subplot(2,1,1);
plot(tplot, sqrt(sum(((sv(:,1:3)-sv_analytical(:,1:3))).^2,2))*lunit,'LineWidth',1.2);
% xlabel('Time [days]','FontWeight','bold');
ylabel('RSS Pos [km]','FontWeight','bold');
grid on;
subplot(2,1,2);
plot(tplot, sqrt(sum((sv(:,4:6)-sv_analytical(:,4:6)).^2,2))*vunit,'LineWidth',1.2);
xlabel('Time [days]','FontWeight','bold');
ylabel('RSS Vel [km/s]','FontWeight','bold');
grid on;

% short-period variations of Kepler elements (osculating - mean)
% ANALYTICAL
figure;
tiledlayout(3,2);
nexttile;
plot(tplot, (Kepler_an(:,1)-Kepler_mean(:,1))*lunit);
xlabel('Time [days]','FontWeight','bold');
ylabel('[km]','FontWeight','bold');
title('a osc - a mean [km]','FontWeight','bold');
grid on;

nexttile;
plot(tplot, Kepler_an(:,2)-Kepler_mean(:,2));
xlabel('Time [days]','FontWeight','bold');
title('e osc - e mean','FontWeight','bold');
grid on;

nexttile;
plot(tplot, (Kepler_an(:,3)-Kepler_mean(:,3))*180/pi);
xlabel('Time [days]','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('i osc - i mean [deg]','FontWeight','bold');
grid on;

nexttile;
plot(tplot, (Kepler_an(:,4)-Kepler_mean(:,4))*180/pi);
xlabel('Time [days]','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('omega osc - omega mean [deg]','FontWeight','bold');
grid on;

nexttile;
plot(tplot, (Kepler_an(:,5)-Kepler_mean(:,5))*180/pi);
xlabel('Time [days]','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('RAAN osc - RAAN mean [deg]','FontWeight','bold');
grid on;

nexttile;
plot(tplot, (Kepler_an(:,6)-Kepler_mean(:,6))*180/pi);
xlabel('Time [days]','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('MA osc - MA mean [deg]','FontWeight','bold');
grid on;

% NUMERICAL
figure;
tiledlayout(3,2);
nexttile;
plot(tplot, (Kepler_num(:,1)-Kepler_mean(:,1))*lunit);
xlabel('Time [days]','FontWeight','bold');
ylabel('[km]','FontWeight','bold');
title('a osc - a mean [km]','FontWeight','bold');
grid on;

nexttile;
plot(tplot, Kepler_num(:,2)-Kepler_mean(:,2));
xlabel('Time [days]','FontWeight','bold');
title('e osc - e mean','FontWeight','bold');
grid on;

nexttile;
plot(tplot, (Kepler_num(:,3)-Kepler_mean(:,3))*180/pi);
xlabel('Time [days]','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('i osc - i mean [deg]','FontWeight','bold');
grid on;

nexttile;
plot(tplot, (Kepler_num(:,4)-Kepler_mean(:,4))*180/pi);
xlabel('Time [days]','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('omega osc - omega mean [deg]','FontWeight','bold');
grid on;

nexttile;
plot(tplot, (Kepler_num(:,5)-Kepler_mean(:,5))*180/pi);
xlabel('Time [days]','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('RAAN osc - RAAN mean [deg]','FontWeight','bold');
grid on;

nexttile;
plot(tplot, (Kepler_num(:,6)-Kepler_mean(:,6))*180/pi);
xlabel('Time [days]','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('MA osc - MA mean [deg]','FontWeight','bold');
grid on;

%calculating long-period variations of Delaunay elements
Zonal_LP = zeros(length(tspan),6);
for i = 1:length(tspan)
    % Zonal Long-period variations
    Zonal_LP(i,:) = LPVarJ26(Delaunay_mean(i,:)', mu, re, Jcoeff);
end
%plots of LP variations
figure;
tiledlayout(3,2);
nexttile;
%a_fit = fit(tplot, Zonal_LP(:,1)*180/pi,'poly1');
plot(tplot, Zonal_LP(:,1)*180/pi, '-');
xlabel('Time [days]','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('l long period','FontWeight','bold');
grid on;

nexttile;
plot(tplot, Zonal_LP(:,2)*180/pi, '-');
xlabel('Time [days]','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('g long period','FontWeight','bold');
grid on;

nexttile;
plot(tplot, Zonal_LP(:,3)*180/pi, '-');
xlabel('Time [days]','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('h long period','FontWeight','bold');
grid on;

nexttile;
plot(tplot, Zonal_LP(:,4), '-');
xlabel('Time [days]','FontWeight','bold');
ylabel('[-]','FontWeight','bold');
title('L long period (nondim)','FontWeight','bold');
grid on;

nexttile;
plot(tplot, Zonal_LP(:,5), '-');
xlabel('Time [days]','FontWeight','bold');
ylabel('[-]','FontWeight','bold');
title('G long period (nondim)','FontWeight','bold');
grid on;

nexttile;
plot(tplot, Zonal_LP(:,6), '-');
xlabel('Time [days]','FontWeight','bold');
ylabel('[-]','FontWeight','bold');
title('H long period (nondim)','FontWeight','bold');
grid on;


%% Comparing numerical and analytical

j=1;
r_min = 1800.0;
r_max = 5000.0;
r_step = 200.0;

a_diff = zeros(1,(r_max-r_min)/r_step+1);
e_diff = zeros(1,(r_max-r_min)/r_step+1);
i_diff = zeros(1,(r_max-r_min)/r_step+1);
omega_diff = zeros(1,(r_max-r_min)/r_step+1);
RAAN_diff = zeros(1,(r_max-r_min)/r_step+1);
MA_diff = zeros(1,(r_max-r_min)/r_step+1);
r_diff = zeros(1,(r_max-r_min)/r_step+1);
v_diff = zeros(1,(r_max-r_min)/r_step+1);
initial_rad = r_min:r_step:r_max;

for zero_radius = r_min:r_step:r_max
    r0 = zero_radius;
    %ЧТО-ТО не так с маленькими ecc!!!
    sv0 = Kepl2Dec(mu, [r0/lunit,0.03,50*pi/180,50*pi/180,50*pi/180,0*pi/180]);
    Degree = 2;
    Order = 0;
    
    %Numerical propagation
    tspan = linspace(0,2*pi*sqrt(384748^3/398600.44)/tunit,30000);
    % non-normalized (normalized too)
    [T,sv] = ode113(@(T, sv) GravAcc_original(T, sv, 0, 0,we, Degree, Order, mu,re,Clm,Slm),tspan,sv0,optn);
    
    % numerical Kepler and Delaunay elements
    Kepler_num = zeros(length(T), 6);
    Delaunay_num = zeros(length(T), 6);
    for i = 1:length(T)
        Kepler_num(i,:) = Dec2Kepl(mu, sv(i,:));
        Delaunay_num(i,:) = Kepl2Del(mu, Kepler_num(i,:), false);
    end
    
    % Compute Absolute Initial Mean States
    
    % analytical Kepler and Delaunay elements at time = 0
    Kepler_an0 = Kepler_num(1,:)';
    Delaunay_an0 = Delaunay_num(1,:)';
    
    Delaunay_mean0 = Osc2Mean(Delaunay_an0, 0, 0, 0, Degree, Order, Clm, Slm, mu, re, we, tol, quadtol);
    
    %disp(['Mean IC Computation Done: ' num2str(toc)]);
    
    %Propagate Mean States
    
    Delaunay_mean = MeanProp(tspan, Delaunay_mean0, mu, re, Jcoeff, tol);
    
    Kepler_mean = zeros(length(tspan),6);
    for i = 1:length(tspan)
        Kepler_mean(i,:) = Kepl2Del( mu, Delaunay_mean(i,1:6)', true )'; %mean Kepler elements
    end
    
    %Mean to Osculating Transformation
    
    Delaunay_an = zeros(length(tspan),6);
    
    for i = 1:length(tspan)
        
        Delaunay_an(i,:) = Mean2Osc(Delaunay_mean(i,:)', 0, tspan(i), 0, Degree, Order, Clm, Slm, mu, re, we, tol, quadtol)';
        
    end
    
    %%%%% ДОБАВИТЬ В МАССИВЫ 0 элемент первой строчкой!
    
    sv_analytical = zeros(length(tspan),6);
    Kepler_an = zeros(length(tspan),6);
    for i = 1:length(tspan)
        Kepler_an(i,:) = Kepl2Del( mu, Delaunay_an(i,:)', true )';
        sv_analytical(i,:) = Kepl2Dec(mu, Kepler_an(i,:));
    end
    
    %calculated differences on the last step
    a_diff(j) = Kepler_num(end,1) - Kepler_an(end,1);
    e_diff(j) = Kepler_num(end,2) - Kepler_an(end,2);
    i_diff(j) = Kepler_num(end,3) - Kepler_an(end,3);
    omega_diff(j) = Kepler_num(end,4) - Kepler_an(end,4);
    RAAN_diff(j) = Kepler_num(end,5) - Kepler_an(end,5);
    MA_diff(j) = Kepler_num(end,6) - Kepler_an(end,6);
    r_diff(j) = norm(sv(end,1:3)) - norm(sv_analytical(end,1:3));
    v_diff(j) = norm(sv(end,4:6)) - norm(sv_analytical(end,4:6));
    
    j = j + 1;
end

figure;
tiledlayout(3,2)
nexttile;
plot(initial_rad, a_diff*lunit);
title('delta a, Km');
nexttile;
plot(initial_rad, e_diff);
title('delta e');
nexttile;
plot(initial_rad, i_diff*180/pi);
title('delta i, deg');
nexttile;
plot(initial_rad, omega_diff*180/pi);
title('delta omega, deg');
nexttile;
plot(initial_rad, RAAN_diff*180/pi);
title('delta RAAN, deg');
nexttile;
plot(initial_rad, MA_diff*180/pi);
title('delta MA, deg');

figure;
tiledlayout(2,1)
nexttile;
plot(initial_rad, r_diff*lunit);
title('delta r, Km');
nexttile;
plot(initial_rad, v_diff*vunit);
title('delta v, Km');

%% Comparing with Mahajan's analytic and numerical matrix

sv_analytical_norm = zeros(length(tspan),6);
sv_numerical_norm = zeros(length(tspan),6);
for i = 1:length(tspan)
    %normalizing for comparing with Mahajan's part
    sv_analytical_norm(i,:) = [sv_analytical(i,1:3)/lunit, sv_analytical(i,4:6)/vunit];
    sv_numerical_norm(i,:) = [sv(i,1:3)/lunit, sv(i,4:6)/vunit];
end

load('../NonEarthTheory-Terry/ARV.mat');
load('../NonEarthTheory-Terry/ORV.mat');
comparing_analytical = zeros(length(tspan),6);
comparing_numerical = zeros(length(tspan),6);
% for i = 1:length(tspan)
%     comparing_analytical(i,:) = ARV(i,:) - sv_analytical_norm(i,:);
%     comparing_numerical(i,:) = ORV(i,:) - sv_numerical_norm(i,:);
% end
for i = 1:length(tspan)
    comparing_analytical(i,:) = ARV(i,:) - sv_analytical(i,:);
    comparing_numerical(i,:) = ORV(i,:) - sv(i,:);
end

figure;
tiledlayout(3,2)
nexttile;
plot(tplot, comparing_analytical(:,1)*lunit);
title('delta X_a');
ylabel('Km','FontWeight','bold');
nexttile;
plot(tplot, comparing_analytical(:,2)*lunit);
title('delta Y_a');
ylabel('Km','FontWeight','bold');
nexttile;
plot(tplot, comparing_analytical(:,3)*lunit);
title('delta Z_a');
ylabel('Km','FontWeight','bold');
nexttile;
plot(tplot, comparing_analytical(:,4)*vunit);
title('delta Vx_a');
ylabel('Km/s','FontWeight','bold');
nexttile;
plot(tplot, comparing_analytical(:,5)*vunit);
title('delta Vy_a');
ylabel('Km/s','FontWeight','bold');
nexttile;
plot(tplot, comparing_analytical(:,6)*vunit);
title('delta Vz_a');
ylabel('Km/s','FontWeight','bold');

figure;
tiledlayout(3,2)
nexttile;
plot(tplot, comparing_numerical(:,1)*lunit);
title('delta X_n');
ylabel('Km','FontWeight','bold');
nexttile;
plot(tplot, comparing_numerical(:,2)*lunit);
title('delta Y_n');
ylabel('Km','FontWeight','bold');
nexttile;
plot(tplot, comparing_numerical(:,3)*lunit);
title('delta Z_n');
ylabel('Km','FontWeight','bold');
nexttile;
plot(tplot, comparing_numerical(:,4)*vunit);
title('delta Vx_n');
ylabel('Km/s','FontWeight','bold');
nexttile;
plot(tplot, comparing_numerical(:,5)*vunit);
title('delta Vy_n');
ylabel('Km/s','FontWeight','bold');
nexttile;
plot(tplot, comparing_numerical(:,6)*vunit);
title('delta Vz_n');
ylabel('Km/s','FontWeight','bold');

%% Local Functions

function J = ComputeError(tspan, Em0, Eto, GMST0, t0, AnICDegree, AnICOrder, Clm, Slm, mu, Re, we, tol, quadtol)

Jcoeff = -Clm(2:end,1);

%  mean elements
Em = MeanProp(tspan, Em0, mu, Re, Jcoeff, tol);

% % mean 2 osc
% Eo = zeros(length(tspan),6);
% for ctr = 1:length(tspan)
%     Eo(ctr,:) = Mean2Osc(Em(ctr,:)', GMST0, tspan(ctr), t0, AnICDegree, AnICOrder, Clm, Slm, mu, Re, we, tol, quadtol);
% end

% Error
J = 0;
for ctr = 4:6
    J = J + sum((Eto(:,ctr) - Em(:,ctr)).^2);
end

end


function Em = MeanProp(tspan, Em0, mu, Re, Jcoeff, tol)

m = length(tspan);

Emdot = SecularRates(Em0, mu, Re,  Jcoeff);

Em = zeros(m,6);

for ctr = 1:6
    Em(:,ctr) = Em0(ctr) + tspan*Emdot(ctr);
end

% optn = odeset('RelTol',tol,'AbsTol',tol*1e-3);

% numericall propagate short-period averaged mean elements
% [T,Em] = ode45(@(T, Em) SPAvgDiffEq(T, Em, mu,Re,Jcoeff),tspan,Em0,optn);

end

function Eo = Mean2Osc(Em, GMST0, t, t0, AnICDegree, AnICOrder, Clm, Slm, mu, Re, we, tol, quadtol)

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

end

function Em = Osc2Mean(Eo, GMST0, t, t0, AnICDegree, AnICOrder, Clm, Slm, mu, Re, we, tol, quadtol)

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
%Here I do 3 approximations in a row and not only 1
function Em = Osc2Mean_approx(Eo, GMST0, t, t0, AnICDegree, AnICOrder, Clm, Slm, mu, Re, we, tol, quadtol)

% ONLY IF TESSERALS ARE 0 (FOR NOW)

Jcoeff = -Clm(2:end,1);
GMST = GMST0 + we*(t - t0);
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

error = 1000;
j = 0;
% Substracting zonal Long-period
X_lp_k = X_sp_k_ - LPVarJ26(X_sp_k_, mu, Re, Jcoeff);
while error > 10^-12
    j = j + 1;
    
    % ONLY IF TESSERALS ARE 0
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
%             coe = [a,e,i,h,g,f]';
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