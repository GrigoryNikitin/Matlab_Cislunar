clear;
clc;

% ODE tolerances
tol = 1e-13;
quadtol = 1e-6;
optn = odeset('RelTol',tol,'AbsTol',tol*1e-3);

%units for nondimensionalization
runit = 384400;
tunit = 4.348*24*60*60;
vunit = runit/tunit;
mu  = 0.0121505856;

%initial conditions in the rotating nondim ref frame
sv1 = [0.906540675867110, 0, 0, 0, -0.237268185992377, -0.460110729423766];
%sv1 = [43000/runit, 0, 0, 0, sqrt(398600.44/43000)/vunit, 0];
T1 = 2*2.009450749673591;
sv2 = [0.945800000000000, 0, 0.163004475105074, 0, -0.038800217236384, 0];
T2 = 2*1.464302315134892;
sv3 = [1.021591576694260, 0, 0, 0, 0.801438324933024, 0];
T3 = 2*2.458091271152401;

%propagations
number_of_points = 5000;

timespan = linspace(0,T1,number_of_points);
[tspan,rot_sv_nondim1] = ode113(@(tspan, rot_sv_nondim1) derivatives_CR3BP_nondim(tspan, rot_sv_nondim1),timespan,sv1,optn);

timespan = linspace(0,T2,number_of_points);
[tspan,rot_sv_nondim2] = ode113(@(tspan, rot_sv_nondim2) derivatives_CR3BP_nondim(tspan, rot_sv_nondim2),timespan,sv2,optn);

timespan = linspace(0,T3,number_of_points);
[tspan,rot_sv_nondim3] = ode113(@(tspan, rot_sv_nondim3) derivatives_CR3BP_nondim(tspan, rot_sv_nondim3),timespan,sv3,optn);


%calc inertial coords
rE1 = zeros(number_of_points, 3);
rM1 = zeros(number_of_points, 3);
rE2 = zeros(number_of_points, 3);
rM2 = zeros(number_of_points, 3);
rE3 = zeros(number_of_points, 3);
rM3 = zeros(number_of_points, 3);
for i = 1:number_of_points
    [rE1(i,:), rM1(i,:), inert_sv1_nond(i,:)] = rot2inert(tspan(i),mu,rot_sv_nondim1(i,:));
    [rE2(i,:), rM2(i,:), inert_sv2_nond(i,:)] = rot2inert(tspan(i),mu,rot_sv_nondim2(i,:));
    [rE3(i,:), rM3(i,:), inert_sv3_nond(i,:)] = rot2inert(tspan(i),mu,rot_sv_nondim3(i,:));
end
%adding dimensions
for i = 1:number_of_points
    rot_sv1(i, 1:3) = rot_sv_nondim1(i,1:3)*runit;
    rot_sv1(i, 4:6) = rot_sv_nondim1(i,4:6)*vunit;
    rot_sv2(i, 1:3) = rot_sv_nondim2(i,1:3)*runit;
    rot_sv2(i, 4:6) = rot_sv_nondim2(i,4:6)*vunit;
    rot_sv3(i, 1:3) = rot_sv_nondim3(i,1:3)*runit;
    rot_sv3(i, 4:6) = rot_sv_nondim3(i,4:6)*vunit;
    
    rE1(i,:) = rE1(i,:)*runit;
    rM1(i,:) = rM1(i,:)*runit;
    inert_sv1(i, 1:3) = inert_sv1_nond(i,1:3)*runit;
    inert_sv1(i, 4:6) = inert_sv1_nond(i,4:6)*vunit;
    rE2(i,:) = rE2(i,:)*runit;
    rM2(i,:) = rM2(i,:)*runit;
    inert_sv2(i, 1:3) = inert_sv2_nond(i,1:3)*runit;
    inert_sv2(i, 4:6) = inert_sv2_nond(i,4:6)*vunit;
    rE3(i,:) = rE3(i,:)*runit;
    rM3(i,:) = rM3(i,:)*runit;
    inert_sv3(i, 1:3) = inert_sv3_nond(i,1:3)*runit;
    inert_sv3(i, 4:6) = inert_sv3_nond(i,4:6)*vunit;
end


%% Plots
figure;
tiledlayout(2,1);
nexttile;
plot(rot_sv1(:,1), rot_sv1(:,2));
hold on;
plot(-mu*runit, 0, 'o');
plot((1-mu)*runit, 0, 'o');
%axis equal;
hold off;
nexttile
plot(inert_sv1(:,1), inert_sv1(:,2));
hold on;
plot(rE1(:,1),rE1(:,2));
plot(rM1(:,1),rM1(:,2));
axis equal;
hold off;

figure;
tiledlayout(2,1);
nexttile;
plot(rot_sv2(:,1), rot_sv2(:,2));
hold on;
plot(-mu*runit, 0, 'o');
plot((1-mu)*runit, 0, 'o');
%axis equal;
hold off;
nexttile
plot(inert_sv2(:,1), inert_sv2(:,2));
hold on;
plot(rE2(:,1),rE2(:,2));
plot(rM2(:,1),rM2(:,2));
axis equal;
hold off;

figure;
tiledlayout(2,1);
nexttile;
plot(rot_sv3(:,1), rot_sv3(:,2));
hold on;
plot(-mu*runit, 0, 'o');
plot((1-mu)*runit, 0, 'o');
%axis equal;
hold off;
nexttile
plot(inert_sv3(:,1), inert_sv3(:,2));
hold on;
plot(rE3(:,1),rE3(:,2));
plot(rM3(:,1),rM3(:,2));
axis equal;
hold off;