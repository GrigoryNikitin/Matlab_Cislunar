function error = compare_calculations(Delaunay,Kepler,Jcoeff)
%COMPARE_CALCULATIONS Summary of this function goes here
%   Comparing the calculation of domega/dt mine and Bharat's
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

mu = 1;
a_ = 221.3827754323089;
n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);

l = Delaunay(1);
g = Delaunay(2);
h = Delaunay(3);
L = Delaunay(4);
G = Delaunay(5);
H = Delaunay(6);

a = Kepler(1);
e = Kepler(2);
inc = Kepler(3);
omega = Kepler(4);
RAAN = Kepler(5);
MA = Kepler(6);

J(2:6) = Jcoeff(2:6);
%J(2) = 0;
%J(4) = 0;
%J(6) = 0;
R_e = 1;
R__e = 1;
J_2 = 0.000203223560871;
J_4 = -9.70434253277103*10^-6;
J_6 = -1.37675387844709*10^-5;
n = sqrt(mu/a^3);

%third body part from secular rates:
domegadt_1 = -3*mu_*n_^2/8/mu^2*L^2*(5*H^2*L^2/G^3*(cos(2*g)-1)-5*G*cos(2*g)+G);
%third body part from frozen orb:
%CHANGED THE SIGN!
domegadt_2 = -3*mu_*n_^2/(8*n)*(5*(cos(2*omega)-1)/sqrt(1-e^2)*cos(inc)^2 - 5*sqrt(1-e^2)*cos(2*omega)+sqrt(1-e^2));

%the grav part from Bharat:
domegadt_3 = -0.2625e4 / 0.2048e4 * ((G ^ 10) + ((-27 * H ^ 2 - 6 * L ^ 2) * G ^ 8) + ((99 * H ^ 4) + (154 * H ^ 2 * L ^ 2) + 0.33e2 / 0.5e1 * (L ^ 4)) * (G ^ 6) + (-0.429e3 / 0.5e1 * (H ^ 6) - (546 * H ^ 4 * L ^ 2) - 0.819e3 / 0.5e1 * (H ^ 2) * (L ^ 4)) * (G ^ 4) + ((462 * H ^ 6 * L ^ 2 + 567 * H ^ 4 * L ^ 4) * G ^ 2) - 0.11781e5 / 0.25e2 * (H ^ 6) * (L ^ 4)) * mu ^ 8 / (L ^ 7) / (G ^ 18) * J(6) * R__e ^ 6 + 0.15e2 / 0.128e3 / (L ^ 5) * mu ^ 6 / (G ^ 12) * (9 * G ^ 6 - 126 * G ^ 4 * H ^ 2 - 21 * G ^ 4 * L ^ 2 + 189 * G ^ 2 * H ^ 4 + 270 * G ^ 2 * H ^ 2 * L ^ 2 - 385 * H ^ 4 * L ^ 2) * J(4) * R__e ^ 4 - 0.3e1 / 0.4e1 / (G ^ 6) / (L ^ 3) * (G ^ 2 - 5 * H ^ 2) * mu ^ 4 * J(2) * R__e ^ 2;

%the grav part I got:
%CHANGED THE SIGNS HERE (???) TOO
%NEED TO REPLACE J6
t = cos(inc)^2;
domegadt_4_correct = - J_2*3*n*R_e^2/(4*a^2*(1-e^2)^2)*(1-5*t) ...
    + J_4*15*n*R_e^4/(128*a^4*(1-e^2)^4)*(9*(1-e^2)-21*(6*(1-e^2)*t+1)+27*t*(7*(1-e^2)*t+10)-385*t^2) ...
    + J_6*5*n*R_e^6/(2048*(1-e^2)^6*a^6)*((5-105*t+315*t^2-231*t^3)*(60*(1-e^2)^2-140*(1-e^2)) ...
    + (210*t-1260*t^2+1386*t^3)*(63+15*(1-e^2)^2-70*(1-e^2)) ...
    - 11*(5-105*t+315*t^2-231*t^3)*(63+15*(1-e^2)^2-70*(1-e^2)));
%CHANGED THE SIGN
%CORRECT J6:
domegadt_4 = -J_6*5*mu^8*R_e^6/(2048*G^12*L^3)*((5-105*H^2/G^2+315*H^4/G^4-231*H^6/G^6)*(60*G^4/L^4-140*G^2/L^2) ...
    + (210*H^2/G^2-1260*H^4/G^4+1386*H^6/G^6)*(63+15*G^4/L^4-70*G^2/L^2) ...
    - 11*(5-105*H^2/G^2+315*H^4/G^4-231*H^6/G^6)*(63+15*G^4/L^4-70*G^2/L^2));
domegadt_4_ = -J_6*5*n*R_e^6/(2048*(1-e^2)^6*a^6)*((5-105*t+315*t^2-231*t^3)*(60*(1-e^2)^2-140*(1-e^2)) ...
    + (210*t-1260*t^2+1386*t^3)*(63+15*(1-e^2)^2-70*(1-e^2)) ...
    - 11*(5-105*t+315*t^2-231*t^3)*(63+15*(1-e^2)^2-70*(1-e^2)));

error = 'abc';
end

