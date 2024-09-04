function [Mean] = Osc2Mean_thirdbody_new_one(Osc, Kepler_thirdbody, successive_approx)
%OSC2MEAN_THIRDBODY takes Osculating Delaunay elements and calculates Mean
%elements in order to remove a 3rd body perturbations

% Spacecraft parameters
l = Osc(1);
g = Osc(2);
h = Osc(3);
L = Osc(4);
G = Osc(5);
H = Osc(6);
e = sqrt(1-G^2/L^2);
E = Kepler_Eqn_solver(l, e, 10^-9);

% Earth's Kepler elements and other
a_ = Kepler_thirdbody(1);
e_ = Kepler_thirdbody(2);
inc_ = Kepler_thirdbody(3);
omega_ = Kepler_thirdbody(4);
RAAN_ = Kepler_thirdbody(5);
MA_ = Kepler_thirdbody(6);
E_ = Kepler_Eqn_solver(MA_, e_, 10^-10);
r_ = a_*(1-e_*cos(E_));
l_ = MA_;

% Calculating partials of W short per
delta = 10^-7;
eps = 1;
mu = 1;
%1)
W1_SP_plus = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, l+delta, g, h, L, G, H);
W1_SP_minus = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, l-delta, g, h, L, G, H);
dWdl_num = (W1_SP_plus-W1_SP_minus)/(2*delta);
%2)
W1_SP_plus = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, l, g+delta, h, L, G, H);
W1_SP_minus = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, l, g-delta, h, L, G, H);
dWdg_num = (W1_SP_plus-W1_SP_minus)/(2*delta);
%3)
W1_SP_plus = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, l, g, h+delta, L, G, H);
W1_SP_minus = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, l, g, h-delta, L, G, H);
dWdh_num = (W1_SP_plus-W1_SP_minus)/(2*delta);
%4)
W1_SP_plus = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, l, g, h, L+delta, G, H);
W1_SP_minus = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, l, g, h, L-delta, G, H);
dWdL_num = (W1_SP_plus-W1_SP_minus)/(2*delta);
%5)
W1_SP_plus = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, l, g, h, L, G+delta, H);
W1_SP_minus = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, l, g, h, L, G-delta, H);
dWdG_num = (W1_SP_plus-W1_SP_minus)/(2*delta);
%6)
W1_SP_plus = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, l, g, h, L, G, H+delta);
W1_SP_minus = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, l, g, h, L, G, H-delta);
dWdH_num = (W1_SP_plus-W1_SP_minus)/(2*delta);

dWdx_SP = [dWdl_num, dWdg_num, dWdh_num, dWdL_num, dWdG_num, dWdH_num];

% Short Period averaging
Mean_sp(1:3, 1) = Osc(1:3) + dWdx_SP(4:6)';
Mean_sp(4:6, 1) = Osc(4:6) - dWdx_SP(1:3)';

% Calculating partials of W long per
dWdx_LP(1) = 0;
dWdx_LP(2) = dWdg_LP(a_, r_, l_, Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6));
dWdx_LP(3) = dWdh_LP(a_, r_, l_, Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6));
dWdx_LP(4) = dWdL_LP(a_, r_, l_, Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6));
dWdx_LP(5) = dWdG_LP(a_, r_, l_, Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6));
dWdx_LP(6) = dWdH_LP(a_, r_, l_, Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6));

% Long Period averaging
Mean(1:3, 1) = Mean_sp(1:3) + dWdx_LP(4:6)';
Mean(4:6, 1) = Mean_sp(4:6) - dWdx_LP(1:3)';

% if successive_approx then approximates below
if successive_approx == 1
    error = 1000;
    k = 0;
    % Short-period
    %can use the partials from b4 here
    Mean_sp(1:3, 1) = Osc(1:3) + dWdx_SP(4:6)';
    Mean_sp(4:6, 1) = Osc(4:6) - dWdx_SP(1:3)';
    while (error > 10^-11) && (k < 100)
        k = k + 1;
        
        %1)
        W1_SP_plus = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, Mean_sp(1)+delta, Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6));
        W1_SP_minus = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, Mean_sp(1)-delta, Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6));
        dWdx_SP(1) = (W1_SP_plus-W1_SP_minus)/(2*delta);
        %2)
        W1_SP_plus = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, Mean_sp(1), Mean_sp(2)+delta, Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6));
        W1_SP_minus = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, Mean_sp(1), Mean_sp(2)-delta, Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6));
        dWdx_SP(2) = (W1_SP_plus-W1_SP_minus)/(2*delta);
        %3)
        W1_SP_plus = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, Mean_sp(1), Mean_sp(2), Mean_sp(3)+delta, Mean_sp(4), Mean_sp(5), Mean_sp(6));
        W1_SP_minus = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, Mean_sp(1), Mean_sp(2), Mean_sp(3)-delta, Mean_sp(4), Mean_sp(5), Mean_sp(6));
        dWdx_SP(3) = (W1_SP_plus-W1_SP_minus)/(2*delta);
        %4)
        W1_SP_plus = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, Mean_sp(1), Mean_sp(2), Mean_sp(3), Mean_sp(4)+delta, Mean_sp(5), Mean_sp(6));
        W1_SP_minus = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, Mean_sp(1), Mean_sp(2), Mean_sp(3), Mean_sp(4)-delta, Mean_sp(5), Mean_sp(6));
        dWdx_SP(4) = (W1_SP_plus-W1_SP_minus)/(2*delta);
        %5)
        W1_SP_plus = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, Mean_sp(1), Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5)+delta, Mean_sp(6));
        W1_SP_minus = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, Mean_sp(1), Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5)-delta, Mean_sp(6));
        dWdx_SP(5) = (W1_SP_plus-W1_SP_minus)/(2*delta);
        %6)
        W1_SP_plus = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, Mean_sp(1), Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6)+delta);
        W1_SP_minus = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, Mean_sp(1), Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6)-delta);
        dWdx_SP(6) = (W1_SP_plus-W1_SP_minus)/(2*delta);
        
        Mean_sp_(1:3,1) = Osc(1:3) + dWdx_SP(4:6)';
        Mean_sp_(4:6,1) = Osc(4:6) - dWdx_SP(1:3)';
        
        error = norm(Mean_sp_ - Mean_sp); %calculating error
        Mean_sp = Mean_sp_; %kth iteration equals (k+1)th
    end
    
    error = 1000;
    j = 0;
    % Long period
    dWdx_LP(1) = 0;
    dWdx_LP(2) = dWdg_LP(a_, r_, l_, Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6));
    dWdx_LP(3) = dWdh_LP(a_, r_, l_, Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6));
    dWdx_LP(4) = dWdL_LP(a_, r_, l_, Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6));
    dWdx_LP(5) = dWdG_LP(a_, r_, l_, Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6));
    dWdx_LP(6) = dWdH_LP(a_, r_, l_, Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6));
    Mean_lp(1:3, 1) = Mean_sp(1:3, 1) + dWdx_LP(4:6)';
    Mean_lp(4:6, 1) = Mean_sp(4:6, 1) - dWdx_LP(1:3)';
    while error > 10^-11 && (j < 50)
        j = j + 1;
        
        dWdx_LP(1) = 0;
        dWdx_LP(2) = dWdg_LP(a_, r_, l_, Mean_lp(2), Mean_lp(3), Mean_lp(4), Mean_lp(5), Mean_lp(6));
        dWdx_LP(3) = dWdh_LP(a_, r_, l_, Mean_lp(2), Mean_lp(3), Mean_lp(4), Mean_lp(5), Mean_lp(6));
        dWdx_LP(4) = dWdL_LP(a_, r_, l_, Mean_lp(2), Mean_lp(3), Mean_lp(4), Mean_lp(5), Mean_lp(6));
        dWdx_LP(5) = dWdG_LP(a_, r_, l_, Mean_lp(2), Mean_lp(3), Mean_lp(4), Mean_lp(5), Mean_lp(6));
        dWdx_LP(6) = dWdH_LP(a_, r_, l_, Mean_lp(2), Mean_lp(3), Mean_lp(4), Mean_lp(5), Mean_lp(6));
        Mean_lp_(1:3, 1) = Mean_sp(1:3, 1) + dWdx_LP(4:6)';
        Mean_lp_(4:6, 1) = Mean_sp(4:6, 1) - dWdx_LP(1:3)';
        
        error = norm(Mean_lp_ - Mean_lp); %calculating error
        Mean_lp = Mean_lp_; %kth iteration equals (k+1)th
    end
    
    Mean = Mean_lp_;
end

end

%% Local Functions
% Long Period Averaging
function dWdL = dWdL_LP(a_, r_, l_, g, h, L, G, H)
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

mu = 1;
n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);
e = sqrt(1-G^2/L^2);

koeff = mu_*n_/(2*mu^2);

dWdL = koeff*((a_/r_)^3*L/8*(3*G^2-10*L^2)*(2*l_*(1-3*H^2/G^2)+3*(H^2/G^2-1)*sin(2*l_-2*h)) ...
    + (a_/r_)^3*15/4*L/G^2*cos(2*g)*(l_-h)*(G^2-H^2)*(2*L^2-G^2) ...
    + (a_/r_)^3*15/8*L/G^2*(2*L^2-G^2)*cos(l_-h)*((G+H)^2*sin(l_-h-2*g)+(G-H)^2*sin(l_-h+2*g)) ...
    - l_/4*L/G^2*(G^2-3*H^2)*(3*G^2-10*L^2) ...
    +15/4*l_*L/G^2*(G^2-H^2)*(G^2-2*L^2)*cos(2*g));
end

function dWdG = dWdG_LP(a_, r_, l_, g, h, L, G, H)
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

mu = 1;
n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);
e = sqrt(1-G^2/L^2);

koeff = mu_*n_/(2*mu^2);

dWdG = koeff*((a_/r_)^3*3*L^2/8/G^3*(2*l_*(G^4-5*H^2*L^2)-(3*G^4-5*H^2*L^2)*sin(2*l_-2*h)) ...
    + (a_/r_)^3*15/4*L^2/G^3*cos(2*g)*(h-l_)*(G^4-H^2*L^2) ...
    - (a_/r_)^3*15/8*L^2/G*cos(l_-h)*((G+H)^2*sin(l_-h-2*g)+(G-H)^2*sin(l_-h+2*g)) ...
    - (a_/r_)^3*15/8*L^4/G^3*e^2*cos(l_-h)*((G+H)^2*sin(l_-h-2*g)+(G-H)^2*sin(l_-h+2*g)) ...
    + (a_/r_)^3*15/16*L^4/G^2*e^2*cos(l_-h)*(2*(G+H)*sin(l_-h-2*g)+2*(G-H)*sin(l_-h+2*g)) ...
    - 3/4*l_*L^2/G^3*(G^4-5*H^2*L^2) ...
    + 15/4*l_*L^2/G^3*(G^4-H^2*L^2)*cos(2*g));
end

function dWdH = dWdH_LP(a_, r_, l_, g, h, L, G, H)
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

mu = 1;
n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);
e = sqrt(1-G^2/L^2);

koeff = mu_*n_/(2*mu^2);

dWdH = koeff*H*L^4/G^2*((a_/r_)^3*3/8*(5-3*G^2/L^2)*(2*l_-sin(2*l_-2*h)) ...
    + (a_/r_)^3*15/4*e^2*cos(2*g)*(h-l_) ...
    + (a_/r_)^3*15/8/H*e^2*cos(l_-h)*((G+H)*sin(l_-h-2*g)-(G-H)*sin(l_-h+2*g)) ...
    + 3/4*l_*(3*G^2/L^2-5) ...
    + 15/4*l_*e^2*cos(2*g));
end

function dWdg = dWdg_LP(a_, r_, l_, g, h, L, G, H)
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

mu = 1;
n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);
e = sqrt(1-G^2/L^2);

koeff = mu_*n_/(2*mu^2);

dWdg = koeff*L^4*e^2*((a_/r_)^3*15/4*sin(2*g)*(h-l_)*(1-H^2/G^2) ...
    - (a_/r_)^3*15/8/G^2*cos(l_-h)*((G+H)^2*cos(l_-h-2*g)-(G-H)^2*cos(l_-h+2*g)) ...
    + l_*15/4*(1-H^2/G^2)*sin(2*g));
end

function dWdh = dWdh_LP(a_, r_, l_, g, h, L, G, H)
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

mu = 1;
n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);
e = sqrt(1-G^2/L^2);

koeff = mu_*n_/(2*mu^2);

dWdh = koeff*(a_/r_)^3*3/8*L^4*((H^2/G^2-1)*(5-3*G^2/L^2)*cos(2*l_-2*h) ...
    - 5*e^2*cos(2*g)*(1-H^2/G^2) ...
    + 5/2/G^2*e^2*(sin(l_-h)*((G+H)^2*sin(l_-h-2*g)+(G-H)^2*sin(l_-h+2*g))-cos(l_-h)*((G+H)^2*cos(l_-h-2*g)+(G-H)^2*cos(l_-h+2*g))));
end