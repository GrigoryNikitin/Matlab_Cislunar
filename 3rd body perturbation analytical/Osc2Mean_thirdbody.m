function [Mean] = Osc2Mean_thirdbody(Osc, Kepler_thirdbody, successive_approx)
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
dWdx_SP(1) = dWdl_SP(a_, r_, l_, g, h, L, G, H, E);
dWdx_SP(2) = dWdg_SP(a_, r_, l_, g, h, L, G, H, E);
dWdx_SP(3) = dWdh_SP(a_, r_, l_, g, h, L, G, H, E);
dWdx_SP(4) = dWdL_SP(a_, r_, l_, g, h, L, G, H, E);
dWdx_SP(5) = dWdG_SP(a_, r_, l_, g, h, L, G, H, E);
dWdx_SP(6) = dWdH_SP(a_, r_, l_, g, h, L, G, H, E);

% Short Period averaging
Mean_sp(1:3, 1) = Osc(1:3) - dWdx_SP(4:6)';
Mean_sp(4:6, 1) = Osc(4:6) + dWdx_SP(1:3)';


% Calculating partials of W long per
dWdx_LP(1) = 0;
dWdx_LP(2) = dWdg_LP(a_, r_, l_, Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6));
dWdx_LP(3) = dWdh_LP(a_, r_, l_, Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6));
dWdx_LP(4) = dWdL_LP(a_, r_, l_, Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6));
dWdx_LP(5) = dWdG_LP(a_, r_, l_, Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6));
dWdx_LP(6) = dWdH_LP(a_, r_, l_, Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6));

% Long Period averaging
Mean(1:3, 1) = Mean_sp(1:3) - dWdx_LP(4:6)';
Mean(4:6, 1) = Mean_sp(4:6) + dWdx_LP(1:3)';

% if successive_approx then approximates below
if successive_approx == 1
    error = 1000;
    k = 0;
    % Short-period
    %can use the partials from b4 here
    Mean_sp(1:3, 1) = Osc(1:3) - dWdx_SP(4:6)';
    Mean_sp(4:6, 1) = Osc(4:6) + dWdx_SP(1:3)';
    while error > 10^-12
        k = k + 1;
        
        E = Kepler_Eqn_solver(Mean_sp(1), e, 10^-9);
        dWdx_SP(1) = dWdl_SP(a_, r_, l_, Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6), E);
        dWdx_SP(2) = dWdg_SP(a_, r_, l_, Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6), E);
        dWdx_SP(3) = dWdh_SP(a_, r_, l_, Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6), E);
        dWdx_SP(4) = dWdL_SP(a_, r_, l_, Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6), E);
        dWdx_SP(5) = dWdG_SP(a_, r_, l_, Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6), E);
        dWdx_SP(6) = dWdH_SP(a_, r_, l_, Mean_sp(2), Mean_sp(3), Mean_sp(4), Mean_sp(5), Mean_sp(6), E);
        
        Mean_sp_(1:3,1) = Osc(1:3) - dWdx_SP(4:6)';
        Mean_sp_(4:6,1) = Osc(4:6) + dWdx_SP(1:3)';
        
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
    Mean_lp(1:3, 1) = Mean_sp(1:3, 1) - dWdx_LP(4:6)';
    Mean_lp(4:6, 1) = Mean_sp(4:6, 1) + dWdx_LP(1:3)';
    while error > 10^-12
        j = j + 1;
        
        dWdx_LP(1) = 0;
        dWdx_LP(2) = dWdg_LP(a_, r_, l_, Mean_lp(2), Mean_lp(3), Mean_lp(4), Mean_lp(5), Mean_lp(6));
        dWdx_LP(3) = dWdh_LP(a_, r_, l_, Mean_lp(2), Mean_lp(3), Mean_lp(4), Mean_lp(5), Mean_lp(6));
        dWdx_LP(4) = dWdL_LP(a_, r_, l_, Mean_lp(2), Mean_lp(3), Mean_lp(4), Mean_lp(5), Mean_lp(6));
        dWdx_LP(5) = dWdG_LP(a_, r_, l_, Mean_lp(2), Mean_lp(3), Mean_lp(4), Mean_lp(5), Mean_lp(6));
        dWdx_LP(6) = dWdH_LP(a_, r_, l_, Mean_lp(2), Mean_lp(3), Mean_lp(4), Mean_lp(5), Mean_lp(6));
        Mean_lp_(1:3, 1) = Mean_sp(1:3, 1) - dWdx_LP(4:6)';
        Mean_lp_(4:6, 1) = Mean_sp(4:6, 1) + dWdx_LP(1:3)';
        
        error = norm(Mean_lp_ - Mean_lp); %calculating error
        Mean_lp = Mean_lp_; %kth iteration equals (k+1)th
    end
    
    Mean = Mean_lp_;
end

end

%% Local Functions

% Short Period Averaging
function dWdL = dWdL_SP(a_, r_, l_, g, h, L, G, H, E)
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

mu = 1;
n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);
a = L^2/mu;
if abs(G) > abs(L)
    e = 0;
else
    e = sqrt(1-G^2/L^2);
end
r = a*(1-e*cos(E));

alpha = cos(g)*cos(l_-h) + sin(g)*sin(l_-h)*H/G;
beta = -sin(g)*cos(l_-h) + cos(g)*sin(l_-h)*H/G;

koeff = mu_*n_^2/(2*mu^4)*(a_/r_)^3*1;

dWdL = koeff*L^6*(1/3*(1-3*beta^2)*sin(3*E)*e*(7/4-G^2/L^2) ...
    + 3*sin(E)*(4*alpha^2-beta^2-1)*(7/4-G^2/L^2)*e ...
    + 15/4*sin(2*E)*(2*alpha^2+beta^2-1)*(7/5-G^2/L^2) ...
    + alpha*beta*G/L*(6*cos(2*E)*(G^2/L^2-3)+5/2/e*(15*cos(E)+cos(3*E))*(6/5-G^2/L^2)) ...
    + 3/2/e*(27*alpha^2-3*beta^2-8)*sin(E)*(G^2/L^2-7/6) ...
    + 3/2/e*(alpha^2-beta^2)*sin(3*E)*(G^2/L^2-7/6) ...
    + 21/4*(alpha^2-beta^2)*sin(2*E));

dWdE = koeff*L^7*(3*(1/12-1/4*beta^2)*e^3*cos(3*E) ...
    + (3*alpha^2-3/4*beta^2-3/4)*e^3*cos(E) ...
    + (3*alpha^2+3/2*beta^2-3/2)*e^2*cos(2*E) ...
    + 3*(e^2+1)*alpha*beta*G/L*sin(2*E) - 3/2*alpha*beta*e*(5*sin(E)+sin(3*E))*G/L ...
    + 1/4*(-27*alpha^2 + 3*beta^2+8)*cos(E)*e ...
    + 3/4*(beta^2-alpha^2)*cos(3*E)*e ...
    + 3/2*(alpha^2-beta^2)*cos(2*E));

dEdL = a/r*G^2/e/L^3*sin(E);

dWdL = dWdL + dWdE*dEdL;
end

function dWdG = dWdG_SP(a_, r_, l_, g, h, L, G, H, E)
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

mu = 1;
n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);
a = L^2/mu;
if abs(G) > abs(L)
    e = 0;
else
    e = sqrt(1-G^2/L^2);
end
r = a*(1-e*cos(E));

alpha = cos(g)*cos(l_-h) + sin(g)*sin(l_-h)*H/G;
beta = -sin(g)*cos(l_-h) + cos(g)*sin(l_-h)*H/G;
dalphadG = -H/G^2*sin(g)*sin(l_-h);
dbetadG = -H/G^2*cos(g)*sin(l_-h);

koeff = mu_*n_^2/(2*mu^4)*(a_/r_)^3;

dWdG = koeff*L^5*(1/4*(3*beta^2-1)*G*sin(3*E)*e ...
    - 9/4*G*sin(E)*(4*alpha^2-beta^2-1)*e ...
    -3/2*G*sin(2*E)*(2*alpha^2+beta^2-1) ...
    +9/2*alpha*beta*L*cos(2*E)*(G^2/L^2-2/3)+1/e*alpha*beta*L*(15*cos(E)+cos(3*E))*(1/2-G^2/L^2) ...
    - G/4/e*(-27*alpha^2+3*beta^2+8)*sin(E) ...
    +G/4/e*(alpha^2-beta^2)*sin(3*E));

dWdalpha = koeff*L^7*(6*alpha*sin(E)*e^3 + 3*alpha*sin(2*E)*e^2 ...
    - 3/2*beta*G/L*cos(2*E)*(1+e^2) + 1/2*beta*G/L*(15*cos(E)+cos(3*E))*e ...
    - 27/2*alpha*sin(E)*e - 1/2*alpha*sin(3*E)*e + 3/2*alpha*sin(2*E));

dWdbeta = koeff*L^7*(-1/2*beta*sin(3*E)*e^3 - 3/2*beta*sin(E)*e^3 ...
    +3/2*beta*sin(2*E)*e^2 - 3/2*alpha*G/L*cos(2*E)*(2-G^2/L^2) ...
    + 1/2*alpha*G/L*(15*cos(E)+cos(3*E))*e + 3/2*beta*sin(E)*e ...
    +1/2*beta*sin(3*E)*e - 3/2*beta*sin(2*E));

dWdE = koeff*L^7*(3*(1/12-1/4*beta^2)*e^3*cos(3*E) ...
    + (3*alpha^2-3/4*beta^2-3/4)*e^3*cos(E) ...
    + (3*alpha^2+3/2*beta^2-3/2)*e^2*cos(2*E) ...
    + 3*(e^2+1)*alpha*beta*G/L*sin(2*E) - 3/2*alpha*beta*e*(5*sin(E)+sin(3*E))*G/L ...
    + 1/4*(-27*alpha^2 + 3*beta^2+8)*cos(E)*e ...
    + 3/4*(beta^2-alpha^2)*cos(3*E)*e ...
    + 3/2*(alpha^2-beta^2)*cos(2*E));

dEdG = -a/r*G/e/L^2*sin(E);

dWdG = dWdG + dWdalpha*dalphadG + dWdbeta*dbetadG + dWdE*dEdG;
end

function dWdH = dWdH_SP(a_, r_, l_, g, h, L, G, H, E)
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

mu = 1;
n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);
if abs(G) > abs(L)
    e = 0;
else
    e = sqrt(1-G^2/L^2);
end

alpha = cos(g)*cos(l_-h) + sin(g)*sin(l_-h)*H/G;
beta = -sin(g)*cos(l_-h) + cos(g)*sin(l_-h)*H/G;
dalphadH = 1/G*sin(g)*sin(l_-h);
dbetadH = 1/G*cos(g)*sin(l_-h);

koeff = mu_*n_^2/(2*mu^4)*(a_/r_)^3;

dWdalpha = koeff*L^7*(6*alpha*sin(E)*e^3 + 3*alpha*sin(2*E)*e^2 ...
    - 3/2*beta*G/L*cos(2*E)*(1+e^2) + 1/2*beta*G/L*(15*cos(E)+cos(3*E))*e ...
    - 27/2*alpha*sin(E)*e - 1/2*alpha*sin(3*E)*e + 3/2*alpha*sin(2*E));

dWdbeta = koeff*L^7*(-1/2*beta*sin(3*E)*e^3 - 3/2*beta*sin(E)*e^3 ...
    +3/2*beta*sin(2*E)*e^2 - 3/2*alpha*G/L*cos(2*E)*(2-G^2/L^2) ...
    + 1/2*alpha*G/L*(15*cos(E)+cos(3*E))*e + 3/2*beta*sin(E)*e ...
    +1/2*beta*sin(3*E)*e - 3/2*beta*sin(2*E));

dWdH = 0;

dWdH = dWdH + dWdalpha*dalphadH + dWdbeta*dbetadH;
end

function dWdl = dWdl_SP(a_, r_, l_, g, h, L, G, H, E)
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

mu = 1;
n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);
if abs(G) > abs(L)
    e = 0;
else
    e = sqrt(1-G^2/L^2);
end

alpha = cos(g)*cos(l_-h) + sin(g)*sin(l_-h)*H/G;
beta = -sin(g)*cos(l_-h) + cos(g)*sin(l_-h)*H/G;

koeff = mu_*n_^2/(2*mu^4)*(a_/r_)^3;

dWdE = koeff*L^7*(3*(1/12-1/4*beta^2)*e^3*cos(3*E) ...
    + (3*alpha^2-3/4*beta^2-3/4)*e^3*cos(E) ...
    + (3*alpha^2+3/2*beta^2-3/2)*e^2*cos(2*E) ...
    + 3*(e^2+1)*alpha*beta*G/L*sin(2*E) - 3/2*alpha*beta*e*(5*sin(E)+sin(3*E))*G/L ...
    + 1/4*(-27*alpha^2 + 3*beta^2+8)*cos(E)*e ...
    + 3/4*(beta^2-alpha^2)*cos(3*E)*e ...
    + 3/2*(alpha^2-beta^2)*cos(2*E));

dWdl = 1/(1-e*cos(E))*dWdE;
end

function dWdg = dWdg_SP(a_, r_, l_, g, h, L, G, H, E)
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

mu = 1;
n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);
if abs(G) > abs(L)
    e = 0;
else
    e = sqrt(1-G^2/L^2);
end

alpha = cos(g)*cos(l_-h) + sin(g)*sin(l_-h)*H/G;
beta = -sin(g)*cos(l_-h) + cos(g)*sin(l_-h)*H/G;
dalphadg = -sin(g)*cos(l_-h)+cos(g)*sin(l_-h)*H/G;
dbetadg = -cos(g)*cos(l_-h)-sin(g)*sin(l_-h)*H/G;

koeff = mu_*n_^2/(2*mu^4)*(a_/r_)^3;

dWdalpha = koeff*L^7*(6*alpha*sin(E)*e^3 + 3*alpha*sin(2*E)*e^2 ...
    - 3/2*beta*G/L*cos(2*E)*(1+e^2) + 1/2*beta*G/L*(15*cos(E)+cos(3*E))*e ...
    - 27/2*alpha*sin(E)*e - 1/2*alpha*sin(3*E)*e + 3/2*alpha*sin(2*E));

dWdbeta = koeff*L^7*(-1/2*beta*sin(3*E)*e^3 - 3/2*beta*sin(E)*e^3 ...
    +3/2*beta*sin(2*E)*e^2 - 3/2*alpha*G/L*cos(2*E)*(2-G^2/L^2) ...
    + 1/2*alpha*G/L*(15*cos(E)+cos(3*E))*e + 3/2*beta*sin(E)*e ...
    +1/2*beta*sin(3*E)*e - 3/2*beta*sin(2*E));

dWdg = 0;

dWdg = dWdg + dWdalpha*dalphadg + dWdbeta*dbetadg;
end

function dWdh = dWdh_SP(a_, r_, l_, g, h, L, G, H, E)
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

mu = 1;
n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);
if abs(G) > abs(L)
    e = 0;
else
    e = sqrt(1-G^2/L^2);
end

alpha = cos(g)*cos(l_-h) + sin(g)*sin(l_-h)*H/G;
beta = -sin(g)*cos(l_-h) + cos(g)*sin(l_-h)*H/G;
dalphadh = cos(g)*sin(l_-h)-sin(g)*cos(l_-h)*H/G;
dbetadh = -sin(g)*sin(l_-h)-cos(g)*cos(l_-h)*H/G;

koeff = mu_*n_^2/(2*mu^4)*(a_/r_)^3;

dWdalpha = koeff*L^7*(6*alpha*sin(E)*e^3 + 3*alpha*sin(2*E)*e^2 ...
    - 3/2*beta*G/L*cos(2*E)*(1+e^2) + 1/2*beta*G/L*(15*cos(E)+cos(3*E))*e ...
    - 27/2*alpha*sin(E)*e - 1/2*alpha*sin(3*E)*e + 3/2*alpha*sin(2*E));

dWdbeta = koeff*L^7*(-1/2*beta*sin(3*E)*e^3 - 3/2*beta*sin(E)*e^3 ...
    +3/2*beta*sin(2*E)*e^2 - 3/2*alpha*G/L*cos(2*E)*(2-G^2/L^2) ...
    + 1/2*alpha*G/L*(15*cos(E)+cos(3*E))*e + 3/2*beta*sin(E)*e ...
    +1/2*beta*sin(3*E)*e - 3/2*beta*sin(2*E));

dWdh = 0;

dWdh = dWdh + dWdalpha*dalphadh + dWdbeta*dbetadh;
end

% Long Period Averaging
function dWdL = dWdL_LP(a_, r_, l_, g, h, L, G, H)
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

mu = 1;
n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);
if abs(G) > abs(L)
    e = 0;
else
    e = sqrt(1-G^2/L^2);
end

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
if abs(G) > abs(L)
    e = 0;
else
    e = sqrt(1-G^2/L^2);
end

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
if abs(G) > abs(L)
    e = 0;
else
    e = sqrt(1-G^2/L^2);
end

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
if abs(G) > abs(L)
    e = 0;
else
    e = sqrt(1-G^2/L^2);
end

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
if abs(G) > abs(L)
    e = 0;
else
    e = sqrt(1-G^2/L^2);
end

koeff = mu_*n_/(2*mu^2);

dWdh = koeff*(a_/r_)^3*3/8*L^4*((H^2/G^2-1)*(5-3*G^2/L^2)*cos(2*l_-2*h) ...
    - 5*e^2*cos(2*g)*(1-H^2/G^2) ...
    + 5/2/G^2*e^2*(sin(l_-h)*((G+H)^2*sin(l_-h-2*g)+(G-H)^2*sin(l_-h+2*g))-cos(l_-h)*((G+H)^2*cos(l_-h-2*g)+(G-H)^2*cos(l_-h+2*g))));
end