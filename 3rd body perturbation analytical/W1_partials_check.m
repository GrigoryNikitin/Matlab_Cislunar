% Here I calculate W1_SP and W1_LP third body partials numerically and compare it to
% the analytical results
%ONLY FOR NONDIMENSIONALIZED EQTS
clear;
clc;

global Mass_Earth;
global Mass_Moon;
global mu_Moon;
global mu_Earth;
Mass_Earth = 5.9722*10^24;
Mass_Moon = 7.3477*10^22;
mu_Moon = 4902.8001224453001;
mu_Earth = 6.674*10^-20*Mass_Earth;

eps = 1;
r_ = 221;
a_ = 221;
l_ = 1;
l = 40*pi/180;
g = 30*pi/180;
h = 50*pi/180;
mu = 1;
L = 1.0177;
G = 1.0172;
H = 0.8809;
%Needed for the derivatives:
e = sqrt(1-G^2/L^2);
E = Kepler_Eqn_solver(l, e, 10^-9);
n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3);
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);
alpha = cos(g)*cos(l_-h) + sin(g)*sin(l_-h)*H/G;
beta = -sin(g)*cos(l_-h) + cos(g)*sin(l_-h)*H/G;
dalphadG = -H/G^2*sin(g)*sin(l_-h);
dbetadG = -H/G^2*cos(g)*sin(l_-h);
%%%

W1_SP_new_one = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, l, g, h, L, G, H);
W1_SP_orig = W1_thirdbody_SP(eps, a_, r_, mu, l_, l, g, h, L, G, H);
W1_LP_orig = W1_thirdbody_LP(eps, a_, r_, mu, l_, g, h, L, G, H);

delta = 10^-7;

% Short Period

%1)
W1_SP_plus = W1_thirdbody_SP(eps, a_, r_, mu, l_, l, g, h, L+delta, G, H);
W1_SP_minus = W1_thirdbody_SP(eps, a_, r_, mu, l_, l, g, h, L-delta, G, H);
dWdL_num = (W1_SP_plus-W1_SP_minus)/(2*delta);

% Analitycal dW/dL:
dWdL_an = dWdL_SP(eps, a_, r_, mu, l_, g, h, L, G, H, E);
%%%

%2)
W1_SP_plus = W1_thirdbody_SP(eps, a_, r_, mu, l_, l, g, h, L, G+delta, H);
W1_SP_minus = W1_thirdbody_SP(eps, a_, r_, mu, l_, l, g, h, L, G-delta, H);
dWdG_num = (W1_SP_plus-W1_SP_minus)/(2*delta);

% Analitycal dW/dG:
dWdG_an = dWdG_SP(eps, a_, r_, mu, l_, g, h, L, G, H, E);
%%%

%3)
W1_SP_plus = W1_thirdbody_SP(eps, a_, r_, mu, l_, l, g, h, L, G, H+delta);
W1_SP_minus = W1_thirdbody_SP(eps, a_, r_, mu, l_, l, g, h, L, G, H-delta);
dWdH_num = (W1_SP_plus-W1_SP_minus)/(2*delta);

% Analitycal dW/dH:
dWdH_an = dWdH_SP(eps, a_, r_, mu, l_, g, h, L, G, H, E);
%%%

%4)
W1_SP_plus = W1_thirdbody_SP(eps, a_, r_, mu, l_, l+delta, g, h, L, G, H);
W1_SP_minus = W1_thirdbody_SP(eps, a_, r_, mu, l_, l-delta, g, h, L, G, H);
dWdl_num = (W1_SP_plus-W1_SP_minus)/(2*delta);

% Analitycal dW/dl:
dWdl_an = dWdl_SP(eps, a_, r_, mu, l_, g, h, L, G, H, E);
%%%

%5)
W1_SP_plus = W1_thirdbody_SP(eps, a_, r_, mu, l_, l, g+delta, h, L, G, H);
W1_SP_minus = W1_thirdbody_SP(eps, a_, r_, mu, l_, l, g-delta, h, L, G, H);
dWdg_num = (W1_SP_plus-W1_SP_minus)/(2*delta);

% Analitycal dW/dg:
dWdg_an = dWdg_SP(eps, a_, r_, mu, l_, g, h, L, G, H, E);
%%%

%6)
W1_SP_plus = W1_thirdbody_SP(eps, a_, r_, mu, l_, l, g, h+delta, L, G, H);
W1_SP_minus = W1_thirdbody_SP(eps, a_, r_, mu, l_, l, g, h-delta, L, G, H);
dWdh_num = (W1_SP_plus-W1_SP_minus)/(2*delta);

% Analitycal dW/dh:
dWdh_an = dWdh_SP(eps, a_, r_, mu, l_, g, h, L, G, H, E);
%%%

deltaL = dWdL_an - dWdL_num;
deltaG = dWdG_an - dWdG_num;
deltaH = dWdH_an - dWdH_num;
deltal = dWdl_an - dWdl_num;
deltag = dWdg_an - dWdg_num;
deltah = dWdh_an - dWdh_num;


% Long Period

%1)
W1_LP_plus = W1_thirdbody_LP(eps, a_, r_, mu, l_, g, h, L+delta, G, H);
W1_LP_minus = W1_thirdbody_LP(eps, a_, r_, mu, l_, g, h, L-delta, G, H);
dWdL_num_LP = (W1_LP_plus-W1_LP_minus)/(2*delta);

% Analitycal dW/dL:
dWdL_an_LP = dWdL_LP(eps, a_, r_, mu, l_, g, h, L, G, H);
%%%

%2)
W1_LP_plus = W1_thirdbody_LP(eps, a_, r_, mu, l_, g, h, L, G+delta, H);
W1_LP_minus = W1_thirdbody_LP(eps, a_, r_, mu, l_, g, h, L, G-delta, H);
dWdG_num_LP = (W1_LP_plus-W1_LP_minus)/(2*delta);

% Analitycal dW/dG:
dWdG_an_LP = dWdG_LP(eps, a_, r_, mu, l_, g, h, L, G, H);
%%%

%3)
W1_LP_plus = W1_thirdbody_LP(eps, a_, r_, mu, l_, g, h, L, G, H+delta);
W1_LP_minus = W1_thirdbody_LP(eps, a_, r_, mu, l_, g, h, L, G, H-delta);
dWdH_num_LP = (W1_LP_plus-W1_LP_minus)/(2*delta);

% Analitycal dW/dH:
dWdH_an_LP = dWdH_LP(eps, a_, r_, mu, l_, g, h, L, G, H);
%%%

%5)
W1_LP_plus = W1_thirdbody_LP(eps, a_, r_, mu, l_, g+delta, h, L, G, H);
W1_LP_minus = W1_thirdbody_LP(eps, a_, r_, mu, l_, g-delta, h, L, G, H);
dWdg_num_LP = (W1_LP_plus-W1_LP_minus)/(2*delta);

% Analitycal dW/dg:
dWdg_an_LP = dWdg_LP(eps, a_, r_, mu, l_, g, h, L, G, H);
%%%

%6)
W1_LP_plus = W1_thirdbody_LP(eps, a_, r_, mu, l_, g, h+delta, L, G, H);
W1_LP_minus = W1_thirdbody_LP(eps, a_, r_, mu, l_, g, h-delta, L, G, H);
dWdh_num_LP = (W1_LP_plus-W1_LP_minus)/(2*delta);

% Analitycal dW/dh:
dWdh_an_LP = dWdh_LP(eps, a_, r_, mu, l_, g, h, L, G, H);
%%%

deltaL_LP = dWdL_an_LP - dWdL_num_LP;
deltaG_LP = dWdG_an_LP - dWdG_num_LP;
deltaH_LP = dWdH_an_LP - dWdH_num_LP;
deltag_LP = dWdg_an_LP - dWdg_num_LP;
deltah_LP = dWdh_an_LP - dWdh_num_LP;

delta_relative_SP = [deltaL/dWdL_num; deltaG/dWdG_num; deltaH/dWdH_num; deltal/dWdl_num; deltag/dWdg_num; deltah/dWdh_num];
delta_relative_LP = [deltaL_LP/dWdL_num_LP; deltaG_LP/dWdG_num_LP; deltaH_LP/dWdH_num_LP; 0; deltag_LP/dWdg_num_LP; deltah_LP/dWdh_num_LP];
%% Local Functions

function [dEdL, dEdG, dalphadG, dbetadG, dalphadH, dbetadH, dalphadg, dbetadg, dalphadh, dbetadh] = small_partials(a, r, e, l_, g, h, L, G, H, E)

dEdL = a/r*G^2/e/L^3*sin(E);
dEdG = -a/r*G/e/L^2*sin(E);

dalphadG = -H/G^2*sin(g)*sin(l_-h);
dbetadG = -H/G^2*cos(g)*sin(l_-h);
dalphadH = 1/G*sin(g)*sin(l_-h);
dbetadH = 1/G*cos(g)*sin(l_-h);

dalphadg = -sin(g)*cos(l_-h)+cos(g)*sin(l_-h)*H/G;
dbetadg = -cos(g)*cos(l_-h)-sin(g)*sin(l_-h)*H/G;
dalphadh = cos(g)*sin(l_-h)-sin(g)*cos(l_-h)*H/G;
dbetadh = -sin(g)*sin(l_-h)-cos(g)*cos(l_-h)*H/G;

end

%% Short Period Averaging
function dWdL = dWdL_SP(eps, a_, r_, mu, l_, g, h, L, G, H, E)
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);
a = L^2/mu;
e = sqrt(1-G^2/L^2);
r = a*(1-e*cos(E));

alpha = cos(g)*cos(l_-h) + sin(g)*sin(l_-h)*H/G;
beta = -sin(g)*cos(l_-h) + cos(g)*sin(l_-h)*H/G;

koeff = mu_*n_^2/(2*mu^4)*(a_/r_)^3*1/eps;

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

function dWdG = dWdG_SP(eps, a_, r_, mu, l_, g, h, L, G, H, E)
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);
a = L^2/mu;
e = sqrt(1-G^2/L^2);
r = a*(1-e*cos(E));

alpha = cos(g)*cos(l_-h) + sin(g)*sin(l_-h)*H/G;
beta = -sin(g)*cos(l_-h) + cos(g)*sin(l_-h)*H/G;
dalphadG = -H/G^2*sin(g)*sin(l_-h);
dbetadG = -H/G^2*cos(g)*sin(l_-h);

koeff = mu_*n_^2/(2*mu^4)*(a_/r_)^3*1/eps;

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

function dWdH = dWdH_SP(eps, a_, r_, mu, l_, g, h, L, G, H, E)
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);
e = sqrt(1-G^2/L^2);

alpha = cos(g)*cos(l_-h) + sin(g)*sin(l_-h)*H/G;
beta = -sin(g)*cos(l_-h) + cos(g)*sin(l_-h)*H/G;
dalphadH = 1/G*sin(g)*sin(l_-h);
dbetadH = 1/G*cos(g)*sin(l_-h);

koeff = mu_*n_^2/(2*mu^4)*(a_/r_)^3*1/eps;

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

function dWdl = dWdl_SP(eps, a_, r_, mu, l_, g, h, L, G, H, E)
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);
e = sqrt(1-G^2/L^2);

alpha = cos(g)*cos(l_-h) + sin(g)*sin(l_-h)*H/G;
beta = -sin(g)*cos(l_-h) + cos(g)*sin(l_-h)*H/G;

koeff = mu_*n_^2/(2*mu^4)*(a_/r_)^3*1/eps;

dWdE = koeff*L^7*(3*(1/12-1/4*beta^2)*e^3*cos(3*E) ...
    + (3*alpha^2-3/4*beta^2-3/4)*e^3*cos(E) ...
    + (3*alpha^2+3/2*beta^2-3/2)*e^2*cos(2*E) ...
    + 3*(e^2+1)*alpha*beta*G/L*sin(2*E) - 3/2*alpha*beta*e*(5*sin(E)+sin(3*E))*G/L ...
    + 1/4*(-27*alpha^2 + 3*beta^2+8)*cos(E)*e ...
    + 3/4*(beta^2-alpha^2)*cos(3*E)*e ...
    + 3/2*(alpha^2-beta^2)*cos(2*E));

dWdl = 1/(1-e*cos(E))*dWdE;
end

function dWdg = dWdg_SP(eps, a_, r_, mu, l_, g, h, L, G, H, E)
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);
e = sqrt(1-G^2/L^2);

alpha = cos(g)*cos(l_-h) + sin(g)*sin(l_-h)*H/G;
beta = -sin(g)*cos(l_-h) + cos(g)*sin(l_-h)*H/G;
dalphadg = -sin(g)*cos(l_-h)+cos(g)*sin(l_-h)*H/G;
dbetadg = -cos(g)*cos(l_-h)-sin(g)*sin(l_-h)*H/G;

koeff = mu_*n_^2/(2*mu^4)*(a_/r_)^3*1/eps;

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

function dWdh = dWdh_SP(eps, a_, r_, mu, l_, g, h, L, G, H, E)
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);
e = sqrt(1-G^2/L^2);

alpha = cos(g)*cos(l_-h) + sin(g)*sin(l_-h)*H/G;
beta = -sin(g)*cos(l_-h) + cos(g)*sin(l_-h)*H/G;
dalphadh = cos(g)*sin(l_-h)-sin(g)*cos(l_-h)*H/G;
dbetadh = -sin(g)*sin(l_-h)-cos(g)*cos(l_-h)*H/G;

koeff = mu_*n_^2/(2*mu^4)*(a_/r_)^3*1/eps;

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

%% Long Period Averaging
function dWdL = dWdL_LP(eps, a_, r_, mu, l_, g, h, L, G, H)
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);
e = sqrt(1-G^2/L^2);

koeff = mu_*n_/(2*mu^2)*1/eps;

dWdL = koeff*((a_/r_)^3*L/8*(3*G^2-10*L^2)*(2*l_*(1-3*H^2/G^2)+3*(H^2/G^2-1)*sin(2*l_-2*h)) ...
    + (a_/r_)^3*15/4*L/G^2*cos(2*g)*(l_-h)*(G^2-H^2)*(2*L^2-G^2) ...
    + (a_/r_)^3*15/8*L/G^2*(2*L^2-G^2)*cos(l_-h)*((G+H)^2*sin(l_-h-2*g)+(G-H)^2*sin(l_-h+2*g)) ...
    - l_/4*L/G^2*(G^2-3*H^2)*(3*G^2-10*L^2) ...
    +15/4*l_*L/G^2*(G^2-H^2)*(G^2-2*L^2)*cos(2*g));
end

function dWdG = dWdG_LP(eps, a_, r_, mu, l_, g, h, L, G, H)
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);
e = sqrt(1-G^2/L^2);

koeff = mu_*n_/(2*mu^2)*1/eps;

dWdG = koeff*((a_/r_)^3*3*L^2/8/G^3*(2*l_*(G^4-5*H^2*L^2)-(3*G^4-5*H^2*L^2)*sin(2*l_-2*h)) ...
    + (a_/r_)^3*15/4*L^2/G^3*cos(2*g)*(h-l_)*(G^4-H^2*L^2) ...
    - (a_/r_)^3*15/8*L^2/G*cos(l_-h)*((G+H)^2*sin(l_-h-2*g)+(G-H)^2*sin(l_-h+2*g)) ...
    - (a_/r_)^3*15/8*L^4/G^3*e^2*cos(l_-h)*((G+H)^2*sin(l_-h-2*g)+(G-H)^2*sin(l_-h+2*g)) ...
    + (a_/r_)^3*15/16*L^4/G^2*e^2*cos(l_-h)*(2*(G+H)*sin(l_-h-2*g)+2*(G-H)*sin(l_-h+2*g)) ...
    - 3/4*l_*L^2/G^3*(G^4-5*H^2*L^2) ...
    + 15/4*l_*L^2/G^3*(G^4-H^2*L^2)*cos(2*g));
end

function dWdH = dWdH_LP(eps, a_, r_, mu, l_, g, h, L, G, H)
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);
e = sqrt(1-G^2/L^2);

koeff = mu_*n_/(2*mu^2)*1/eps;

dWdH = koeff*H*L^4/G^2*((a_/r_)^3*3/8*(5-3*G^2/L^2)*(2*l_-sin(2*l_-2*h)) ...
    + (a_/r_)^3*15/4*e^2*cos(2*g)*(h-l_) ...
    + (a_/r_)^3*15/8/H*e^2*cos(l_-h)*((G+H)*sin(l_-h-2*g)-(G-H)*sin(l_-h+2*g)) ...
    + 3/4*l_*(3*G^2/L^2-5) ...
    + 15/4*l_*e^2*cos(2*g));
end

function dWdg = dWdg_LP(eps, a_, r_, mu, l_, g, h, L, G, H)
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);
e = sqrt(1-G^2/L^2);

koeff = mu_*n_/(2*mu^2)*1/eps;

dWdg = koeff*L^4*e^2*((a_/r_)^3*15/4*sin(2*g)*(h-l_)*(1-H^2/G^2) ...
    - (a_/r_)^3*15/8/G^2*cos(l_-h)*((G+H)^2*cos(l_-h-2*g)-(G-H)^2*cos(l_-h+2*g)) ...
    + l_*15/4*(1-H^2/G^2)*sin(2*g));
end

function dWdh = dWdh_LP(eps, a_, r_, mu, l_, g, h, L, G, H)
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);
e = sqrt(1-G^2/L^2);

koeff = mu_*n_/(2*mu^2)*1/eps;

dWdh = koeff*(a_/r_)^3*3/8*L^4*((H^2/G^2-1)*(5-3*G^2/L^2)*cos(2*l_-2*h) ...
    - 5*e^2*cos(2*g)*(1-H^2/G^2) ...
    + 5/2/G^2*e^2*(sin(l_-h)*((G+H)^2*sin(l_-h-2*g)+(G-H)^2*sin(l_-h+2*g))-cos(l_-h)*((G+H)^2*cos(l_-h-2*g)+(G-H)^2*cos(l_-h+2*g))));
end