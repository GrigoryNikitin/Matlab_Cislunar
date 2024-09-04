function [W1] = W1_thirdbody_SP_new_one(eps, a_, r_, mu, l_, l, g, h, L, G, H)
%W1_THIRDBODY_SP calculates W1 for removing all short period terms (dependent on l)
%from a third body perturbation
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);

koeff = mu_*n_^2*L^7/(2*mu^4)*(a_/r_)^3*1/eps;
alpha = cos(g)*cos(l_-h) + sin(g)*sin(l_-h)*H/G;
beta = -sin(g)*cos(l_-h) + cos(g)*sin(l_-h)*H/G;

e = sqrt(1-G^2/L^2);
E = Kepler_Eqn_solver(l, e, 10^-13);

W1 = koeff*(alpha^2*(3/2*E+6*E*e^2-3*(e^2+11/4)*e*sin(E)+(3/2*e^2+3/4)*sin(2*E)-e/4*sin(3*E)) ...
    + alpha*beta*sqrt(1-e^2)*cos(E)*(e*(cos(2*E)+7)-3*(e^2+1)*cos(E)) ...
    + 1/4*(1-e^2)*beta^2*(6*E-3*e*sin(E)-3*sin(2*E)+e*sin(3*E)) ...
    -(1+3/2*e^2)*E+3/4*e*(4+e^2)*sin(E)-3/4*e^2*sin(2*E)+1/12*e^3*sin(3*E) ...
    -((1+3/2*e^2)*(3/2*(alpha^2+beta^2)-1)+15/4*e^2*(alpha^2-beta^2))*l);

end

