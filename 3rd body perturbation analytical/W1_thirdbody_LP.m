function [W1] = W1_thirdbody_LP(eps, a_, r_, mu, l_, g, h, L, G, H)
%W1_THIRDBODY_LP calculates W1 for removing all short period terms (dependent on l_)
%from a third body perturbation
global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;

n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %only works for nondimensionalized equations!
mu_ = Mass_Earth/(Mass_Moon + Mass_Earth);

koeff = mu_*n_/(2*mu^2)*1/eps;

e = sqrt(1-G^2/L^2);

W1 = koeff*(L^4*(a_/r_)^3*(1+3/2*e^2)*(3/4*H^2/G^2*(l_-1/2*sin(2*l_-2*h))+3/8*sin(2*l_-2*h)-l_/4) ...
    + L^4*(a_/r_)^3*15/4*e^2*1/2*(1-H^2/G^2)*(l_-h)*cos(2*g) ...
    + L^4*(a_/r_)^3*15/4*e^2*1/(4*G^2)*cos(l_-h)*((G+H)^2*sin(l_-h-2*g)+(G-H)^2*sin(l_-h+2*g)) ...
    - L^4*l_/8*(3*H^2/G^2-1)*(2+3*e^2) ...
    - L^4*l_/8*15*(1-H^2/G^2)*e^2*cos(2*g));

end

