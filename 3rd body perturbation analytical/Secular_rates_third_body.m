function [derivative] = Secular_rates_third_body(t, Delaunay)
%SECULAR_RATES_THIRD_BODY Calculates derivatives on the double averaged
%Hamiltonian for secular Delaunay elements rates

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

e = sqrt(1-(G/L)^2);
inc = acos(H/G);

%The derivatives of the double averaged R4 are caclulated below
K2 = 9*mu_*n_^2/(65.536*a_^2*mu^2);
C1 = 144+320*cos(2*inc)+560*cos(4*inc);
C3 = 1680+2240*cos(2*inc)-3920*cos(4*inc);
C6 = 4410-5880*cos(2*inc)+1470*cos(4*inc);

dR4dL_part = K2*4*L^3*(C1+5*C1*e^2+C3*e^2*cos(2*g)+15/8*C1*e^4+1/2*C3*e^4*cos(2*g)+C6*e^4*cos(4*g));
dR4de = K2*L^4*(10*C1*e+2*C3*e*cos(2*g)+15/2*C1*e^3+2*C3*e^3*cos(2*g)+4*C6*e^3*cos(4*g));
dedL = G^2/e/L^3;
dedG = -G/e/L^2;
dR4dC1 = K2*L^4*(1+5*e^2+15/8*e^4);
dR4dC3 = K2*L^4*(e^2*cos(2*g)+1/2*e^4*cos(2*g));
dR4dC6 = K2*L^4*(e^4*cos(4*g));
dC1di = -640*sin(2*inc)-2240*sin(4*inc);
dC3di = -4480*sin(2*inc)+15680*sin(4*inc);
dC6di = 11760*sin(2*inc)-5880*sin(4*inc);
didG = H/G^2/sin(inc);
didH = -1/G/sin(inc);

dR4dg = 0;%K2*L^4*(-2*C3*e^2*sin(2*g)-C3*e^4*sin(2*g)-4*C6*e^4*sin(4*g));
dR4dL = 0;%dR4dL_part + dR4de*dedL;
dR4dG = 0;%dR4de*dedG + dR4dC1*dC1di*didG + dR4dC3*dC3di*didG + dR4dC6*dC6di*didG;
dR4dH = 0;%dR4dC1*dC1di*didH + dR4dC3*dC3di*didH + dR4dC6*dC6di*didH;
%
 
% The sign are changed from my paper solution and also in the osc2mean
% transformation and back as well
%the first term here has to be added/removed, depends!!! mu^2/L^3
derivative(1,1) = + mu_*n_^2/8/mu^2*L*((3*H^2/G^2-1)*(3*G^2-10*L^2)-15*(H^2/G^2-1)*(G^2-2*L^2)*cos(2*g)) - dR4dL;
%derivative(1,1) = 0.072671835808406;
derivative(2,1) = -3*mu_*n_^2/8/mu^2*L^2*(5*H^2*L^2/G^3*(cos(2*g)-1)-5*G*cos(2*g)+G) - dR4dG;
%derivative(2,1) = -1.403400445022995e-04;
derivative(3,1) = -3*mu_*n_^2/8/mu^2*H*L^2*(5*L^2/G^2*(1-cos(2*g))+5*cos(2*g)-3) - dR4dH;
%derivative(3,1) = 6.712655014468600e-05;
derivative(4,1) = 0;
derivative(5,1) = -15*mu_*n_^2/8/mu^2*L^2*(H^2/G^2-1)*(G^2-L^2)*sin(2*g) + dR4dg;
%derivative(5,1) = 8.932042992936895e-07;
derivative(6,1) = 0;
end

