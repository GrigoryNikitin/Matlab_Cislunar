function [sv] = Kepl2Dec(mu,Kepler_elements)
%KEPL2DEC takes [a,e,i,omega,RAAN,MA] and returns vector [XYZVxVyVz]

a = Kepler_elements(1);
e = Kepler_elements(2);
inc = Kepler_elements(3);
omega = Kepler_elements(4);
RAAN = Kepler_elements(5);
MA = Kepler_elements(6);

EA = Kepler_Eqn_solver(MA, e, 10e-9);
TA = atan2(sqrt(1 - e^2) * sin(EA) / (1 - e * cos(EA)), (cos(EA) - e) / (1 - e * cos(EA)));
if (TA < 0)
    TA = TA + 2 * pi;
end

p = a * (1 - e * e);
r = p / (1 + e * cos(TA));
V_r = sqrt(mu / p) * e * sin(TA);
V_m = sqrt(mu / p) * (1 + e * cos(TA));
u = omega + TA;
X = (cos(RAAN) * cos(u) - sin(RAAN) * sin(u) * cos(inc)) * r;
Y = (sin(RAAN) * cos(u) + cos(RAAN) * sin(u) * cos(inc)) * r;
Z = (sin(u) * sin(inc)) * r;
Vx = (cos(RAAN) * cos(u) - sin(RAAN) * sin(u) * cos(inc)) * V_r + (-cos(RAAN) * sin(u) - sin(RAAN) * cos(u) * cos(inc)) * V_m;
Vy = (sin(RAAN) * cos(u) + cos(RAAN) * sin(u) * cos(inc)) * V_r + (-sin(RAAN) * sin(u) + cos(RAAN) * cos(u) * cos(inc)) * V_m;
Vz = (sin(u) * sin(inc)) * V_r + (cos(u) * sin(inc)) * V_m;

sv = [X, Y, Z, Vx, Vy, Vz];

end

