function [Kepler_elements] = Dec2Kepl(mu, sv)
%DEC2KEPL returns vector of [a,e,i,omega,RAAN,MA]

X = sv(1);
Y = sv(2);
Z = sv(3);
Vx = sv(4);
Vy = sv(5);
Vz = sv(6);

r = sqrt(X^2 + Y^2 + Z^2);

h_x = Y * Vz - Z * Vy;
h_y = -X * Vz + Z * Vx;
h_z = X * Vy - Y * Vx;
h = sqrt(h_x^2 + h_y^2 + h_z^2);
i = acos(h_z / h);                     %inc КА
if i == 0
    i = 0.0000000001;
end
if i < 0.000000001
    RAAN = 0.0000000001;               %RAAN КА
else
    RAAN = atan2(h_x, -h_y);
end
if RAAN < 0
    RAAN = RAAN + 2 * pi;
end

Lap_x = -mu * X / r - Vz * (Z * Vx - X * Vz) + Vy * (X * Vy - Y * Vx);
Lap_y = -mu * Y / r - Vx * (X * Vy - Y * Vx) + Vz * (Y * Vz - Z * Vy);
Lap_z = -mu * Z / r - Vy * (Y * Vz - Z * Vy) + Vx * (Z * Vx - X * Vz);
Lap = sqrt(Lap_x^2 + Lap_y^2 + Lap_z^2);

p = h * h / mu;                                       %параметр орбиты
e = Lap / mu;                                        %Эксцентриситет орбиты КА
a = p / (1 - e^2);                   %Большая полуось орбиты КА

if e <= 10^-7
    omega = 0;
else
    omega = atan2(Lap_z / (Lap * sin(i)), (Lap_x * cos(RAAN) + Lap_y * sin(RAAN)) / Lap); %Аргумент перицентра КА
end
if omega < 0
    omega = omega + 2 * pi;
end
        
u = atan2(Z / (r * sin(i)), (X * cos(RAAN) + Y * sin(RAAN)) / r);   %Аргумент широты КА
if u < 0
    u = u + 2 * pi;
end
    
if e <= 10^-7
    f = u - omega; %Истинная аномалия
    if f < 0
        f = f + 2 * pi;
    end
else
    f = atan2((X * Vx + Y * Vy + Z * Vz) / (r * sqrt(mu/p)*e), (p / r - 1)/e);
    if abs(f) < 10^-10
        f = 0;
    end
    if f < 0
        f = f + 2 * pi;
    end
end
                
if e < 1
    EA = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(0.5 * f));
    if abs(EA) < 10^-9
        EA = 0;
    end
    if EA < 0                                           %EccAnomaly and MeanAnomaly
        EA = EA + 2 * pi;
    end
    MA = EA - e * sin(EA);
else
    HA = atanh(sqrt(e^2 - 1) * sin(f) / (e + cos(f))); %HypAnomaly and MeanAnomaly
    MA = e * sinh(HA) - HA;
end
            
%tau = self.time - sqrt(abs(self.a * self.a * self.a) / mu) * self.MA;            %Время прохождения перицентра

%Kepler_elements = [a,e,i,omega,RAAN,f,MA];
Kepler_elements = [a,e,i,omega,RAAN,MA]; % deleted TA from the output
end

