function [t, y] = RK4_onestep(t, y, h, M0_moon, include_Earth)
%RK4_ONESTEP integrates "derivatives" for 1 time step

if include_Earth == 1
    k1 = derivatives_cm(t,y, M0_moon);
    k2 = derivatives_cm(t+h/2, y'+k1*h/2, M0_moon);
    k3 = derivatives_cm(t+h/2, y'+k2*h/2, M0_moon);
    k4 = derivatives_cm(t+h, y'+k3*h, M0_moon);
else
    k1 = derivatives_moon(t,y);
    k2 = derivatives_moon(t+h/2, y'+k1*h/2);
    k3 = derivatives_moon(t+h/2, y'+k2*h/2);
    k4 = derivatives_moon(t+h, y'+k3*h);
end

y = y' + (k1 + 2*k2 + 2*k3 + k4) * h / 6;
t = t + h;

end

