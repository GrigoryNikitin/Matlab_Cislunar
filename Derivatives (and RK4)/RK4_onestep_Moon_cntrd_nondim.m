function [t, y] = RK4_onestep_Moon_cntrd_nondim(t, y, h, M0_moon, Degree, Order, Clm, Slm)
%RK4_ONESTEP_MOON_CNTRD integrates nondimensionalized "derivatives" for 1 time step

k1 = derivatives_Moon_cntrd_nondim(t,y, M0_moon, Degree, Order, Clm, Slm);
k2 = derivatives_Moon_cntrd_nondim(t+h/2, y+k1'*h/2, M0_moon, Degree, Order, Clm, Slm);
k3 = derivatives_Moon_cntrd_nondim(t+h/2, y+k2'*h/2, M0_moon, Degree, Order, Clm, Slm);
k4 = derivatives_Moon_cntrd_nondim(t+h, y+k3'*h, M0_moon, Degree, Order, Clm, Slm);

y = y' + (k1 + 2*k2 + 2*k3 + k4) * h / 6;
t = t + h;
end

