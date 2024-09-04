function EA = Kepler_Eqn_solver(Mean_Anomaly, ecc, tolerance)
%KEPLER_EQN_SOLVER solves Kepler's equation for E

EA_next = Mean_Anomaly + ecc * sin(Mean_Anomaly); 
k = 0;
while 1
    k = k+1;
    EA_ = EA_next;
    EA_next = EA_ - (EA_ - ecc * sin(EA_) - Mean_Anomaly) / (1 - ecc * cos(EA_));
    if abs(EA_next - EA_) < tolerance || k >= 1000
        break
    end
end
EA = EA_next;
end

