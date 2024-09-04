function inc = find_inc(e,a,a_,mu,mu_,R_e)
%This function finds solutions for inclination, where domega/dt = 0
% All equations here assume that omega = 90 deg!
% Works for specific e, a, mu and mu_

%Coefficients from Bharat's paper Table 1
J_2 = -0.000203223560871;
J_4 = 9.70434253277103*10^-6;
J_6 = 1.37675387844709*10^-5;

n = sqrt(mu/a^3);
n_ = sqrt(mu/a_^3);
gamma_1 = 3*mu_*n_^2/(8*n)*sqrt(1-e^2);
gamma_2 = 3*mu*R_e^2/(4*n*a^5*(1-e^2)^2);
J_2_star = J_2;
J_4_star = J_4*R_e^2/(16*a^2);
J_6_star = J_6*5*R_e^4/(512*a^4);

%Coefficients of the equation A+Bt+Ct^2+Dt^3 = 0:
A = -4*gamma_1 + gamma_2*J_2_star + 5/2*gamma_2*J_4_star*(12+9*e^2)/(1-e^2)^2 + 5*gamma_2*J_6_star*(56+140*e^2+35*e^4)/(1-e^2)^4;
B = 5*(2*gamma_1/(1-e^2)-gamma_2*J_2_star-gamma_2*J_4_star*(72+27*e^2)/(1-e^2)^2-gamma_2*J_6_star*(1288+3500*e^2+945*e^4)/(1-e^2)^4);
C = 5*(gamma_2*J_4_star*(84+21*e^2)/(2*(1-e^2)^2) + gamma_2*J_6_star*(4200+12180*e^2+3465*e^4)/(1-e^2)^4);
D = -gamma_2*J_6_star*(16632+50820*e^2+15015*e^4)/(1-e^2)^4;

syms t
eqn = A+B*t+C*t^2+D*t^3 == 0;
solution = vpasolve(eqn,t);

%checking how many solutions in range [0,1]
k = 0;
for i = 1:length(solution)
    if solution(i) > 0 && solution(i) < 1
        k = k+1;
        j = i;
    end
end
if k == 0
    disp('ERROR: zero solutions for inclination are found');
    inc = NaN;
elseif k == 1
    disp('A solution for inclination has been found');
    inc = double(acos(sqrt(solution(j))));
else
    disp('ERROR: more than 1 solution for inclination has been found');
    inc = NaN;
end


%comparison with Y. Matsumoto's equation (25):
k_M = mu_/(mu+mu_);
%k_M = mu_;
r_m = 1;
A_M = 3/8*k_M*n_^2/n;
B_M = 3*n*J_2*r_m^2/a^2;
inc_Matsumoto = asin(sqrt( (8*(2+3*e^2)*(1-e^2)^1.5*A_M+4*B_M)/(40*A_M*(1-e^2)^1.5+5*B_M) ));

%comparison with Nie's equation (60):
k_M = mu_/(mu+mu_);
buba = -(k_M*n_^2*a^5*(1-e^2)^1.5*(6*e^2-1)+3*J_2*mu*r_m) / (5*(a^5*(1-e^2)^1.5*k_M*n_^3+J_2*mu*r_m));
inc_Nie = 1/2*acos(buba);

%inc = inc_Nie;

end

