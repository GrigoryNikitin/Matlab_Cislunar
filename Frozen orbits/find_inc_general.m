function inc = find_inc_general(omega,e,a,a_,mu,mu_,R_e)
%This function finds solutions for inclination, where domega/dt = 0
% For any omega!
% Works for specific omega, e, a, mu and mu_

global mu_Earth;
global mu_Moon;
global Mass_Earth;
global Mass_Moon;
%Coefficients from Bharat's paper Table 1
%I forgot the sign change here too!
J_2 = 0.000203223560871;
J_4 = -9.70434253277103*10^-6;
J_6 = -1.37675387844709*10^-5;

n = sqrt(mu/a^3);
n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3); %that's should be correct
%n_ = 2.66e-06; %this is dimensionalised 1/sec


%SITUATION NOW:
%I changed the signs back (as in Overleaf) in domega/dt ONLY for the Grav
% part. The 3rd body part is still with another sign (-1*Overleaf).
%Changed the signs in J2,J4,J6 coefficients as they were wrong.

syms t
eqn = - 3*mu_*n_^2/(8*n)*(5*(cos(2*omega)-1)/sqrt(1-e^2)*t - 5*sqrt(1-e^2)*cos(2*omega)+sqrt(1-e^2)) ...
    - J_2*3*n*R_e^2/(4*a^2*(1-e^2)^2)*(1-5*t) ...
    + J_4*15*n*R_e^4/(128*a^4*(1-e^2)^4)*(9*(1-e^2)-21*(6*(1-e^2)*t+1)+27*t*(7*(1-e^2)*t+10)-385*t^2) ...
    + J_6*5*n*R_e^6/(2048*(1-e^2)^6*a^6)*((5-105*t+315*t^2-231*t^3)*(60*(1-e^2)^2-140*(1-e^2)) ...
    + (210*t-1260*t^2+1386*t^3)*(63+15*(1-e^2)^2-70*(1-e^2)) ...
    - 11*(5-105*t+315*t^2-231*t^3)*(63+15*(1-e^2)^2-70*(1-e^2))) == 0;

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
    %just taking a smaller one?
    %always the second?
    inc = double(acos(sqrt(solution(2))));
    %inc = NaN;
end

%checking that the derivative is zero:
t = cos(inc)^2;
check = - 3*mu_*n_^2/(8*n)*(5*(cos(2*omega)-1)/sqrt(1-e^2)*t - 5*sqrt(1-e^2)*cos(2*omega)+sqrt(1-e^2)) ...
    - J_2*3*n*R_e^2/(4*a^2*(1-e^2)^2)*(1-5*t) ...
    + J_4*15*n*R_e^4/(128*a^4*(1-e^2)^4)*(9*(1-e^2)-21*(6*(1-e^2)*t+1)+27*t*(7*(1-e^2)*t+10)-385*t^2) ...
    + J_6*5*n*R_e^6/(2048*(1-e^2)^6*a^6)*((5-105*t+315*t^2-231*t^3)*(60*(1-e^2)^2-140*(1-e^2)) ...
    + (210*t-1260*t^2+1386*t^3)*(63+15*(1-e^2)^2-70*(1-e^2)) ...
    - 11*(5-105*t+315*t^2-231*t^3)*(63+15*(1-e^2)^2-70*(1-e^2)));
check_1 = - 3*mu_*n_^2/(8*n)*(5*(cos(2*omega)-1)/sqrt(1-e^2)*t - 5*sqrt(1-e^2)*cos(2*omega)+sqrt(1-e^2));
check_2 = - J_2*3*n*R_e^2/(4*a^2*(1-e^2)^2)*(1-5*t) ...
    + J_4*15*n*R_e^4/(128*a^4*(1-e^2)^4)*(9*(1-e^2)-21*(6*(1-e^2)*t+1)+27*t*(7*(1-e^2)*t+10)-385*t^2) ...
    + J_6*5*n*R_e^6/(2048*(1-e^2)^6*a^6)*((5-105*t+315*t^2-231*t^3)*(60*(1-e^2)^2-140*(1-e^2)) ...
    + (210*t-1260*t^2+1386*t^3)*(63+15*(1-e^2)^2-70*(1-e^2)) ...
    - 11*(5-105*t+315*t^2-231*t^3)*(63+15*(1-e^2)^2-70*(1-e^2)));

%Their J2 is negative???
J_2 = 0.000203223560871;

%comparison with Y. Matsumoto's equation (25):
%k_M = mu_/(mu+mu_);
k_M = Mass_Earth/(Mass_Earth+Mass_Moon);
%k_M = mu_;
r_m = 1;
A_M = 3/8*k_M*n_^2/n;
B_M = 3*n*J_2*r_m^2/a^2;
inc_Matsumoto = asin(sqrt( (8*(2+3*e^2)*(1-e^2)^1.5*A_M+4*B_M)/(40*A_M*(1-e^2)^1.5+5*B_M) ));

%comparison with Nie's equation (60):
%k_M = mu_/(mu+mu_);
k_M = Mass_Earth/(Mass_Earth+Mass_Moon);
buba = -(k_M*n_^2*a^5*(1-e^2)^1.5*(6*e^2-1)+3*J_2*mu*r_m^2) / (5*(a^5*(1-e^2)^1.5*k_M*n_^2+J_2*mu*r_m^2));
inc_Nie = 1/2*acos(buba);

%inc = inc_Nie;

end

