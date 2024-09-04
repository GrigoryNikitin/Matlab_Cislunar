%Here are the plots of differences in frozen orbit inclinations are
%calculated

clear;
clc;

Mass_Earth = 5.9722*10^24;
Mass_Moon = 7.3477*10^22;
mu_Moon = 4902.8001224453001;
G = 6.674*10^-20;
mu_Earth = G*Mass_Earth;

% Units for nondimensionalization 
runit = 1738.0;
tunit = sqrt(1738.0^3/(mu_Moon));
vunit = runit/tunit;

%Coefficients from Bharat's paper Table 1
J_2 = 0.000203223560871;
J_4 = -9.70434253277103*10^-6;
J_6 = -1.37675387844709*10^-5;

%normilized r_Moon and mu_Moon
re = 1;
R_e = 1;
r_m = 1;
mu = 1;
mu_ = mu_Earth/(mu_Moon+mu_Earth);
a_ = 384748/runit;
n_ = sqrt((mu+mu_Earth/mu_Moon)/a_^3);

% SET OMEGA HERE
omega = 90*pi/180;

%dimensions of the matrix of differences
dimensions = 101;
ecc = linspace(0,0.9,dimensions);
a_dim = linspace(2000,10000,dimensions);

for i = 1:dimensions
    e = ecc(i);
    for j = 1:dimensions
        a = a_dim(j)/runit;
        n = sqrt(mu/a^3);
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
        for l = 1:length(solution)
            if solution(l) > 0 && solution(l) < 1
                k = k+1;
                m = l;
            end
        end
        if k == 0
            disp('ERROR: zero solutions for inclination are found');
            inc = NaN;
        elseif k == 1
            disp('A solution for inclination has been found');
            disp([i j]);
            inc = double(acos(sqrt(solution(m))));
        else
            disp('ERROR: more than 1 solution for inclination has been found');
            %just taking a smaller one?
            %always the second?
            inc = double(acos(sqrt(solution(2))));
            %inc = NaN;
        end
        
        %Different equations for diff omega!
        %comparison with Nie's equation (60):
        k_M = Mass_Earth/(Mass_Earth+Mass_Moon);
        buba = -(k_M*n_^2*a^5*(1-e^2)^1.5*(6*e^2-1)+3*J_2*mu*r_m^2) / (5*(a^5*(1-e^2)^1.5*k_M*n_^2+J_2*mu*r_m^2));
        inc_Nie = 1/2*acos(buba);
        
        difference_arr(i,j) = inc - inc_Nie;
    end
end
difference_arr_deg = difference_arr.*180/pi;

%% Plots

%difference_arr_deg_flip = flip(difference_arr_deg,1); %reverted the columns
imagesc([2000 10000], [0 0.9], difference_arr_deg);
set(gca, 'YDir', 'normal'); % Correct the Y-axis direction
hold on;
xlabel('Mean a, km','FontWeight','bold');
ylabel('Mean ecc','FontWeight','bold');

% Set the colormap and color limits
colormap;
caxis([min(difference_arr_deg(:)), max(difference_arr_deg(:))]);

cb = colorbar;
cb.Label.String = 'Difference in inc, deg';
cb.Label.FontWeight = 'bold';
cb.Direction = 'reverse';

% Define the x values for plotting the equation
x = linspace(2000, 10000, 1000); % Generate 1000 points between 2000 and 10000
y = 1 - 1738 ./ x;

% Plot the line
plot(x, y, 'r--', 'LineWidth', 2); % Plot a red line with a specified width

hold off;


% Calculate the indices of the matrix data that correspond to the plot coordinates
xImage = linspace(2000, 10000, size(difference_arr_deg, 2));
yImage = linspace(0, 0.9, size(difference_arr_deg, 1));
[XImage, YImage] = meshgrid(xImage, yImage);

% Determine which matrix points are above the line
aboveLine = YImage > (1 - 1738 ./ XImage);

difference_arr_deg_new = difference_arr_deg;

% Set the values above the line to 0 in the data matrix
difference_arr_deg_new(aboveLine) = 0.1;

% % Create a custom blue-to-green-to-yellow-to-red colormap
% nColors = 1024; % Number of colors in the colormap
% % Define the transition points
% nBlueToGreen = floor(nColors / 3);
% nGreenToYellow = floor(nColors / 3);
% nYellowToRed = nColors - nBlueToGreen - nGreenToYellow;
% % Create the blue to green gradient
% blueToGreen = [linspace(0, 1, nBlueToGreen)', linspace(0, 0, nBlueToGreen)', linspace(1, 0, nBlueToGreen)'];
% % Create the green to yellow gradient
% greenToYellow = [ones(nGreenToYellow, 1), linspace(0, 1, nGreenToYellow)', zeros(nGreenToYellow, 1)];
% % Create the yellow to red gradient
% yellowToRed = [linspace(1, 0, nYellowToRed)', ones(nYellowToRed, 1), zeros(nYellowToRed, 1)];
% % Combine the gradients
% blueToRed = [blueToGreen; greenToYellow; yellowToRed];
% 
% % Add white to the end of the colormap
% customColormap = [blueToRed; 1 1 1];

colormap(parula); % Use the default parula colormap
customColormap = [colormap; 1 1 1];

figure;
imagesc([2000 10000], [0 0.9], difference_arr_deg_new);
set(gca, 'YDir', 'normal'); % Correct the Y-axis direction
hold on;
xlabel('Mean a, km','FontWeight','bold');
ylabel('Mean ecc','FontWeight','bold');

% Set the colormap and color limits
colormap(customColormap);
caxis([min(difference_arr_deg_new(:)), max(difference_arr_deg_new(:))]);

cb = colorbar;
cb.Label.String = 'Difference in inc, deg';
cb.Label.FontWeight = 'bold';
cb.Direction = 'reverse';

% Plot the line
plot(x, y, 'r--', 'LineWidth', 2); % Plot a red line with a specified width

hold off;
