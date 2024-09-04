%% Creating new arrays to get info from files
array_of_r_0 = array_of_r;
Error_Kepler_of_r_0 = Error_Kepler_of_r;
Error_r_of_r_0 = Error_r_of_r;
Error_v_of_r_0 = Error_v_of_r;
LP_of_r_0 = LP_of_r;

%% PLOTS
figure;
tiledlayout(3,2);
nexttile;
plot(array_of_r, LP_of_r_0(:,1)*180/pi, '-', 'LineWidth', 2);
xlabel('a_0, km','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('l long period','FontWeight','bold');
grid on;
hold on;
plot(array_of_r, LP_of_r_30(:,1)*180/pi, '-', 'LineWidth', 2);
plot(array_of_r, LP_of_r_70(:,1)*180/pi, '-', 'LineWidth', 2);
plot(array_of_r, LP_of_r_90(:,1)*180/pi, '-', 'LineWidth', 2);
hold off;
legend('i=1 deg', 'i=30 deg', 'i=70 deg', 'i=90 deg');
xlim([1800 6000]);

nexttile;
plot(array_of_r, LP_of_r_0(:,2)*180/pi, '-', 'LineWidth', 2);
xlabel('a_0, km','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('g long period','FontWeight','bold');
grid on;
hold on;
plot(array_of_r, LP_of_r_30(:,2)*180/pi, '-', 'LineWidth', 2);
plot(array_of_r, LP_of_r_70(:,2)*180/pi, '-', 'LineWidth', 2);
plot(array_of_r, LP_of_r_90(:,2)*180/pi, '-', 'LineWidth', 2);
hold off;
legend('i=1 deg', 'i=30 deg', 'i=70 deg', 'i=90 deg', 'Location','southeast');
xlim([1800 6000]);

nexttile;
plot(array_of_r, LP_of_r_0(:,3)*180/pi, '-', 'LineWidth', 2);
xlabel('a_0, km','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('h long period','FontWeight','bold');
grid on;
hold on;
plot(array_of_r, LP_of_r_30(:,3)*180/pi, '-', 'LineWidth', 2);
plot(array_of_r, LP_of_r_70(:,3)*180/pi, '-', 'LineWidth', 2);
plot(array_of_r, LP_of_r_90(:,3)*180/pi, '-', 'LineWidth', 2);
hold off;
legend('i=1 deg', 'i=30 deg', 'i=70 deg', 'i=90 deg', 'Location','southeast');
xlim([1800 6000]);

nexttile;
plot(array_of_r, LP_of_r_0(:,4), '-', 'LineWidth', 2);
xlabel('a_0, km','FontWeight','bold');
ylabel('[-]','FontWeight','bold');
title('L long period','FontWeight','bold');
grid on;
hold on;
plot(array_of_r, LP_of_r_30(:,4), '-', 'LineWidth', 2);
plot(array_of_r, LP_of_r_70(:,4), '-', 'LineWidth', 2);
plot(array_of_r, LP_of_r_90(:,4), '-', 'LineWidth', 2);
hold off;
legend('i=1 deg', 'i=30 deg', 'i=70 deg', 'i=90 deg');
xlim([1800 6000]);

nexttile;
plot(array_of_r, LP_of_r_0(:,5), '-', 'LineWidth', 2);
xlabel('a_0, km','FontWeight','bold');
ylabel('[-]','FontWeight','bold');
title('G long period','FontWeight','bold');
grid on;
hold on;
plot(array_of_r, LP_of_r_30(:,5), '-', 'LineWidth', 2);
plot(array_of_r, LP_of_r_70(:,5), '-', 'LineWidth', 2);
plot(array_of_r, LP_of_r_90(:,5), '-', 'LineWidth', 2);
hold off;
legend('i=1 deg', 'i=30 deg', 'i=70 deg', 'i=90 deg');
xlim([1800 6000]);

nexttile;
plot(array_of_r, LP_of_r_0(:,6), '-', 'LineWidth', 2);
xlabel('a_0, km','FontWeight','bold');
ylabel('[-]','FontWeight','bold');
title('H long period','FontWeight','bold');
grid on;
hold on;
plot(array_of_r, LP_of_r_30(:,6), '-', 'LineWidth', 2);
plot(array_of_r, LP_of_r_70(:,6), '-', 'LineWidth', 2);
plot(array_of_r, LP_of_r_90(:,6), '-', 'LineWidth', 2);
hold off;
legend('i=1 deg', 'i=30 deg', 'i=70 deg', 'i=90 deg');
xlim([1800 6000]);


%%
figure;
tiledlayout(2,1)

nexttile
plot(array_of_r, Error_r_of_r_0, '-', 'LineWidth', 2);
xlabel('a_0, km','FontWeight','bold');
ylabel('[km]','FontWeight','bold');
title('\Delta{\itr}');
grid on;
hold on;
plot(array_of_r, Error_r_of_r_30, '-', 'LineWidth', 2);
plot(array_of_r, Error_r_of_r_70, '-', 'LineWidth', 2);
plot(array_of_r, Error_r_of_r_90, '-', 'LineWidth', 2);
hold off;
legend('i=1 deg', 'i=30 deg', 'i=70 deg', 'i=90 deg', 'Location','northeast');
xlim([1800 3000]);

nexttile
plot(array_of_r, Error_v_of_r_0, '-', 'LineWidth', 2);
xlabel('a_0, km','FontWeight','bold');
ylabel('[km/s]','FontWeight','bold');
title('\Delta{\itv}');
grid on;
hold on;
plot(array_of_r, Error_v_of_r_30, '-', 'LineWidth', 2);
plot(array_of_r, Error_v_of_r_70, '-', 'LineWidth', 2);
plot(array_of_r, Error_v_of_r_90, '-', 'LineWidth', 2);
hold off;
legend('i=1 deg', 'i=30 deg', 'i=70 deg', 'i=90 deg', 'Location','northeast');
xlim([1800 3000]);


%%
figure;
tiledlayout(2,3)

nexttile
plot(array_of_r, Error_Kepler_of_r_0(:,1), '-', 'LineWidth', 2);
xlabel('a_0, km','FontWeight','bold');
ylabel('[km]','FontWeight','bold');
title('\Delta{\ita}');
grid on;
hold on;
plot(array_of_r, Error_Kepler_of_r_30(:,1), '-', 'LineWidth', 2);
plot(array_of_r, Error_Kepler_of_r_70(:,1), '-', 'LineWidth', 2);
plot(array_of_r, Error_Kepler_of_r_90(:,1), '-', 'LineWidth', 2);
%plot(array_of_r, Error_Kepler_of_r_40(:,1), '-', 'LineWidth', 2);
hold off;
legend('i=1 deg', 'i=30 deg', 'i=70 deg', 'i=90 deg', 'i=40 deg');
xlim([1800 3000]);

nexttile
plot(array_of_r, Error_Kepler_of_r_0(:,2), '-', 'LineWidth', 2);
xlabel('a_0, km','FontWeight','bold');
ylabel('[-]','FontWeight','bold');
title('\Delta{\ite}');
grid on;
hold on;
plot(array_of_r, Error_Kepler_of_r_30(:,2), '-', 'LineWidth', 2);
plot(array_of_r, Error_Kepler_of_r_70(:,2), '-', 'LineWidth', 2);
plot(array_of_r, Error_Kepler_of_r_90(:,2), '-', 'LineWidth', 2);
%plot(array_of_r, Error_Kepler_of_r_40(:,2), '-', 'LineWidth', 2);
hold off;
legend('i=1 deg', 'i=30 deg', 'i=70 deg', 'i=90 deg', 'i=40 deg', 'Location','southeast');
xlim([1800 3000]);

nexttile
plot(array_of_r, rad2deg(Error_Kepler_of_r_0(:,3)), '-', 'LineWidth', 2);
xlabel('a_0, km','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('\Delta{\iti}');
grid on;
hold on;
plot(array_of_r, rad2deg(Error_Kepler_of_r_30(:,3)), '-', 'LineWidth', 2);
plot(array_of_r, rad2deg(Error_Kepler_of_r_70(:,3)), '-', 'LineWidth', 2);
plot(array_of_r, rad2deg(Error_Kepler_of_r_90(:,3)), '-', 'LineWidth', 2);
%plot(array_of_r, rad2deg(Error_Kepler_of_r_40(:,3)), '-', 'LineWidth', 2);
hold off;
legend('i=1 deg', 'i=30 deg', 'i=70 deg', 'i=90 deg', 'i=40 deg');
xlim([1800 3000]);

nexttile
%plot(array_of_r, rad2deg(Error_Kepler_of_r_0(:,4)), '-', 'LineWidth', 2);
xlabel('a_0, km','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('\Delta{\omega}');
grid on;
hold on;
plot(array_of_r, rad2deg(Error_Kepler_of_r_30(:,4)), '-', 'LineWidth', 2);
plot(array_of_r, rad2deg(Error_Kepler_of_r_70(:,4)), '-', 'LineWidth', 2);
plot(array_of_r, rad2deg(Error_Kepler_of_r_90(:,4)), '-', 'LineWidth', 2);
%plot(array_of_r, rad2deg(Error_Kepler_of_r_40(:,4)), '-', 'LineWidth', 2);
hold off;
legend('i=30 deg', 'i=70 deg', 'i=90 deg', 'i=40 deg');
xlim([1800 3000]);

nexttile
%plot(array_of_r, rad2deg(Error_Kepler_of_r_0(:,5)), '-', 'LineWidth', 2);
xlabel('a_0, km','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('\Delta{\Omega}');
grid on;
hold on;
plot(array_of_r, rad2deg(Error_Kepler_of_r_30(:,5)), '-', 'LineWidth', 2);
plot(array_of_r, rad2deg(Error_Kepler_of_r_70(:,5)), '-', 'LineWidth', 2);
plot(array_of_r, rad2deg(Error_Kepler_of_r_90(:,5)), '-', 'LineWidth', 2);
%plot(array_of_r, rad2deg(Error_Kepler_of_r_40(:,5)), '-', 'LineWidth', 2);
hold off;
legend('i=30 deg', 'i=70 deg', 'i=90 deg', 'i=40 deg', 'Location','southeast');
xlim([1800 3000]);

nexttile
plot(array_of_r, rad2deg(Error_Kepler_of_r_0(:,6)), '-', 'LineWidth', 2);
xlabel('a_0, km','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('\Delta{\itM}');
grid on;
hold on;
plot(array_of_r, rad2deg(Error_Kepler_of_r_30(:,6)), '-', 'LineWidth', 2);
plot(array_of_r, rad2deg(Error_Kepler_of_r_70(:,6)), '-', 'LineWidth', 2);
plot(array_of_r, rad2deg(Error_Kepler_of_r_90(:,6)), '-', 'LineWidth', 2);
%plot(array_of_r, rad2deg(Error_Kepler_of_r_40(:,6)), '-', 'LineWidth', 2);
hold off;
legend('i=1 deg', 'i=30 deg', 'i=70 deg', 'i=90 deg', 'i=40 deg', 'Location','southeast');
xlim([1800 3000]);

figure;
plot(array_of_r, rad2deg(Error_Kepler_of_r_30(:,4)+Error_Kepler_of_r_30(:,6)), '-', 'LineWidth', 2);
xlabel('a_0, km','FontWeight','bold');
ylabel('[deg]','FontWeight','bold');
title('\Delta({\omega+\itM})');
grid on;
hold on;
plot(array_of_r, rad2deg(Error_Kepler_of_r_70(:,4)+Error_Kepler_of_r_70(:,6)), '-', 'LineWidth', 2);
plot(array_of_r, rad2deg(Error_Kepler_of_r_90(:,4)+Error_Kepler_of_r_90(:,6)), '-', 'LineWidth', 2);
%plot(array_of_r, rad2deg(Error_Kepler_of_r_40(:,4)+Error_Kepler_of_r_40(:,6)), '-', 'LineWidth', 2);
hold off;
legend('i=30 deg', 'i=70 deg', 'i=90 deg', 'i=40 deg');
xlim([1800 3000]);

%% Creating new arrays to get info from files
array_of_r_0 = array_of_r;
Error_Kepler_of_r_30_without = Error_Kepler_of_r;
Error_r_of_r_30_without = Error_r_of_r;
Error_v_of_r_30_without = Error_v_of_r;
LP_of_r_30_without = LP_of_r;

%%
figure;
tiledlayout(2,1)

nexttile
plot(array_of_r, Error_r_of_r_30_without, '--', 'LineWidth', 2);
title('\Delta{\itr}({\ita}0), km (inc = 30 deg)');
grid on;
hold on;
plot(array_of_r, Error_r_of_r_30_with, '-', 'LineWidth', 2);
hold off;
legend('without Earth gravity effects', 'with Earth gravity effects', 'Location','southeast');
xlim([1800 6000]);

nexttile
plot(array_of_r, Error_v_of_r_30_without, '--', 'LineWidth', 2);
title('\Delta{\itv}({\ita}0), km/s (inc = 30 deg)');
grid on;
hold on;
plot(array_of_r, Error_v_of_r_30_with, '-', 'LineWidth', 2);
hold off;
legend('without Earth gravity effects', 'with Earth gravity effects', 'Location','southeast');
xlim([1800 6000]);