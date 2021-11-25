
%% Exercises on Antennas

%% 1. Exc 9.17

clc;

f = 1;
c = physconst('lightspeed');
lambda = c/f;

k = 2*pi/lambda;

d = lambda/2; % Array element distance


l = 1*lambda;
e_p = [-0.2; 0; 0.2];

beta = (2*pi*l/lambda).*(e_p./(1+e_p));

theta = 0:0.01:pi; % Scan angle

Psi(1,:) = k*d*cos(theta) + beta(1);
Psi(2,:) = k*d*cos(theta) + beta(2);
Psi(3,:) = k*d*cos(theta) + beta(3);

N = 5;

AFn = sin(N*Psi/2)./(N*sin(Psi/2)); % Normalized array factor

[m, i] = max(AFn.');

t_min = theta(i(1));
t_max = theta(i(3));

figure();
plot(rad2deg(theta), AFn, 'LineWidth', 1.3);
grid on, grid minor; 
xlim(rad2deg([0, pi]));
title("Normalized array factor vs scan angle",'Interpreter','latex', 'FontSize', 15);
xlabel("Scan angle $\theta$ [deg]",'Interpreter','latex', 'FontSize', 15);
ylabel("Normalized array factor",'Interpreter','latex', 'FontSize', 15);
xline(90, '--r', "$90^\circ$", 'LineWidth', 0.8, 'Interpreter','latex', 'FontSize', 12, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'horizontal')
xline(rad2deg(t_min), '--b', "$60^\circ$", 'LineWidth', 0.8, 'Interpreter','latex', 'FontSize', 12, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'horizontal')
xline(rad2deg(t_max), '--', "$109.4^\circ$", 'Color', '#EDB120', 'LineWidth', 1, 'Interpreter','latex', 'FontSize', 12, 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'horizontal')
legend("$f = f_0-20\%$","$f = f_0$", "$f = f_0+20\%$" ,'Interpreter','latex', 'FontSize', 10);



f_e = f+f*ef;

theta_0_r = asin((l/d)*(1-(f./f_e)));
theta_0_d = rad2deg(theta_0_r)




