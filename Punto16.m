%% ----------------- Problem 1.6 --------------------------
clc;

syms A(Pt)
syms Pt positive

Smin = 3e-13; %Receiver minimum detectable signal
c = physconst('LightSpeed'); %Speed of light
lambda = c/(1230e6); %Wavelength
R = 1852 * 200; %Range in meters
sigma = 2; %Target cross section
rho_e = 0.6; %Antenna aperture efficiency

%Definition of the area as a function of the power
A(Pt) = ((Smin*4*pi*lambda^2*R^4)/(Pt*sigma*rho_e^2))^(1/2);

%Definition of the total cost
Ctotal = 1400*A + 2.20*Pt + 1e6;

%Derivate of the total cost
dCtotal = diff(Ctotal, Pt);

Pt_opt = double(solve(dCtotal == 0, Pt)); %Optimal power
Ctot_opt = (vpa(Ctotal(Pt_opt), 8)); %Optimal total cost

%Plotting the variables of interest
figure();
fplot(Ctotal, [1e-9, 1e12],'LineWidth', 1.5);
set(gca, 'XScale','log')
grid on, grid minor;
%ylabel("Total Cost, C_tot [$]");
ylabel('Total Cost, $C_{tot}$ [\$]','Interpreter','latex','FontSize', 15);
xlabel('Transmitter peak power, $P_t$ [W]' ,'Interpreter','latex', 'FontSize', 15);
title('Total Cost vs Transmitter peak power' ,'Interpreter','latex', 'FontSize', 15);
%ylim([0, 1e9])
%hold on;
%xline(Pt_opt, '-r');


figure();
subplot(2, 1, 1);
fplot(Ctotal, [1e3, 0.5e6],'LineWidth', 1.5);
%set(gca, 'XScale','log')
ylim([0, 3e6])
grid on, grid minor;
hold on;
ax = gca;
ax.XAxis.Exponent = 3;
ax.YAxis.Exponent = 6;
xline(Pt_opt, '--r', 'P_t = 84 kW', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
yline(double(Ctotal(Pt_opt)), '-.r', '$1,554,411', 'LineWidth', 1.2, 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'center');
ylabel('Total Cost, $C_{tot}$ [\$]','Interpreter','latex','FontSize', 15);
xlabel('Transmitter peak power, $P_t$ [W]' ,'Interpreter','latex', 'FontSize', 15);
title('Total Cost vs Transmitter peak power' ,'Interpreter','latex', 'FontSize', 15);

subplot(2, 1, 2);
fplot(dCtotal, [1e3, 0.5e6],'LineWidth', 1.5);
%set(gca, 'XScale','log')
ylim([-10, 5])
grid on, grid minor;
hold on;
xline(Pt_opt, '--r', 'P_t = 84 kW', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
yline(0, '-')
ax = gca;
ax.XAxis.Exponent = 3;
%ax.YAxis.Exponent = 6;
title('Total Cost derivative vs Transmitter peak power' ,'Interpreter','latex', 'FontSize', 15);
ylabel('Total Cost derivativev, $\frac{d\, C_{tot}}{d\, P_t}$ [$\frac{\$}{W}$]','Interpreter','latex','FontSize', 15);
xlabel('Transmitter peak power, $P_t$ [W]' ,'Interpreter','latex', 'FontSize', 15);




