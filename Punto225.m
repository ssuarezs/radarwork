%% ----------------- Problem 2.25 --------------------------
clc;

%Definition of variables used in the problem
f = 2.8e9; %Operation frequency
c = physconst('LightSpeed'); %Speed of light
lambda = c/f; %Wavelength
thetaH = 1.35; %deg
G = 33; %dB Gain
rotation = 12.8; %rpm
Pt = 1.4e6; %W Power
tp = 0.6e-6; %pulse-width
fp = 1040; %Hz Pulse repetition rate
Fn = 4; %dB Receiver noise figure
B = 15e6; % Receiver bandwidth
Ls = 12; %dB System losses
Tfa = 20*60; %sec Time between false alarm

%Calculation of the Probability of false alarm with the Pulse repetition
%rate and the time between false alarm
Pfa = 1/(20*60*fp);

sigma = 2; %m^2 Target cross section

%Calculation of the number of pulses in each scan
thetaS = rotation*360/60; %deg/s

%Calculation of the value of n
n = thetaH*fp/thetaS; %% 18 approx

% Using the reference graphic, we choose an average value
nEi = 10.5;

omega = 1; %azimuth * sin(THelevation)
k = 1.38e-23; %Boltzman constant
To = 290; %K Temperature
kTo = 4e-21; 

%Take the variables into the natural scale
G = 10^(G/10);
Fn = 10^(Fn/10);
Ls = 10^(Ls/10);

%Calculation of the average power
Pav = Pt*tp*fp;

%Solution of the subsection a.)

Pd = (1/3):0.001:0.99; %Sweep of the Probability of detection

%Calculation of the necessary constants for the signal-to-noise ratio
A = log(0.62/Pfa);
B = log(Pd./(1-Pd));

%Calculation of the signal-to-noise ratio
SN = A + 0.12*A.*B + B;
SN = 10.^(SN/10);

%Calculation of the maximum range as a function of the Probability of
%detection
Rmax4 = ( Pt*tp*G*lambda^2*sigma*nEi ) / ( (4*pi)^2*kTo*Fn*Ls*fp*omega ) ./ SN;
Rmax_a = (Rmax4).^(1/4)*0.005399568; % nmi

%Plotting the graphic required in this subsection
figure();
plot(Rmax_a, Pd,'LineWidth', 1.5)
%xlim([5.5, 14])
grid on, grid minor;
hold on;
title('Probability of Detection as a function of Range' ,'Interpreter','latex', 'FontSize', 15);
ylabel('Probability of Detection','Interpreter','latex','FontSize', 15);
xlabel('Maximum Range, $R_{max}$ [nmi]' ,'Interpreter','latex', 'FontSize', 15);
xlim([5, 35]);
ylim([0.3, 1]);

%Solution of the subsection b.)

Pd = (2/3):0.001:0.99; %Sweep of the Probability of detection

%Calculation of the necessary constants for the signal-to-noise ratio
A = log(0.62/Pfa);
B = log(Pd./(1-Pd));

%Calculation of the signal-to-noise ratio
SN = A + 0.12*A.*B + B;
SN = 10.^(SN/10);

%Calculation of the maximum range as a function of the Probability of
%detection
Rmax4 = ( Pt*tp*G*lambda^2*sigma*nEi ) / ( (4*pi)^2*kTo*Fn*Ls*fp*omega ) ./ SN;
Rmax_b = (Rmax4).^(1/4)*0.005399568; % nmi

%Plotting the graphic required in this subsection
figure();
plot(Rmax_b, Pd,'LineWidth', 1.5)
%xlim([5.5, 14])
grid on, grid minor;
hold on;
title('Probability of Detection as a function of Range' ,'Interpreter','latex', 'FontSize', 15);
ylabel('Probability of Detection','Interpreter','latex','FontSize', 15);
xlabel('Maximum Range, $R_{max}$ [nmi]' ,'Interpreter','latex', 'FontSize', 15);
xlim([5, 35]);
ylim([0.3, 1]);

%Next, we will made the graphic of the Fluctuaction losses, which are
%necessary to use the Shnidman formulae
N = 1; 

Pd = 0.1:0.001:0.99;
Lf = zeros(1, length(Pd));

for ii = 1:length(Pd)
    Lf(ii) = shnidman(Pd(ii), Pfa, N, 1) - shnidman(Pd(ii), Pfa, N, 0);
end

figure();
hold on;
plot(Pd,Lf, 'LineWidth', 1.5);
grid on, grid minor;
ylabel('$L_f$, dB','Interpreter','latex','FontSize', 15);
xlabel('Probability of detection, $P_d$' ,'Interpreter','latex', 'FontSize', 15);
title('Fluctuaction loss (dB) vs Probability of detection' ,'Interpreter','latex', 'FontSize', 15);

%It is important to say that in order to solve the problems in the 
%subsections c and d, we used the toolbox Phases array systems

%Solution of the subsection c.)
Pd = (1/3):0.001:0.99; %Sweep of the Probability of detection

%Calculation of the necessary constants for the signal-to-noise ratio
A = log(0.62/Pfa);
B = log(Pd./(1-Pd));

%Calculation of the signal-to-noise ratio
SN = A + 0.12*A.*B + B;
SN = 10.^(SN/10);

%Calculation of the maximum range as a function of the Probability of
%detection which is the same calculated in the subsection a.
Rmax4 = ( Pt*tp*G*lambda^2*sigma*nEi ) / ( (4*pi)^2*kTo*Fn*Ls*fp*omega ) ./ SN;

%Calculation of the maximum range for a fluctuating target, using the
%Shnidman empirical formulae
Lf = zeros(1, length(Pd));
for ii = 1:length(Pd)
    Lf(ii) = shnidman(Pd(ii), Pfa, N, 1) - shnidman(Pd(ii), Pfa, N, 0);
end
Lf = 10.^(Lf/10);

Rmax4_Lf = Rmax4./Lf;
Rmax_Lf = (Rmax4_Lf).^(1/4)*0.005399568; % nmi

%Plotting the graphic required in this subsection
figure();
plot(Rmax_a, Pd,'LineWidth', 1.5)
hold on;
plot(Rmax_Lf, Pd,'LineWidth', 1.5)
grid on, grid minor;
legend("No fluctuaction", "Swerling I")
title('Probability of Detection as a function of Range' ,'Interpreter','latex', 'FontSize', 15);
ylabel('Probability of Detection','Interpreter','latex','FontSize', 15);
xlabel('Maximum Range, $R_{max}$ [nmi]' ,'Interpreter','latex', 'FontSize', 15);
xlim([5, 35]);
ylim([0.3, 1]);

%Solution of the subsection d.)
Pd = (2/3):0.001:0.99; %Sweep of the Probability of detection

%Calculation of the necessary constants for the signal-to-noise ratio
A = log(0.62/Pfa);
B = log(Pd./(1-Pd));

%Calculation of the signal-to-noise ratio
SN = A + 0.12*A.*B + B;
SN = 10.^(SN/10);

%Calculation of the maximum range as a function of the Probability of
%detection which is the same calculated in the subsection b.
Rmax4 = ( Pt*tp*G*lambda^2*sigma*nEi ) / ( (4*pi)^2*kTo*Fn*Ls*fp*omega ) ./ SN;

%Calculation of the maximum range for a fluctuating target, using the
%Shnidman empirical formulae
Lf = zeros(1, length(Pd));
for ii = 1:length(Pd)
    Lf(ii) = shnidman(Pd(ii), Pfa, N, 1) - shnidman(Pd(ii), Pfa, N, 0);
end
Lf = 10.^(Lf/10);

Rmax4_Lf = Rmax4./Lf;
Rmax_Lf = (Rmax4_Lf).^(1/4)*0.005399568; % nmi

%Plotting the graphic required in this subsection
figure();
plot(Rmax_b, Pd,'LineWidth', 1.5)
hold on;
plot(Rmax_Lf, Pd,'LineWidth', 1.5)
grid on, grid minor;
legend("No fluctuaction", "Swerling I")
title('Probability of Detection as a function of Range' ,'Interpreter','latex', 'FontSize', 15);
ylabel('Probability of Detection','Interpreter','latex','FontSize', 15);
xlabel('Maximum Range, $R_{max}$ [nmi]' ,'Interpreter','latex', 'FontSize', 15);
xlim([5, 35]);
ylim([0.3, 1]);

%Solution of the subsection e.)

%Calculation of the maximum unambigous range with the pulse repetition rate
R_un = c/(2*fp*1852); % nmi