%% ----------------- Problem 2.24 --------------------------
clc;

%Definition of variables used in the problem
f = 9400e6; %Operation frequency
c = physconst('LightSpeed'); %Speed of light
lambda = c/f; %Wavelength
thetaH = 0.8; %deg
G = 33; %dB Gain
rotation = 20; %rpm
Pt = 25e3; %W Power
tp = 0.15e-6; %pulse-width
fp = 4000; %Hz pulse repetition rate
Fn = 5; %dB noise figure
B = 15e6; % Reciver bandwidth
Ls = 12; %dB system losses
Tfa = 14*60*60; %seg time between false alarm

sigma = 10; %m^2 Target cross section

%Calculation of the number of pulses in each scan
thetaS = rotation*360/60; %deg/s

%Calculation of the value of n
n = thetaH*fp/thetaS; %% 26 approx

% Using the reference graphic, we choose an average value
nEi = 14;

omega = 2*pi*sin(pi/12); %azimuth * sin(THelevation)
k = 1.38e-23; %Boltzman constant 
To = 290; %K Temperature
kTo = 4e-21;

G = 10^(G/10);
Fn = 10^(Fn/10);
Ls = 10^(Ls/10);

Pfa = 1/(Tfa*B); %Calculation of the Probability of false alarm
Pd = 0.3:0.001:0.99; %Sweep of the Probability of detection

%Calculation of the necessary constants for the signal-to-noise ratio
A = log(0.62/Pfa);
B = log(Pd./(1-Pd));

%Calculation of the signal-to-noise ratio
SN = A + 0.12*A.*B + B;
SN = 10.^(SN/10);

%Calculation of the maximum range as a function of the Probability of
%detection
Rmax4 = ( Pt*tp*G*lambda^2*sigma*nEi ) / ( (4*pi)^2*kTo*Fn*Ls*fp*omega ) ./ SN;
Rmax = (Rmax4).^(1/4)*0.005399568; % nmi

%Plotting the graphic required in this subsection
figure();
plot(Rmax, Pd,'LineWidth', 1.5)
xlim([0.5, 2.2])
grid on, grid minor;
hold on;
title('Probability of Detection as a function of Range' ,'Interpreter','latex', 'FontSize', 15);
ylabel('Probability of Detection','Interpreter','latex','FontSize', 15);
xlabel('Maximum Range, $R_{max}$ [nmi]' ,'Interpreter','latex', 'FontSize', 15);

% Swerling fluctuation model (Phased array system toolbox required)
SNR_Sw1 = zeros(1, length(Pd));
%SNR_Sw3 = zeros(1, length(Pd));
SNR_Sw0 = zeros(1, length(Pd));

N = 1; 

for ii = 1:length(Pd)
    SNR_Sw1(ii) = shnidman(Pd(ii), Pfa, N, 1);
    %SNR_Sw3(ii) = shnidman(Pd(ii), Pfa, N, 3);
    SNR_Sw0(ii) = shnidman(Pd(ii), Pfa, N, 0);
end

Lf = -SNR_Sw0 + SNR_Sw1;

%Plotting the graphic of the Fluctuaction loss as a function of the
%Probability of detection
figure();
hold on;
plot(Pd, Lf, 'LineWidth', 1.5);
grid on, grid minor;
ylabel('$L_f$, dB','Interpreter','latex','FontSize', 15);
xlabel('Probability of detection, $P_d$' ,'Interpreter','latex', 'FontSize', 15);
title('Fluctuaction loss (dB) vs Probability of detection' ,'Interpreter','latex', 'FontSize', 15);

Lf = 10.^(Lf/10);

Rmax4_Lf = Rmax4./Lf;
Rmax_Lf = (Rmax4_Lf).^(1/4)*0.005399568; % nmi

%Plotting the Probability of detection for both cases
figure();
plot(Rmax, Pd,'LineWidth', 1.5)
hold on;
plot(Rmax_Lf, Pd,'LineWidth', 1.5)
grid on, grid minor;
legend("No fluctuaction", "Swerling I")
title('Probability of Detection as a function of Range' ,'Interpreter','latex', 'FontSize', 15);
ylabel('Probability of Detection','Interpreter','latex','FontSize', 15);
xlabel('Maximum Range, $R_{max}$ [nmi]' ,'Interpreter','latex', 'FontSize', 15);