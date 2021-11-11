%% ----------------- Problem 1.5 --------------------------

R = 152.4; % 500ft
Smin = 5e-13; %Minimum detectable signal
Ae = 0.0557; %Area of the antenna
G = 648.005; %Gain in natural units
sigma = 10; %Radar cross section

Pt = Smin * ((4*pi)^2*R^4)/(G*Ae*sigma); %Calculation of the average power