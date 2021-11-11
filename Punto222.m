%% ----------------- Problem 2.22 --------------------------
%This script is just an implementation of the equation shown in the inform

n = 5;
nb = 3;

suma = 0;
for k=1:1:2
	suma = suma + exp( (-5.55*k^2) / (nb - 1)^2 );
end

BSL = n / ( 1 + 2*suma ); % Beam Shape Loss (dB)
twoWayBSL = 2*BSL; % Two Way Beam Shape Loss (dB)
 

f = 9375e6;
l = c/f;