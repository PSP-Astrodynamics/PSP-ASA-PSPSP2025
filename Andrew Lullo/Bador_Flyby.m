% ***********************************
% Constants/Plotting Settup
% ***********************************

% Bador Plot
e_bador = 0;
r_bador = 1738.2;
thetastar_bador = linspace(0, deg2rad(360), 200);
w_bador = deg2rad(90);



% ********************
% Calculations
% ********************

ra_trans2 = a_trans2 * (1 + e_trans2);
va_trans2 = sqrt(mu_moon * ((2 / ra_trans2) - (1 / a_trans2)));
((2 / ra_trans2) - (1 / a_trans2))

va_trans2^2 - (2 * mu_moon / rp_hyp)
vinf_hyp = sqrt(va_trans2^2 - (2 * mu_moon / rp_hyp));

a_hyp = -1 * mu_moon / vinf_hyp^2

% **********
% Plots
% **********

% Plot Bador
orbitplot2D(r_bador, e_bador, thetastar_bador, w_bador, "Bador", r_scale = AU); hold on
% Plot Bador Sphere of Influence
orbitplot2D(r_msoi, e_msoi, thetastar_msoi, w_msoi, "Bador Sphere of Influence", r_scale = AU);
hold off
grid on
axis equal
title("Plot of Transfer between Kuat and Endor Orbits")
xlabel("X [AU]")
ylabel("Y [AU]")
legend()