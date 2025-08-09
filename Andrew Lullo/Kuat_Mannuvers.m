% ***********************************
% Constants/Plotting Settup
% ***********************************

% The moon Bador is the same as our moon
r_kaut = 10000;
a_moon = 384400;
mu_moon = 4902.8005821478;
r_soi = a_kuat * (mu_kuat / mu_Sun)^(2/5);
r_msoi = a_moon * (mu_moon / mu_kuat)^(2/5);

% Kuat Plot
e_kuat_planet = 0;
thetastar_kuat_planet = linspace(0, deg2rad(360), 200);
w_kuat_planet = deg2rad(90);

% Kuatistationary Plot
e_kuat_geo = 0;
thetastar_kuat_geo = linspace(0, deg2rad(360), 200);
w_kuat_geo = deg2rad(90);

% Bador Plot
e_moon = 0;
thetastar_moon = linspace(0, deg2rad(360), 200);
w_moon = deg2rad(90);

% Kuat Sphere of Influence Plot
e_ksoi = 0;
thetastar_ksoi = linspace(0, deg2rad(360), 200);
w_ksoi = deg2rad(90);

% Moon Sphere of Influence Plot
e_msoi = 0;
thetastar_msoi = linspace(0, deg2rad(360), 200);
w_msoi = deg2rad(90);

% Hohmann in
new_thetastar_trans1 = linspace(deg2rad(179), deg2rad(180), 10);

% Hohmann to Moon
thetastar_trans2 = linspace(0, deg2rad(180), 100);
w_trans2 = deg2rad(90);

% ********************
% Calculations
% ********************

rp_hyp = r_bador + 500;
a_trans2 = (hgeo_kuat + a_moon + rp_hyp) / 2;
e_trans2 = -((hgeo_kuat / a_trans2) - 1)



% **********
% Plots
% **********

% Plot Kuat
orbitplot2D(r_kaut, e_kuat_planet, thetastar_kuat_planet, w_kuat_planet, "Kuat", r_scale = AU); hold on
% Plot Kuatistationary
orbitplot2D(hgeo_kuat, e_kuat_geo, thetastar_kuat_geo, w_kuat_geo, "Kuatistationary", r_scale = AU);
% Plot Bador
orbitplot2D(a_moon, e_moon, thetastar_moon, w_moon, "Moon", r_scale = AU);
% Plot Kuat Sphere of Influence
orbitplot2D(r_soi, e_ksoi, thetastar_ksoi, w_ksoi, "Kuat Sphere of Influence", r_scale = AU);
% Plot Bador Sphere of Influence
orbitplot2D(r_msoi, e_msoi, thetastar_msoi, w_msoi, "Bador Sphere of Influence", r_scale = AU, origin = [0; -a_moon]);
% Plot transfer in
orbitplot2D(a_trans1, e_trans1, new_thetastar_trans1, w_trans1, "Transfer 1", r_scale = AU, LineStyle = "-.", origin = [0; a_kuat]);
% Plot transfer in
orbitplot2D(a_trans2, e_trans2, thetastar_trans2, w_trans2, "Transfer 2", r_scale = AU, LineStyle = "-.");
hold off
grid on
axis equal
title("Plot of Transfer between Kuat and Endor Orbits")
xlabel("X [AU]")
ylabel("Y [AU]")
legend()