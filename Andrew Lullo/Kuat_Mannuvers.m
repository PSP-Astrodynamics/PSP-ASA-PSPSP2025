% ***********************************
% Constants/Plotting Settup
% ***********************************

% The moon Bador is the same as our moon
r_kaut = 10000;
a_moon = 384400;
mu_moon = 4902.8005821478 * 10;
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

% Hohmann In
new_thetastar_trans1 = linspace(deg2rad(179), deg2rad(180), 10);

% Hohmann to Moon
a_trans2 = (1.6 * a_moon) / 2;
e_trans2 = -((hgeo_kuat / a_trans2) - 1);
thetastar_trans2 = linspace(0, deg2rad(360), 200);
w_trans2 = deg2rad(70);

% Kuat Out
e_trans3 = 2;
a_trans3 =  -a_moon / (e_trans3 - 1);
thetastar_trans3 = linspace(deg2rad(-80), deg2rad(80), 100);
w_trans3 = deg2rad(270);

% Given Tatooine Variables
gamma_plus_kuat = 1.7075e-06;
v_inf_plus_kuat = 12.6285;

% ********************
% Calculations
% ********************
a_hyp_kuat = -mu_kuat / v_inf_plus_kuat^2
v_plus_moon = sqrt(mu_kuat * (2 / a_moon - 1 / a_hyp_kuat))

gamma_plus_moon = 0; % The moon is at Rp for the exit hyperbola right?
v_moon = sqrt(mu_kuat / a_moon);
v_inf_plus_moon = v_plus_moon - v_moon

v_minus_moon = sqrt(mu_kuat * (2 / a_moon - 1 / a_trans2))

thetastar_intersect = acosd((a_trans2 * (1 - e_trans2^2) / a_moon - 1) / e_trans2)
gamma_minus_moon = atand(e_trans2 * sind(thetastar_intersect) / (1 + e_trans2 * cosd(thetastar_intersect)))
v_inf_minus_moon = sqrt(v_minus_moon^2 + v_moon^2 - 2 * v_moon * v_minus_moon * cosd(gamma_minus_moon))


delta_expected

dV_needed = sqrt(v_minus_moon^2 + v_plus_moon^2 - 2 * v_plus_moon * v_minus_moon * cosd(gamma_minus_moon - gamma_plus_moon))
delta_needed = acosd((v_inf_plus_moon^2 + v_inf_minus_moon^2 - dV_needed^2) / (2 * v_inf_plus_moon * v_inf_minus_moon))

a_hyp_moon_needed = -mu_moon / v_inf_plus_moon^2
a_hyp_moon_actual = -mu_moon / v_inf_minus_moon^2

[r_pfb, obval] = fsolve(@(r_p) asind(1 / (1 + r_p * v_inf_minus_moon^2 / mu_moon)) + asind(1 / (1 + r_p * v_inf_plus_moon^2 / mu_moon)) - delta_needed, 1);
r_pfb

vp_actual = sqrt(mu_kuat * (2 / r_pfb - 1 / a_hyp_moon_actual))
vp_needed = sqrt(mu_kuat * (2 / r_pfb - 1 / a_hyp_moon_needed))

impulse = vp_needed - vp_actual


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
% Plot transfer In
%orbitplot2D(a_trans1, e_trans1, new_thetastar_trans1, w_trans1, "Transfer 1", r_scale = AU, LineStyle = "-.", origin = [0; a_kuat]);
% Plot transfer to Bador
orbitplot2D(a_trans2, e_trans2, thetastar_trans2, w_trans2, "Transfer 2", r_scale = AU, LineStyle = "-.");
% Plot transfer Out
orbitplot2D(a_trans3, e_trans3, thetastar_trans3, w_trans3, "Transfer 3", r_scale = AU, LineStyle = "-.");
hold off
grid on
axis equal
title("Kuat Mannuvers")
xlabel("X [AU]")
ylabel("Y [AU]")
legend()