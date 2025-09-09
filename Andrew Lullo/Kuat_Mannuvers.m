% ***********************************
% Helper Functions
% ***********************************
function [v] = vis_viva(mu, a, r)
    v = sqrt(mu * (2 ./ r - 1 ./ a));
end

function [gamma] = flight_path_angle(h, r, v)
    gamma = acosd(h / (r * v));
    if imag(gamma) < 1e-5
        gamma = real(gamma);
    end
end

function [c] = law_of_cosines(a, b, gamma)
    c = sqrt(a ^ 2 + b ^ 2 - 2 * a * b * cosd(gamma));
end

function [e, a, v, p, h] = hyperbola_parameters(mu, r_p, v_inf)
    e = 1 + r_p * v_inf ^ 2 / mu;
    a = -mu / v_inf ^ 2;
    v = vis_viva(mu, a, r_p);
    p = a * (1 - e ^ 2);
    h = sqrt(mu * p);
end

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

% Hohmann to Moon
e_trans2 = .85
a_trans2 = hgeo_kuat / (1 - e_trans2)
thetastar_trans2 = linspace(0, deg2rad(360), 200);
w_trans2 = deg2rad(43);

% Kuat Out
thetastar_trans3 = linspace(deg2rad(-79), deg2rad(79), 200);
w_trans3 = deg2rad(270);

% Given Tatooine Variables
gamma_plus_kuat = 0;
v_inf_plus_kuat = 12.6285;

% ********************
% Calculations
% ********************
a_hyp_kuat = -mu_kuat / v_inf_plus_kuat^2
v_plus_moon = sqrt(mu_kuat * (2 / a_moon - 1 / a_hyp_kuat))

v_orb_kuati = sqrt(mu_kuat / hgeo_kuat);
vp_kuat = sqrt(mu_kuat * (2 / hgeo_kuat - 1 / a_hyp_kuat));
dv_direct = vp_kuat - v_orb_kuati
vp_trans2 = sqrt(mu_kuat * (2 / hgeo_kuat - 1 / a_trans2));
dv_trans2 = vp_trans2 - v_orb_kuati

gamma_plus_moon = 39.554; % Gotten from Travis' analysis
v_moon = sqrt(mu_kuat / a_moon);
v_inf_plus_moon = law_of_cosines(v_moon, v_plus_moon, gamma_plus_moon)

v_minus_moon = sqrt(mu_kuat * (2 / a_moon - 1 / a_trans2))

thetastar_intersect = acosd((a_trans2 * (1 - e_trans2^2) / a_moon - 1) / e_trans2);
gamma_minus_moon = atand(e_trans2 * sind(thetastar_intersect) / (1 + e_trans2 * cosd(thetastar_intersect)))
v_inf_minus_moon = sqrt(v_minus_moon^2 + v_moon^2 - 2 * v_moon * v_minus_moon * cosd(gamma_minus_moon))

dV_needed = sqrt(v_minus_moon^2 + v_plus_moon^2 - 2 * v_plus_moon * v_minus_moon * cosd(gamma_minus_moon - gamma_plus_moon))
delta_needed = acosd((v_inf_plus_moon^2 + v_inf_minus_moon^2 - dV_needed^2) / (2 * v_inf_plus_moon * v_inf_minus_moon))

a_hyp_moon_needed = -mu_moon / v_inf_plus_moon^2
a_hyp_moon_actual = -mu_moon / v_inf_minus_moon^2

[r_pfb, obval] = fsolve(@(r_p) asind(1 / (1 + r_p * v_inf_minus_moon^2 / mu_moon)) + asind(1 / (1 + r_p * v_inf_plus_moon^2 / mu_moon)) - delta_needed, 1);
r_pfb

delta_minus = 2 * asind(1 / (1 + r_pfb * v_inf_minus_moon^2 / mu_moon))
delta_plus = 2 * asind(1 / (1 + r_pfb * v_inf_plus_moon^2 / mu_moon))

a_plus = -mu_moon / v_inf_plus_moon^2;
a_minus = -mu_moon / v_inf_minus_moon^2;

vp_plus = vis_viva(mu_moon, a_plus, r_pfb)
vp_minus = vis_viva(mu_moon, a_minus, r_pfb)

dv_kick = vp_plus - vp_minus

% ***********************************
% Final Plotting Settup
% ***********************************
h_exit = a_moon * v_plus_moon * sin(gamma_plus_moon);
p_exit = h_exit^2 / mu_kuat;
e_exit = sqrt(-(p_exit / a_hyp_kuat - 1));


% **********
% Plots
% **********

% Plot Bador
orbitplot2D(a_moon, e_moon, thetastar_moon, w_moon, "Bador Orbit", r_scale = AU);
% Plot Kuat Sphere of Influence
orbitplot2D(r_soi, e_ksoi, thetastar_ksoi, w_ksoi, "Kuat Sphere of Influence", r_scale = AU, LineStyle='-.');
% Plot Kuat
orbitplot2D(r_kaut, e_kuat_planet, thetastar_kuat_planet, w_kuat_planet, "Kuat Planet", r_scale = AU, color='r'); hold on
% Plot Kuatistationary
orbitplot2D(hgeo_kuat, e_kuat_geo, thetastar_kuat_geo, w_kuat_geo, "Kuatistationary Orbit", r_scale = AU);
% Plot Bador Sphere of Influence
orbitplot2D(r_msoi, e_msoi, thetastar_msoi, w_msoi, "Bador Sphere of Influence", r_scale = AU, origin = [-100000; -a_moon + 10000], LineStyle='-.');
% Plot transfer to Bador
orbitplot2D(a_trans2, e_trans2, thetastar_trans2, w_trans2, "Transfer to Bador", r_scale = AU, LineStyle = "--");
% Plot transfer Out
orbitplot2D(a_hyp_kuat, e_exit, thetastar_trans3, w_trans3, "Exit Orbit", r_scale = AU, LineStyle = "--");

% **********
% Trajectory
% **********

thetastar_trans2 = linspace(0, deg2rad(-thetastar_intersect), 200);
thetastar_trans3 = linspace(deg2rad(-16.8), deg2rad(-79), 100);

% Plot transfer to Bador
orbitplot2D(a_trans2, e_trans2, thetastar_trans2, w_trans2, "Spacecraft Initial Trajectory", r_scale = AU, LineStyle = "-", color='y');
% Plot transfer Out
orbitplot2D(a_hyp_kuat, e_exit, thetastar_trans3, w_trans3, "Spacecraft Final Trajectory", r_scale = AU, LineStyle = "-", color='y');

hold off
grid on
axis equal
title("Kuat Mannuvers")
xlabel("X [AU]")
ylabel("Y [AU]")
legend()