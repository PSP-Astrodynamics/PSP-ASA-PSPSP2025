% Constants
mu_Sun = 132712440017.99; % [km3 / s2] Sun gravitational parameter

AU = 149597870.69;

a_vulcan = 2.1 * AU;
mu_vulcan = 2.8307e7;

a_kuat = 0.6 * AU;
mu_kuat = 1.2748e6;

kuat_day = 60 * 60 * 20;
hgeo_kuat = ((mu_kuat * (kuat_day)^2) / (4 * pi^2))^(1/3);

% Add helper function folder to path (for plotting)
addpath("Helper Functions") 

%% 

v_orb_Vlc = sqrt(mu_Sun / a_vulcan);
v_orb_Kut = sqrt(mu_Sun / a_kuat);

a_trans1 = (a_vulcan + a_kuat + hgeo_kuat) / 2;
e_trans1 = -((a_vulcan / a_trans1) - 1)
P_trans1 = 2 * pi * sqrt(a_trans1^3 / mu_Sun)

v_trans1 = sqrt((2 / a_vulcan) - (1 / a_trans1));
v_trans2 = sqrt((2 / a_kuat) - (1 / a_trans1));

dV1 = v_trans1 - v_orb_Vlc
dV2 = v_orb_Kut - v_trans2

%%

% Plot
thetastar_trans1 = linspace(0, deg2rad(180), 100);
w_trans1 = deg2rad(90);

e_kuat = 0;
% Create a true anomaly vector from 0 to 180 degrees with 100 elements
thetastar_kuat = linspace(0, deg2rad(360), 200);
w_kuat = deg2rad(40);

e_vulcan = 0;
% Create a true anomaly vector from 0 to 270 degrees with 300 elements
thetastar_vulcan = linspace(0, deg2rad(360), 300);
w_vulcan = deg2rad(90);


% Plot Kuat
orbitplot2D(a_kuat, e_kuat, thetastar_kuat, w_kuat, "Kuat", r_scale = AU); hold on
% Plot Vulcan
orbitplot2D(a_vulcan, e_vulcan, thetastar_vulcan, w_vulcan, "Vulcan", r_scale = AU, LineStyle = "-")
% Plot transfer
orbitplot2D(a_trans1, e_trans1, thetastar_trans1, w_trans1, "Transfer", r_scale = AU, LineStyle = "-.")
hold off
grid on
axis equal
title("Plot of Transfer between Kuat and Endor Orbits")
xlabel("X [AU]")
ylabel("Y [AU]")
legend()