% assume both start at periapsis

AU = 149597898; % [km]
mu_sun = 132712440017.99; % [km3 / s2]

initial_orbit.a = a_kuat * AU; % [km]
initial_orbit.e = 0;
initial_orbit.thetastar = 0; % [rad]
initial_orbit.w = 0; % [rad]
initial_orbit.mu = mu_kuat; % [km3 / s2]

target_orbit.a = a_endor * AU; % [km]
target_orbit.e = 0;
target_orbit.thetastar = 0; % [rad]
target_orbit.w = 0; % [rad]
target_orbit.mu = mu_endor; % [km3 / s2]

N_depart = 1e2;
N_arrive = 1e2;

%%
P_0 = period(initial_orbit, mu_sun);
P_f = period(target_orbit, mu_sun);

P_syn = 1 / abs(1 / P_0 - 1 / P_f); % Synodic period 

t_depart = linspace(0, P_syn, N_depart);
t_arrive = linspace(0, P_syn, N_arrive);

thetastar_depart = time_to_thetastar(t_depart, initial_orbit, mu_sun);
thetastar_arrive = time_to_thetastar(t_arrive, target_orbit, mu_sun);

xy_depart = kepler2D_to_cartestian(initial_orbit, thetastar_depart);
xy_arrive = kepler2D_to_cartestian(target_orbit, thetastar_arrive);

%%
r_scale = AU;

plot(xy_depart(1, :) / r_scale, xy_depart(2, :) / r_scale, LineStyle="-", DisplayName = "Kuat"); hold on
plot(xy_arrive(1, :) / r_scale, xy_arrive(2, :) / r_scale, LineStyle="-", DisplayName = "Endor");
hold off
title("Test Orbits")
xlabel("X")
ylabel("Y")
grid on
axis equal

%%
c = sqrt(xy_depart)
s = (r1 + r2 + c) / 2;

for i = N_depart
    for j = N_arrive
        lambertSolverSMA(t_arrive(t_depart(i) - t_arrive(j), c, s, mu)
    end
end

%% Helper Functions
function [r] = orbit_equation(orbit_struct)
    r = orbit_struct.a * (1 - orbit_struct.e ^ 2) / (1 + orbit_struct.e * cos(orbit_struct.thetastar - orbit_struct.w));
end

function [P] = period(orbit_struct, mu)
    P = 2 * pi * sqrt(orbit_struct.a ^ 3 / mu);
end

function [xy] = kepler2D_to_cartestian(orbit_struct, thetastar)
    r = orbit_equation(orbit_struct);
    x = r .* cos(thetastar + orbit_struct.w);
    y = r .* sin(thetastar + orbit_struct.w);

    xy = [x; y];
end

function [thetastar] = time_to_thetastar(time_period, orbit_struct, mu)
    mean_anomaly = sqrt(mu / orbit_struct.a ^ 3) * time_period;
    
    thetastar = zeros([1, numel(time_period)]);
    for i = 1 : numel(time_period)
        thetastar(i) = eccentric_to_true_anomaly(mean_to_eccentric_anomaly(mean_anomaly(i), orbit_struct.e), orbit_struct.e);
    end
end