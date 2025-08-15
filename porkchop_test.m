% assume both start at periapsis

AU = 149597898; % [km]
mu_sun = 132712440017.99; % [km3 / s2]

initial_orbit.a = a_kuat * AU; % [km]
initial_orbit.e = 0;
initial_orbit.thetastar = 0; % [rad]
initial_orbit.w = 0; % [rad]
initial_orbit.mu = mu_kuat; % [km3 / s2]

target_orbit.a = 1 * AU; % [km]
target_orbit.e = 0;
target_orbit.thetastar = 0; % [rad]
target_orbit.w = 0; % [rad]
target_orbit.mu = mu_endor; % [km3 / s2]

N_depart = 1e2;
N_arrive = 1e2;

% Compute semilatus rectum
p_1 = initial_orbit.a * (1 - initial_orbit.e ^ 2);
p_2 = target_orbit.a * (1 - target_orbit.e ^ 2);

%%
P_0 = period(initial_orbit, mu_sun);
P_f = period(target_orbit, mu_sun);

P_syn = 1 / abs(1 / P_0 - 1 / P_f); % Synodic period 

t_depart = linspace(0, P_f / 2, N_depart);
t_arrive = linspace(P_f / 2, P_f, N_arrive);

thetastar_depart = time_to_thetastar(t_depart, initial_orbit, mu_sun);
thetastar_arrive = time_to_thetastar(t_arrive, target_orbit, mu_sun);

xy_depart = kepler2D_to_cartestian(initial_orbit, thetastar_depart);
xy_arrive = kepler2D_to_cartestian(target_orbit, thetastar_arrive);

%%
r_scale = AU;

plot(xy_depart(1, :) / r_scale, xy_depart(2, :) / r_scale, LineStyle="-", DisplayName = "Kuat"); hold on
plot(xy_arrive(1, :) / r_scale, xy_arrive(2, :) / r_scale, LineStyle="-", DisplayName = "Endor");
scatter(xy_depart(1, 1) / r_scale, xy_depart(2, 1) / r_scale, 40, "green", "filled")
scatter(xy_arrive(1, 1) / r_scale, xy_arrive(2, 1) / r_scale, 40, "green", "filled")
scatter(xy_depart(1, end) / r_scale, xy_depart(2, end) / r_scale, 40, "red", "x")
scatter(xy_arrive(1, end) / r_scale, xy_arrive(2, end) / r_scale, 40, "red", "x")
hold off
title("Synodic Demonstrator")
xlabel("X")
ylabel("Y")
grid on
axis equal

%% Lambert Solver
dV = zeros(N_depart, N_arrive);
for i = 1 : N_depart
    for j = 1 : N_arrive
        % Create space triangles
        r_1_vec = xy_depart(:, i);
        r_2_vec = xy_arrive(:, j);
        r_1 = norm(r_1_vec);
        r_2 = norm(r_2_vec);
        c = vecnorm(r_1_vec - r_2_vec);
        s = (r_1 + r_2 + c) / 2;

        % Compute lambert solutions
        a_lambertsols = lambertSolverSMA(t_arrive(j) - t_depart(i), c, s, mu_sun);

        valid = false;

        % Iterate over Lambert transfer types to find min dV
        dV_total = ones([4, 6]) * 1e9;
        for l = 1 : 6
            a_trans = a_lambertsols{2, l};

            if a_trans > 1e11
                continue
            end

            p_solution = lambertSolverP(a_trans, c, s, r_1, r_2);

            for p = 1 : 4
                p_trans = p_solution{2, p};

                if real(p_trans) > 1 && imag(p_trans) == 0
                    % Compute specific angular momentum
                    h_1 = sqrt(mu_sun * p_1);
                    h_trans = sqrt(mu_sun * p_trans);
                    h_2 = sqrt(mu_sun * p_2);
        
                    % dV_minus
                    v_minus_1 = vis_viva(initial_orbit.a, r_1, mu_sun);
                    gamma_1_minus = flightpath_angle(h_1, r_1, v_minus_1);
                    v_plus_1 = vis_viva(a_trans, r_1, mu_sun);
                    gamma_1_plus = flightpath_angle(h_1, r_1, v_plus_1);
                    dv_1 = law_of_cosines(v_minus_1, v_plus_1, gamma_1_plus - gamma_1_minus);
                    
                    % dV_plus
                    v_minus_2 = vis_viva(a_trans, r_2, mu_sun);
                    gamma_2_minus = flightpath_angle(h_2, r_2, v_minus_2);
                    v_plus_2 = vis_viva(target_orbit.a, r_2, mu_sun);
                    gamma_2_plus = flightpath_angle(h_2, r_2, v_plus_2);
                    dv_2 = law_of_cosines(v_minus_2, v_plus_2, gamma_2_plus - gamma_2_minus);
            
                    % Record plus
                    dV_total_p = dv_1 + dv_2;

                    if imag(dV_total_p) == 0
                        dV_total(p, l) = dV_total_p;
                        valid = true;
                        fprintf("Sucess\n")
                    end
                end
            end
        end
        
        % dV_total
        if valid
            [min_dV, min_i_flat] = min(dV_total, [], "all");
            best_dV(i, j) = min_dV;
            min_i = mod(min_i_flat - 1, 6) + 1;
            best_type(i, j) = {a_lambertsols{1, min_i}};
        end

        fprintf("Completed: %g, %g\n", i, j)
    end
end

%% Plot Porkchop
X = repmat(t_depart, N_arrive, 1);
Y = repmat(t_arrive', 1, N_depart);

Z = best_dV;
%Z = log(best_dV + 1);
Z(Z == 0) = nan;
Z(Z > exp(5)) = nan;

figure
histogram(Z(:))

figure
pcolor(X, Y, Z, FaceColor="flat")
colorbar

% 
% zmin = floor(min(Z(:))); 
% zmax = ceil(max(Z(:)));
% zinc = (zmax - zmin) / 10;
% zlevs = zmin:zinc:zmax;
% 
% figure
% contour(Z,zlevs)

%% Helper Functions
function [r] = orbit_equation(orbit_struct)
    r = orbit_struct.a * (1 - orbit_struct.e ^ 2) / (1 + orbit_struct.e * cos(orbit_struct.thetastar - orbit_struct.w));
end

function [P] = period(orbit_struct, mu)
    P = 2 * pi * sqrt(orbit_struct.a ^ 3 / mu);
end

function [v] = vis_viva(a, r, mu)
    v = sqrt(mu * (2 / r - 1 / a));
end

function [gamma] = flightpath_angle(h, r, v)
    gamma = acosd(h / (r * v));
end

function [s3] = law_of_cosines(s1, s2, angle)
    s3 = sqrt(s1 ^ 2 + s2 ^ 2 - 2 * s1 * s2 * cosd(angle));
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