AU = 149597898; % [km]
mu_sun = 132712440017.99; % [km3 / s2]


a_kuat = 0.6;
a_vulcan = 2.1;
a_endor = 3.9;
a_tatooine = 4.5;

mu_kuat = 1.2748e6;
mu_vulcan = 2.8303e7;
mu_endor = 1.3957e6;
mu_tatooine = 3.2e5;

initial_orbit.a = a_kuat * AU; % [km]
initial_orbit.e = 0.5;
initial_orbit.thetastar0 = 0; % [rad]
initial_orbit.w = 0.6; % [rad]
initial_orbit.mu = mu_kuat; % [km3 / s2]

target_orbit.a = 1 * AU; % [km]
target_orbit.e = 0.01673163;
target_orbit.thetastar0 = 1; % [rad]
target_orbit.w = 0; % [rad]
target_orbit.mu = mu_endor; % [km3 / s2]

N_depart = 5e2;
N_arrive = 5e2;

% Compute semilatus rectum
p_1 = initial_orbit.a * (1 - initial_orbit.e ^ 2);
p_2 = target_orbit.a * (1 - target_orbit.e ^ 2);

%%
P_0 = period(initial_orbit, mu_sun);
P_f = period(target_orbit, mu_sun);

P_syn = 1 / abs(1 / P_0 - 1 / P_f); % Synodic period - convenient way to plot

multiplier = 2;
t_depart = linspace(0, P_syn, N_depart) * multiplier;
t_arrive = linspace(P_syn, 2 * P_syn, N_arrive) * multiplier + P_syn * 0.1;

thetastar_depart = time_to_thetastar(t_depart, initial_orbit, mu_sun);
thetastar_arrive = time_to_thetastar(t_arrive, target_orbit, mu_sun);

[xy_depart, v_depart] = kepler2D_to_cartestian(initial_orbit, thetastar_depart, mu_sun);
[xy_arrive, v_arrive_3m] = kepler2D_to_cartestian(target_orbit, thetastar_arrive, mu_sun);

%% Lambert Solver

%enter path of the the dll directory with all required files including .bin (with slash at end) 
dllDirectory_Path = convertStringsToChars(string(cd) + "\Helper Functions\Lambert Solvers\ivLamV2p41_738416p65617\matlabInterface\lib\");  %at distribution in this file near the driver, otherwise change here.

addpath(dllDirectory_Path) %add the path where the .dll resides

%load the dll and initialize the lambert routines
iflag=ivLam_initializeDLL(dllDirectory_Path);
if(iflag~=0)
    return
else
    disp('coef path and dll path appear correct, data loaded ok!')
end

[X1, Y1] = meshgrid(xy_depart(1, :), xy_arrive(1, :));
[X2, Y2] = meshgrid(xy_depart(2, :), xy_arrive(2, :));
[t1, t2] = meshgrid(t_depart, t_arrive);
tof = t2 - t1;
[Vd1, Va1] = meshgrid(v_depart(1, :), v_arrive_3m(1, :));
[Vd2, Va2] = meshgrid(v_depart(2, :), v_arrive_3m(2, :));

mu = 1;
mu_star = mu_sun;
l_star = AU;
t_star = sqrt(l_star ^ 3 / mu_star);
v_star = l_star / t_star;

z = zeros([N_depart * N_arrive, 1]);

[v1vec0_pos,v2vec0_pos,~,~] = ivLam_zeroRev_multipleInputDLL(N_depart * N_arrive, [X1(:), X2(:), z]' / l_star, [Y1(:), Y2(:), z]' / l_star, tof(:) / t_star, ones([N_depart * N_arrive, 1]));
[v1vec0_neg,v2vec0_neg,~,~] = ivLam_zeroRev_multipleInputDLL(N_depart * N_arrive, [X1(:), X2(:), z]' / l_star, [Y1(:), Y2(:), z]' / l_star, tof(:) / t_star, -ones([N_depart * N_arrive, 1]));

dV_pos = reshape(vecnorm([Va1(:), Va2(:), z]' / v_star - v2vec0_pos, 2, 1) + vecnorm([Vd1(:), Vd2(:), z]' / v_star - v1vec0_pos, 2, 1), N_depart, N_arrive) * v_star;
dV_neg = reshape(vecnorm([Va1(:), Va2(:), z]' / v_star - v2vec0_neg, 2, 1) + vecnorm([Vd1(:), Vd2(:), z]' / v_star - v1vec0_neg, 2, 1), N_depart, N_arrive) * v_star;
dV_best = min(dV_pos, dV_neg);
pos_sign_choice = dV_best == dV_pos;

% Get best transfer
[min_dV, min_dV_i] = min(dV_best, [], "all");

r_depart_best = xy_depart(:, ceil(min_dV_i / N_arrive));
r_arrive_best = xy_arrive(:, mod(min_dV_i - 1, N_depart) + 1);

t_depart_best = t_depart(ceil(min_dV_i / N_arrive));
t_arrive_best = t_arrive(mod(min_dV_i - 1, N_depart) + 1);
tof_best = t_arrive_best - t_depart_best;

sign_best = pos_sign_choice(min_dV_i);
v1vec0_best = v1vec0_pos(1:2, min_dV_i) * pos_sign_choice(min_dV_i) + v1vec0_neg(1:2, min_dV_i) * ~pos_sign_choice(min_dV_i);
v2vec0_best = v2vec0_pos(1:2, min_dV_i) * pos_sign_choice(min_dV_i) + v2vec0_neg(1:2, min_dV_i) * ~pos_sign_choice(min_dV_i);

[a_dep, e_dep, thetastar_dep, w_dep] = cartesian_to_kepler_2D(r_depart_best, v1vec0_best * v_star, mu_sun);
[a_arr, e_arr, thetastar_arr, w_arr] = cartesian_to_kepler_2D(r_arrive_best, v2vec0_best * v_star, mu_sun);

% Double check the correct transfer was selected - it was 
[v1vecA,v2vecA,~,~,~] = ivLam_singleN_withDetailsDLL([r_depart_best; 0]/l_star,[r_arrive_best; 0]/l_star,tof_best / t_star,2 * pos_sign_choice(min_dV_i) - 1,0);
% 
% % Create transfer with 3 months time of flight
% three_months = 60 * 60 * 24 * 30 * 5;
% 
% thetastar_3m = time_to_thetastar(three_months, target_orbit, mu_sun);
% 
% [xy_arrive_3m, v_arrive_3m] = kepler2D_to_cartestian(target_orbit, thetastar_3m, mu_sun);
% [v1vec_3m,v2vec_3m,~,~,~] = ivLam_singleN_withDetailsDLL([r_depart_best; 0]/l_star,[xy_arrive_3m; 0]/l_star,three_months / t_star, 1, 0);
% [a_3m_dep, e_3m_dep, thetastar_3m_dep, w_3m_dep] = cartesian_to_kepler_2D(r_depart_best, v1vec_3m(1:2) * v_star, mu_sun);
% % want final theta too, just be lazy and call function at final position
% [a_3m_arr, e_3m_arr, thetastar_3m_arr, w_3m_arr] = cartesian_to_kepler_2D(xy_arrive_3m, v2vec_3m(1:2) * v_star, mu_sun);

%unload the dll and clear memory from the lambert routines
iflag= ivLam_unloadDataDLL();

%% Plot Porkchop
X = repmat(t_depart, N_arrive, 1);
Y = repmat(t_arrive', 1, N_depart);

Z = dV_best;
%Z(Z > prctile(Z,95,"all")) = nan; % Get rid of highest 5% of data


figure
pcolor(X / P_f, Y / P_f, Z, FaceColor="interp",EdgeAlpha=0,HandleVisibility="off"); hold on
scatter(t_depart_best / P_f, t_arrive_best / P_f, 40, "red", "filled", DisplayName="Best")
cb = colorbar();
cb.Label.String = 'Total Relative Velocity';

legend()
xlabel("Departure Time [Years]")
ylabel("Arrival Time [Years]")
title("Kuat to Earth (Kuat with 0.5 eccentricity)")
% 
% figure
% pcolor(pos_sign_choice, EdgeAlpha = 0)

% zmin = floor(min(Z(:))); 
% zmax = ceil(max(Z(:)));
% zinc = (zmax - zmin) / 10;
% zlevs = zmin:zinc:zmax;
% 
% figure
% contour(Z,zlevs, "LineWidth",1)

%% Plot best transfer
figure
r_scale = AU;

plot(xy_depart(1, :) / r_scale, xy_depart(2, :) / r_scale, LineStyle="-", DisplayName = "Kuat"); hold on
plot(xy_arrive(1, :) / r_scale, xy_arrive(2, :) / r_scale, LineStyle="-", DisplayName = "Endor");
scatter(r_depart_best(1) / r_scale, r_depart_best(2) / r_scale, 40, "blue", "o")
scatter(r_arrive_best(1) / r_scale, r_arrive_best(2) / r_scale, 40, "blue", "o")
w_offset = 3*pi/2 * ~sign_best + pi * sign_best; % don't know what offset should be when positive sign (also why offset needed?)
orbitplot2D(a_dep, e_arr, linspace(thetastar_dep, thetastar_arr, 100), w_dep + w_offset, "Optimal Transfer", LineStyle = "--", r_scale = r_scale)
hold off
title("Optimal Transfer")
xlabel("X")
ylabel("Y")
grid on
axis equal

%% Helper Functions
function [r] = orbit_equation(orbit_struct, thetastar)
    r = orbit_struct.a * (1 - orbit_struct.e ^ 2) ./ (1 + orbit_struct.e * cos(thetastar));
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

function [r_xy, v_xy] = kepler2D_to_cartestian(orbit_struct, thetastar, mu)
    r = orbit_equation(orbit_struct, thetastar);
    x = r .* cos(thetastar + orbit_struct.w);
    y = r .* sin(thetastar + orbit_struct.w);

    r_xy = [x; y];

    v_xy = zeros([2, numel(thetastar)]);
    for i = 1 : numel(thetastar)
        DCM_RTN_cart = [cos(orbit_struct.w), -sin(orbit_struct.w); sin(orbit_struct.w), cos(orbit_struct.w)];
        
        p = orbit_struct.a * (1 - orbit_struct.e ^ 2);
        hmag = sqrt(p * mu);
        v_RTN = [-(mu / hmag) * sin(thetastar(i)); (mu / hmag) * (orbit_struct.e + cos(thetastar(i)))];

        v_xy(:, i) = DCM_RTN_cart * v_RTN;
    end
end

function [thetastar] = time_to_thetastar(time_period, orbit_struct, mu)
    E_0 = true_to_eccentric_anomaly(orbit_struct.thetastar0, orbit_struct.e);
    time_since_periapsis_0 = sqrt(orbit_struct.a ^ 3 / mu) * (E_0 - orbit_struct.e * sin(E_0));
    mean_anomaly = sqrt(mu / orbit_struct.a ^ 3) * (time_period - time_since_periapsis_0);
    
    thetastar = zeros([1, numel(time_period)]);
    for i = 1 : numel(time_period)
        thetastar(i) = eccentric_to_true_anomaly(mean_to_eccentric_anomaly(mean_anomaly(i), orbit_struct.e), orbit_struct.e);
    end
end

function [a, e, thetastar, w] = cartesian_to_kepler_2D(rvec_2D, vvec_2D, mu)
    rvec = [rvec_2D; 0];
    vvec = [vvec_2D; 0];

    r = norm(rvec);
    v = norm(vvec);

    a = -mu / 2 / (v^2 / 2 - mu / r);

    hvec = cross(rvec, vvec);
    h = norm(hvec);
    evec = 1 / mu * ((v ^ 2 - mu / r) * rvec - dot(rvec, vvec) * vvec);
    e = norm(evec);
    w_sign = 1 - 2 * (evec(3) < 0);
    nvec = cross([1; 0; 0], hvec);
    n = norm(nvec);
    w = acos(dot(nvec, evec) / (n * e)) * w_sign; % QUADRANT? If e(3) > 0 then < pi
    thetastar_sign = 1 - 2 * (dot(rvec, vvec) < 0);
    thetastar = acos(dot(evec, rvec) / (e * r)) * thetastar_sign; % QUADRANT? If dot(r, v) > 0 then < pi
end