function [] = plot_Hohmann_transfer(thetastar_0, ae_sc, orbit_names, wait_times, delta_v, objects, mu, N, title_text, r_scale, units)
%PLOT_HOHMANN_TRANSFER Summary of this function goes here
%   Detailed explanation goes here

n = numel(wait_times);

a_sc = ae_sc(1, :);
e_sc = ae_sc(2, :);

a_obj = objects.a;
e_obj = objects.e;
w_obj = objects.w;
t_max_obj = objects.t_max;
color_obj = objects.c;
name_obj = objects.name;

t = 0;
thetastar = thetastar_0;
t_adj = 0;

% Plot objects
for o = 1 : numel(name_obj)
    [x_orbit, y_orbit] = r_conic(a_obj(o), e_obj(o), w_obj(o), N);

    % Plot underlying orbit
    plot(x_orbit / r_scale, y_orbit / r_scale, LineStyle="--", HandleVisibility='off'); hold on
    
    % Plot section where object travels
    if t_max_obj(o) == -1
        tspan = [0, sum(wait_times)];
    else
        tspan = [0, t_max_obj(o)];
    end
    if tspan(1) ~= tspan(2)
        [x_travel, y_travel, ~] = r_tspan(a_obj(o), e_obj(o), w_obj(o), tspan, mu, N);

        plot(x_travel / r_scale, y_travel / r_scale, LineStyle="-", Color=color_obj(o), DisplayName = name_obj(o), LineWidth=1.2); hold on
    end

    if t_max_obj(o) == -1 || sum(wait_times) < t_max_obj(o)
        scatter(x_travel(end) / r_scale, y_travel(end) / r_scale, color_obj(o), "diamond", LineWidth = 1.2, HandleVisibility='off'); hold on
    end
end

% Plot spacecraft
for m = 1 : n
    if mod(m, 2) == 0 && abs(a_sc(m)) < abs(a_sc(m - 1))
        thetastar = thetastar + pi;
        t_adj = pi * sqrt(a_sc(m) ^ 3 / mu);
    end

    [x_orbit, y_orbit] = r_conic(a_sc(m), e_sc(m), thetastar, N);

    % Plot underlying orbit
    plot(x_orbit / r_scale, y_orbit / r_scale, LineStyle="--", DisplayName = orbit_names(m)); hold on
    
    % Plot section where spacecraft travels
    tspan = [0, wait_times(m)] + t_adj;
    if tspan(1) ~= tspan(2)
        [x_travel, y_travel, thetastar_traveled] = r_tspan(a_sc(m), e_sc(m), thetastar, tspan, mu, N);

        plot(x_travel / r_scale, y_travel / r_scale, LineStyle="-", Color="k", HandleVisibility='off', LineWidth = 1.5); hold on

        thetastar = thetastar + thetastar_traveled;
    end
    if mod(m, 2) == 0 && abs(a_sc(m)) < abs(a_sc(m - 1))
        thetastar = thetastar - pi;
        t_adj = 0;
    end

    % Increment time
    t = t + wait_times(m);

    if mod(m, 2) == 0
        % into transfer burn
        in_burn = [x_travel(2) - x_travel(1), y_travel(2) - y_travel(1)] / norm([x_travel(2) - x_travel(1), y_travel(2) - y_travel(1)])  * delta_v(m - 1);
        quiver(x_travel(1) / r_scale, y_travel(1) / r_scale, in_burn(1), in_burn(2), 0.5, "filled", Color="r", LineWidth=1.2, HandleVisibility="off", MaxHeadSize=1.5);

        % out of transfer burn
        out_burn = [x_travel(end) - x_travel(end - 1), y_travel(end) - y_travel(end - 1)] / norm([x_travel(end) - x_travel(end - 1), y_travel(end) - y_travel(end - 1)]) * delta_v(m);
        quiver(x_travel(end) / r_scale, y_travel(end) / r_scale, out_burn(1), out_burn(2), 0.5, "filled", Color="r", LineWidth=1.2, HandleVisibility="off", MaxHeadSize=1.5);
    end
end
[x_0, y_0, ~] = r_tspan(a_sc(1), e_sc(1), 0, [0,0], mu, 1);
scatter(x_0 / r_scale, y_0 / r_scale, "green", "o", LineWidth = 1.5, DisplayName="Start"); hold on
scatter(x_travel(end) / r_scale, y_travel(end) / r_scale, "red", "x", LineWidth = 1.5, DisplayName="End"); hold on
legend(Location="eastoutside")
xlabel(sprintf("X [%s]", units))
ylabel(sprintf("Y [%s]", units))
title(title_text)

hold off
axis equal
grid on

% 
% figure
% t = t_cont;
% 
% r1 = ;
% r2 = ;
% 
% r = sqrt(r1.^2 + r2.^2);
% 
% spacecraft = animatedline('LineWidth',2);
% mars = animatedline('lineWidth', 1);
% 
% hold on;
% 
% spacecraft.Color = 'blue';
% spacecraft.Visible = 'on';
% mars.Color = 'red';
% mars.Visible = 'on';
% 
% title(sprintf('Spacecraft 1\nTime: %0.2f sec | Radius: %0.2f', t(1),r(1)), 'Interpreter', 'latex','Color','k');
% 
% j_0 = int16((t(1) / t(length(t)))*(length(t))) + 1; % stating time index (integer)
% 
% for j = j_0:length(t) % this foor loop generates the animated plot with position vectors with indexes
%     addpoints(spacecraft,r1(j),r2(j));
%     addpoints(mars,mars_trajectory(1, j), mars_trajectory(2, j), 0);
%     head = scatter(r1(j),r2(j),'filled','MarkerFaceColor','b','LineWidth', 0.8);
%     mars_head = scatter3(mars_trajectory(1, j), mars_trajectory(2, j), 0, "red", "filled");
% 
%     drawnow
%     pause(0.01);
%     delete(head);
%     delete(mars_head);
% 
%     title(sprintf('Spacecraft 1\nTime: %0.0f sec | Mass: %0.2f', t(j), x_cont(5, j)), 'Interpreter', 'Latex', 'Color', 'k');
%     axis equal
%     grid on
% end
% 
end

function [x, y] = r_conic(a, e, w, N)
    if e < 1 % Ellipse
        thetastar_orbit = linspace(0, 2 * pi, N);
    elseif e >= 1 % Parabolic or Hyperbolic
        thetastar_orbit = linspace(-pi * 1 / 2, pi * 1 / 2, N);
    end

    p = a .* (1 - e .^ 2);
    r_orbit = p ./ (1 + e * cos(thetastar_orbit));

    x = r_orbit .* cos(thetastar_orbit + w);
    y = r_orbit .* sin(thetastar_orbit + w);
end

function [x, y, thetastar_traveled] = r_tspan(a, e, w, tspan, mu, N)
    M = sqrt(mu / abs(a) ^ 3) * linspace(tspan(1), tspan(2), N);
    if e < 1 % Ellipse
        thetastar_travel = eccentric_to_true_anomaly(mean_to_eccentric_anomaly(M, e), e);
    elseif e > 1 % Hyperbolic
        thetastar_travel = hyperbolic_to_true_anomaly(mean_to_hyperbolic_anomaly(M, e), e);
    end

    p = a .* (1 - e .^ 2);
    
    r_orbit = p ./ (1 + e * cos(thetastar_travel));

    x = r_orbit .* cos(thetastar_travel + w);
    y = r_orbit .* sin(thetastar_travel + w);

    thetastar_traveled = thetastar_travel(end);
end