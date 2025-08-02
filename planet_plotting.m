a_kuat = 0.6;
a_vulcan = 2.1;
a_endor = 3.9;
a_tatooine = 4.5;

mu_kuat = 1.2748e6;
mu_vulcan = 2.8303e7;
mu_endor = 1.3957e6;
mu_tatooine = 3.2e5;

mu_sun = 132712440017.99;

scale_factor = 5;
scale_kuat = scale_factor * log(mu_kuat);
scale_vulcan = scale_factor * log(mu_vulcan);
scale_endor = scale_factor * log(mu_endor);
scale_tatooine = scale_factor * log(mu_tatooine);
scale_sun = scale_factor * log(mu_sun);

plot_conic(a_kuat, 0, deg2rad(320), 100, "Kuat", scale_kuat, "b"); hold on
plot_conic(a_vulcan, 0, deg2rad(100), 100, "Vulcan", scale_vulcan, "g")
plot_conic(a_endor, 0, deg2rad(30), 100, "Endor", scale_endor, "r")
plot_conic(a_tatooine, 0, deg2rad(190), 100, "Tatooine", scale_tatooine, "#FFA500", scatter_color = [])
plot_conic(0, 0, 0, 100, "Sol", scale_sun, "y")
hold off
axis equal
grid on
xlabel("X [AU]")
ylabel("Y [AU]")
legend()

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

function [] = plot_conic(a, e, w, N, name, size, color, options)
arguments
    a 
    e 
    w 
    N 
    name 
    size 
    color 
    options.scatter_color = color
end
    [x_orbit, y_orbit] = r_conic(a, e, w, N);

    plot(x_orbit, y_orbit, LineStyle="-", HandleVisibility='off', Color=color); hold on
    if ~isempty(options.scatter_color)
        scatter(x_orbit(1), y_orbit(1), size, options.scatter_color, "filled", DisplayName = name, MarkerEdgeColor='k')
    else
        scatter(x_orbit(1), y_orbit(1), size, DisplayName = name, MarkerEdgeColor='k')
    end
end