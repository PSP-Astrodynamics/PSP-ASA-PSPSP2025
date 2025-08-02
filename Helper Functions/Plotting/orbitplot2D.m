function [] = orbitplot2D(a, e, thetastar, w, Name, options)
%ORBITPLOT2D Plot 2D orbit
%   Detailed explanation goes here
arguments
    a % semimajor axis
    e % eccentricity
    thetastar % true anomaly vector (row vector)
    w % argument of periapsis
    Name % name of orbit
    options.r_scale = 1
    options.LineStyle = "-"
end

xy = kepler2D_to_cartestian(a, e, thetastar, w);
plot(xy(1, :) / options.r_scale, xy(2, :) / options.r_scale, DisplayName = Name, LineStyle = options.LineStyle); hold on

end

function [xy] = kepler2D_to_cartestian(a, e, thetastar, w)
    p = a .* (1 - e .^ 2);
    r = p ./ (1 + e * cos(thetastar));
    x = r .* cos(thetastar + w);
    y = r .* sin(thetastar + w);

    xy = [x; y];
end