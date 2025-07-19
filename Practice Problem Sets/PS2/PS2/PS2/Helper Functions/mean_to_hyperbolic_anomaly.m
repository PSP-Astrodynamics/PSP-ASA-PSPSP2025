function H = mean_to_hyperbolic_anomaly(M, e)
    % Solve Kepler's equation M = e sinh(H) - H
    H(1:numel(M)) = M;
    options = optimoptions('fsolve','Display','none');
    for index = 1:numel(M)
        H(index) = fsolve(@(x) M(index) - (e * sinh(x) - x), M, options);
    end
end