function [nu] = hyperbolic_to_true_anomaly(H, e)
%HYPERBOLIC_TO_TRUE_ANOMALY Summary of this function goes here
%   Detailed explanation goes here

nu = 2 * atan(sqrt((e + 1) / (e - 1)) * tanh(H / 2));

end

