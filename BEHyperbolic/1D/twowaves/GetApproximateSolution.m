function result = GetApproximateSolution(type, x, t, k, a1, a2, a12)
    
    if(strcmp(type,'pq') == 1)
        p = x;
        q = t;
        result = -0.2e1 * (k ^ 2) * a12 * (-1 + k ^ 2) * (0.4e1 * a1 ^ 2 * a2 * (k ^ 2) * exp(p - q) -...
            a1 ^ 2 * a2 * exp(p - q) - a2 * a12 ^ 2 * exp(p + q) + a2 * a12 ^ 2 * (k ^ 2) * exp(p + q) - 0.4e1 * a1 * a2 * a12 * exp(p) +...
            0.4e1 * a1 * a2 * a12 * (k ^ 2) * exp(p) - a1 * a12 ^ 2 * exp(0.2e1 * p) + a1 * a12 ^ 2 * (k ^ 2) * exp(0.2e1 * p) +...
            0.4e1 * a1 * a2 ^ 2 * (k ^ 2) - a1 * a2 ^ 2) .* exp(p) ./...
            (0.4e1 * a1 * a2 * (k ^ 2) * exp(p / 0.2e1 - q / 0.2e1) - a1 * a2 * exp(p / 0.2e1 - q / 0.2e1) -...
            a2 * a12 * exp(p / 0.2e1 + q / 0.2e1) + a2 * a12 * (k ^ 2) * exp(p / 0.2e1 + q / 0.2e1) -...
            a1 * a12 * exp(0.3e1 / 0.2e1 * p - q / 0.2e1) + a1 * a12 * (k ^ 2) * exp(0.3e1 / 0.2e1 * p - q / 0.2e1) -...
            a12 ^ 2 * exp(0.3e1 / 0.2e1 * p + q / 0.2e1) + a12 ^ 2 * (k ^ 2) * exp(0.3e1 / 0.2e1 * p + q / 0.2e1)) .^ 2;
    elseif(strcmp(type, 'xt') == 1 || strcmp(type, 'xte') == 1)
        b = 0;
        b = .5171304410;
        result = -0.2e1 * (a1 * (k ^ 2) * exp((k * x) + sqrt((k ^ 2 * (1 - k ^ 2))) * t + b) + a2 * (k ^ 2) * exp(-(k * x) +...
            sqrt((k ^ 2 * (1 - k ^ 2))) * t + b)) ./ ((0.4e1 * a1 * a2 * (k ^ 2) - a1 * a2) / (-a12 + a12 * k ^ 2) + a1 * exp((k * x) +...
            sqrt((k ^ 2 * (1 - k ^ 2))) * t + b) + a2 * exp(-(k * x) + sqrt((k ^ 2 * (1 - k ^ 2))) * t + b) +...
            a12 * exp(0.2e1 * sqrt((k ^ 2 * (1 - k ^ 2))) * t + 0.2e1 * b)) +...
            0.2e1 * (a1 * k * exp((k * x) + sqrt((k ^ 2 * (1 - k ^ 2))) * t + b) - a2 * k * exp(-(k * x) + sqrt((k ^ 2 * (1 - k ^ 2))) * t + b)) .^ 2 ./...
            ((0.4e1 * a1 * a2 * (k ^ 2) - a1 * a2) / (-a12 + a12 * k ^ 2) + a1 * exp((k * x) + sqrt((k ^ 2 * (1 - k ^ 2))) * t + b) +...
            a2 * exp(-(k * x) + sqrt((k ^ 2 * (1 - k ^ 2))) * t + b) + a12 * exp(0.2e1 * sqrt((k ^ 2 * (1 - k ^ 2))) * t + 0.2e1 * b)) .^ 2;
    else
        error('unknown type');
    end
    
end