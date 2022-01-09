function result = GetApproximateSolution(p,q,k,a1,a2,a12)
    result = -0.2e1 * (k ^ 2) * a12 * (-1 + k ^ 2) * (0.4e1 * a1 ^ 2 * a2 * (k ^ 2) * exp(p - q) -...
		a1 ^ 2 * a2 * exp(p - q) - a2 * a12 ^ 2 * exp(p + q) + a2 * a12 ^ 2 * (k ^ 2) * exp(p + q) - 0.4e1 * a1 * a2 * a12 * exp(p) +...
		0.4e1 * a1 * a2 * a12 * (k ^ 2) * exp(p) - a1 * a12 ^ 2 * exp(0.2e1 * p) + a1 * a12 ^ 2 * (k ^ 2) * exp(0.2e1 * p) +...
		0.4e1 * a1 * a2 ^ 2 * (k ^ 2) - a1 * a2 ^ 2) * exp(p) ./...
		(0.4e1 * a1 * a2 * (k ^ 2) * exp(p / 0.2e1 - q / 0.2e1) - a1 * a2 * exp(p / 0.2e1 - q / 0.2e1) -...
		a2 * a12 * exp(p / 0.2e1 + q / 0.2e1) + a2 * a12 * (k ^ 2) * exp(p / 0.2e1 + q / 0.2e1) -...
		a1 * a12 * exp(0.3e1 / 0.2e1 * p - q / 0.2e1) + a1 * a12 * (k ^ 2) * exp(0.3e1 / 0.2e1 * p - q / 0.2e1) -...
		a12 ^ 2 * exp(0.3e1 / 0.2e1 * p + q / 0.2e1) + a12 ^ 2 * (k ^ 2) * exp(0.3e1 / 0.2e1 * p + q / 0.2e1)) ^ 2;
end