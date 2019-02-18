function r = sample_r(times, values, k, X0, sig2, r_low, r_bar)

A = values * (1 / X0 - exp(-k));
B = values * exp(-k) - 1;

min_r = r_low;
max_r = r_bar;
n_grid_deltas = 10000;
delta_r = (max_r - min_r) / n_grid_deltas;
grid_points = min_r:delta_r:max_r;

C = -sum((A .* exp(- times * grid_points) + B).^2) / (2*sig2);
C = C - max(C);
grid_values = exp(C);

r = sample_continuous(grid_points, grid_values);

end
