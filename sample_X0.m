function X0 = sample_X0(times, values, k, r, sig2, X0_low, X0_bar)

A = values .* exp(-r * times);
B = values .* (1 - exp(-r * times)) * exp(-k) - 1;

alpha = A' * A;
beta = 2 * (A' * B);

sig2prim = sig2 / alpha;
muprim = - beta / (2 * alpha);

min_X0 = X0_low;
max_X0 = X0_bar;
n_grid_deltas = 1000;
delta_X0 = (max_X0 - min_X0) / n_grid_deltas;
grid_points = min_X0:delta_X0:max_X0;
C = -(1./grid_points-muprim).^2 / (2*sig2prim);
C = C - max(C);
grid_values = exp(C);

X0 = sample_continuous(grid_points, grid_values);

end
