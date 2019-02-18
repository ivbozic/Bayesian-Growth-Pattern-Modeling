function k = sample_k(times, values, X0, r, sig2, k_bar)

A = values .* (1 - exp(-r * times));
B = values .* (exp(-r * times) / X0) - 1;

alpha = A' * A;
beta = 2 * (A' * B);


sig2prim = sig2 / alpha;
muprim = - beta / (2 * alpha);

min_k = 0;
max_k = k_bar;
n_grid_deltas = 1000;
delta_k = (max_k - min_k) / n_grid_deltas;
grid_points = min_k:delta_k:max_k;
C = -(exp(-grid_points)-muprim).^2 / (2*sig2prim);
C = C - max(C);
grid_values = exp(C);

k = sample_continuous(grid_points, grid_values);

end

