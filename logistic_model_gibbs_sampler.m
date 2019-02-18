function samples = logistic_model_gibbs_sampler(n_samples, n_skip, n_delta, times, values, k_bar, X0_low, X0_bar, r_low, r_bar, kappa, Psi, k_init, X0_init, r_init, sig2_init)

k = k_init;
X0 = X0_init;
r = r_init;
sig2 = sig2_init;

samples = zeros(n_samples, 4);

for i = 1:n_skip+1
    r = sample_r(times, values, k, X0, sig2, r_low, r_bar);
    X0 = sample_X0(times, values, k, r, sig2, X0_low, X0_bar);
    k = sample_k(times, values, X0, r, sig2, k_bar);    
    sig2 = sample_sig2(times, values, k, X0, r, kappa, Psi);
end
samples(1, :) = [k X0 r sig2];

for j = 2:n_samples
    for i = 1:n_delta
        r = sample_r(times, values, k, X0, sig2, r_low, r_bar);
        X0 = sample_X0(times, values, k, r, sig2, X0_low, X0_bar);
        k = sample_k(times, values, X0, r, sig2, k_bar);    
        sig2 = sample_sig2(times, values, k, X0, r, kappa, Psi);
    end
samples(j, :) = [k X0 r sig2];
end

end

