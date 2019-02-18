function sig2 = sample_sig2(times, values, k, X0, r, kappa, Psi)
    
n = length(times);
epsilon = values .* ((1 / X0 - exp(-k)) * exp(- r * times) + exp(-k))  - 1;

kappa_p = kappa + n;
Psi_p = Psi + epsilon' * epsilon;
    
sig2 = iwishrnd(Psi_p, kappa_p);

end