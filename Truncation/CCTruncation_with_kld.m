function [mu, Sigma, KLD] = CCTruncation_with_kld(mu_0, Sigma_0, a, b, z)
%CCTRUNCATION_WITH_KLD Single halfspace CC truncation with a precomputed z.

    q0 = a' * Sigma_0 * a;
    sqrt_q0 = sqrt(q0);

    margin0 = a' * mu_0 - b - z * sqrt_q0;
    if margin0 >= 0
        mu = mu_0;
        Sigma = Sigma_0;
        KLD = 0;
        return;
    end

    bt = (b - a' * mu_0) / sqrt_q0;
    c0 = z;

    sigma = sqrt(c0^2 * bt^2 + 4 * c0^2 + 4) - c0 * bt;
    sigma = sigma / (2 * (c0^2 + 1));
    lambda = -bt - c0 * sigma;

    mu = mu_0 - Sigma_0 * a / sqrt_q0 * lambda;
    Sigma = Sigma_0 - (1 - sigma^2) * (Sigma_0 * a * a' * Sigma_0) / q0;

    sigma = max(sigma, 1e-12);

    KLD = 0.5 * (lambda^2 + sigma^2 - 1 - 2 * log(sigma));
end
