function [mu, Sigma, obj] = PolCCT(mu_0, Sigma_0, A, B, z, mosek_data_list)
%POLCCT Solve the proposed chance-constrained Gaussian truncation subproblem.

    if isempty(A)
        mu = mu_0;
        Sigma = Sigma_0;
        obj = 0;
        return;
    end

    if nargin == 5
        mosek_data_list = [];
    end

    z_i = z(1);
    [n, d] = size(A);

    if isempty(mosek_data_list)
        data = build_solver_data(n, d);
    else
        idx = find([mosek_data_list.d] == d, 1);
        if isempty(idx)
            data = build_solver_data(n, d);
        else
            data = mosek_data_list(idx).data;
        end
    end

    [mu, Sigma, solve_info, ~] = solve_KL_chance_mosek_final( ...
        data, mu_0, Sigma_0, A, B, z_i);

    if ~strcmp(string(solve_info.status), "OPTIMAL")
        obj = -1;
    else
        obj = gaussKLD(mu, mu_0, Sigma, Sigma_0);
    end
end

function [mu, Sigma, solve_info, data] = solve_KL_chance_mosek_final( ...
    data, mu_0, Sigma_0, A, B, z_i)
%SOLVE_KL_CHANCE_MOSEK_FINAL Assemble and solve the MOSEK conic problem.

    [prob, param, data] = assemble_prob(data, A, B, mu_0, Sigma_0, z_i);
    [rcode, res] = mosekopt('minimize echo(0)', prob, param);

    idx_mu = data.index.idx_mu;
    idx_L  = data.index.idx_L;
    n = data.index.n;
    get_L_idx = data.index.get_L_idx;

    xx = res.sol.itr.xx;
    mu = xx(idx_mu);

    L_lower = xx(idx_L);
    L = zeros(n, n);
    for row = 1:n
        for col = 1:row
            L(row, col) = L_lower(get_L_idx(row, col));
        end
    end
    Sigma = L * L';

    solve_info.status = res.sol.itr.solsta;
    solve_info.pobj = res.sol.itr.pobjval - n / 2;
    solve_info.rescode = rcode;
end
