function [prob, param, data] = assemble_prob(data, A, B, mu_0, Sigma_0, z_i)
%ASSEMBLE_PROB Assemble a MOSEK problem from cached template data.

n = data.index.n;
d = data.index.d;
n_x = data.index.n_x;

idx_mu = data.index.idx_mu;
idx_s = data.index.idx_s;
idx_ttrace = data.index.idx_ttrace;
idx_tquad = data.index.idx_tquad;
L_indices_map = data.index.L_indices_map;

if ~isequal(size(A), [n, d])
    error('A must be n-by-d.');
end
if ~isequal(size(Sigma_0), [n, n])
    error('Sigma_0 must be n-by-n.');
end
if numel(mu_0) ~= n
    error('mu_0 must have length n.');
end

if isrow(B), B = B(:); end
if isrow(mu_0), mu_0 = mu_0(:); end
if ~isscalar(z_i) && isrow(z_i), z_i = z_i(:); end

[Sinv, Sinv_sqrt, data] = local_get_sigma_factors(data, Sigma_0);

prob = struct();

prob.c = data.constant.c;
prob.blx = data.constant.blx;
prob.bux = data.constant.bux;
prob.accs = data.constant.accs;

a_pat = data.template.a_pattern;

Va_mu = reshape(A', [], 1);
if isscalar(z_i)
    Va_s = -repmat(z_i, d, 1);
else
    Va_s = -z_i(:);
end

prob.a = sparse( ...
    a_pat.I, ...
    a_pat.J, ...
    [Va_mu; Va_s], ...
    d, n_x);

prob.blc = B;

nk_soc   = data.template.nk_soc;
nk_trace = data.template.nk_trace;
nk_quad  = data.template.nk_quad;

n_rows_exp   = 3 * n;
n_rows_soc   = d * nk_soc;
n_rows_trace = nk_trace;
n_rows_quad  = nk_quad;

nF = n_rows_exp + n_rows_soc + n_rows_trace + n_rows_quad;

F_dense = zeros(nF, n_x);
g_dense = zeros(nF, 1);

row_start = 1;

rows = row_start : (row_start + n_rows_exp - 1);
F_dense(rows, :) = data.template.F_exp;
g_dense(rows) = data.template.g_exp;
row_start = row_start + n_rows_exp;

F_soc_cell = data.template.F_soc_template;

for i = 1:d
    ai = A(:, i);
    Fi = F_soc_cell{i};

    for r = 1:n
        cols_idx = L_indices_map{r};
        coefs = ai(r:n);
        Fi(r + 1, cols_idx) = coefs';
    end

    rows = row_start : (row_start + nk_soc - 1);
    F_dense(rows, :) = Fi;
    row_start = row_start + nk_soc;
end

Fi_trace = zeros(nk_trace, n_x);
gi_trace = zeros(nk_trace, 1);

Fi_trace(1, idx_ttrace) = 1;
gi_trace(2) = 0.5;

for k = 1:n
    rows_range = 2 + (k - 1) * n + (1:n);
    l_vars_idx = L_indices_map{k};
    block_S = Sinv_sqrt(:, k:n);
    Fi_trace(rows_range, l_vars_idx) = block_S;
end

rows = row_start : (row_start + nk_trace - 1);
F_dense(rows, :) = Fi_trace;
g_dense(rows) = gi_trace;
row_start = row_start + nk_trace;

Fi_quad = zeros(nk_quad, n_x);
gi_quad = zeros(nk_quad, 1);

Fi_quad(1, idx_tquad) = 1;
gi_quad(2) = 0.5;
Fi_quad(3:end, idx_mu) = Sinv_sqrt;
gi_quad(3:end) = -(Sinv_sqrt * mu_0(:));

rows = row_start : (row_start + nk_quad - 1);
F_dense(rows, :) = Fi_quad;
g_dense(rows) = gi_quad;
row_start = row_start + nk_quad;

if row_start ~= nF + 1
    error('Internal row counting mismatch in assemble_prob.');
end

prob.f = sparse(F_dense);
prob.g = g_dense;

param = struct();
param.MSK_IPAR_LOG = 0;

prob.userdata.Sinv = Sinv;
prob.userdata.Sinv_sqrt = Sinv_sqrt;
prob.userdata.mu_0 = mu_0;
end

function [Sinv, Sinv_sqrt, data] = local_get_sigma_factors(data, Sigma_0)
% Cache inverse covariance factors for repeated Sigma_0 values.

use_cache = false;
if isfield(data, 'cache') && isfield(data.cache, 'last_sigma_key')
    if ~isempty(data.cache.last_sigma_key) ...
            && isequal(size(data.cache.last_sigma_key), size(Sigma_0)) ...
            && isequal(data.cache.last_sigma_key, Sigma_0)
        use_cache = true;
    end
end

if use_cache
    Sinv = data.cache.last_Sinv;
    Sinv_sqrt = data.cache.last_Sinv_sqrt;
    return;
end

n = size(Sigma_0, 1);

Lc = chol(Sigma_0, 'lower');

Sinv_sqrt = Lc \ eye(n);

Sinv = Sinv_sqrt' * Sinv_sqrt;

data.cache.last_sigma_key = Sigma_0;
data.cache.last_Sinv = Sinv;
data.cache.last_Sinv_sqrt = Sinv_sqrt;

end
