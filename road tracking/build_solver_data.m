function data = build_solver_data(n, d, varargin)
%BUILD_SOLVER_DATA Build reusable MOSEK template data for fixed dimensions.

if nargin >= 3
    eps_diag = varargin{1};
else
    eps_diag = 1e-6;
end

[~, symbcon_raw] = mosekopt('symbcon echo(0)');
symbcon = symbcon_raw.symbcon;

n_L_elem   = n * (n + 1) / 2;
n_mu       = n;
n_s        = d;
n_t_logdet = n;
n_ttrace   = 1;
n_tquad    = 1;

n_x = n_mu + n_L_elem + n_s + n_t_logdet + n_ttrace + n_tquad;

idx_mu     = 1 : n_mu;
idx_L      = n_mu + (1 : n_L_elem);
idx_s      = n_mu + n_L_elem + (1 : n_s);
idx_t      = n_mu + n_L_elem + n_s + (1 : n_t_logdet);
idx_ttrace = n_mu + n_L_elem + n_s + n_t_logdet + 1;
idx_tquad  = n_mu + n_L_elem + n_s + n_t_logdet + 2;

get_L_idx = @(row, col) row .* (row - 1) / 2 + col;
L_diag_idx = arrayfun(@(k) get_L_idx(k, k), 1:n);

L_indices_map = cell(n, 1);
for col = 1:n
    rows = col:n;
    L_indices_map{col} = idx_L(get_L_idx(rows, col));
end

c = zeros(n_x, 1);
c(idx_ttrace) = 0.5;
c(idx_tquad)  = 0.5;
c(idx_t)      = -1.0;

blx = -inf(n_x, 1);
bux =  inf(n_x, 1);

blx(idx_s) = 0;
blx(idx_ttrace) = 0;
blx(idx_tquad) = 0;

for k = 1:n
    blx(idx_L(L_diag_idx(k))) = eps_diag;
end

nk_soc   = n + 1;
nk_trace = 2 + n * n;
nk_quad  = 2 + n;

accs_exp   = repmat([symbcon.MSK_DOMAIN_PRIMAL_EXP_CONE, 3], 1, n);
accs_soc   = repmat([symbcon.MSK_DOMAIN_QUADRATIC_CONE, nk_soc], 1, d);
accs_trace = [symbcon.MSK_DOMAIN_RQUADRATIC_CONE, nk_trace];
accs_quad  = [symbcon.MSK_DOMAIN_RQUADRATIC_CONE, nk_quad];
accs = [accs_exp, accs_soc, accs_trace, accs_quad];

F_exp = zeros(3 * n, n_x);
g_exp = zeros(3 * n, 1);

rows_x0 = 1 : 3 : (3 * n);
cols_L_diag = idx_L(L_diag_idx);
lin_idx_L = (cols_L_diag - 1) * (3 * n) + rows_x0;
F_exp(lin_idx_L) = 1;

rows_x1 = 2 : 3 : (3 * n);
g_exp(rows_x1) = 1.0;

rows_x2 = 3 : 3 : (3 * n);
cols_t = idx_t;
lin_idx_t = (cols_t - 1) * (3 * n) + rows_x2;
F_exp(lin_idx_t) = 1;

Ia_mu = repmat((1:d)', n, 1);
Ja_mu = kron(idx_mu(:), ones(d, 1));

Ia_s = (1:d)';
Ja_s = idx_s(:);

a_pattern = struct();
a_pattern.I = [Ia_mu; Ia_s];
a_pattern.J = [Ja_mu; Ja_s];
a_pattern.num_mu_entries = numel(Ia_mu);
a_pattern.num_s_entries = numel(Ia_s);

F_soc_template = cell(d, 1);
g_soc_template = cell(d, 1);

for i = 1:d
    Fi = zeros(nk_soc, n_x);
    gi = zeros(nk_soc, 1);

    Fi(1, idx_s(i)) = 1;

    F_soc_template{i} = Fi;
    g_soc_template{i} = gi;
end

template = struct();
template.F_exp = F_exp;
template.g_exp = g_exp;
template.accs_exp = accs_exp;

template.F_soc_template = F_soc_template;
template.g_soc_template = g_soc_template;
template.accs_soc = accs_soc;

template.nk_soc = nk_soc;
template.nk_trace = nk_trace;
template.nk_quad = nk_quad;

template.a_pattern = a_pattern;

index = struct();
index.n = n;
index.d = d;
index.n_x = n_x;
index.n_L_elem = n_L_elem;

index.idx_mu = idx_mu;
index.idx_L = idx_L;
index.idx_s = idx_s;
index.idx_t = idx_t;
index.idx_ttrace = idx_ttrace;
index.idx_tquad = idx_tquad;

index.L_diag_idx = L_diag_idx;
index.get_L_idx = get_L_idx;
index.L_indices_map = L_indices_map;

constant = struct();
constant.c = c;
constant.blx = blx;
constant.bux = bux;
constant.accs = accs;
constant.eps_diag = eps_diag;

meta = struct();
meta.symbcon = symbcon;

cache = struct();
cache.last_sigma_key = [];
cache.last_Sinv = [];
cache.last_Sinv_sqrt = [];

data = struct();
data.index = index;
data.constant = constant;
data.meta = meta;
data.template = template;
data.cache = cache;
end
