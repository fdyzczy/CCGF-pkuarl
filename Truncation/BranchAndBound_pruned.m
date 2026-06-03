function [mu_star, Sigma_star, J_ub, info] = BranchAndBound_pruned( ...
    mu_0, Sigma_0, Obstacles, risk_thresholds, mosek_data_list, verbose)
%BRANCHANDBOUND_PRUNED Branch-and-bound solver for polygonal CC truncation.

if nargin < 5 || isempty(mosek_data_list)
    mosek_data_list = [];
end
if nargin < 6 || isempty(verbose)
    verbose = false;
end

max_solves = 1000;
N_c = length(Obstacles);

if isscalar(risk_thresholds)
    risk_thresholds = repmat(risk_thresholds / N_c, 1, N_c);
else
    risk_thresholds = risk_thresholds(:).';
end
z_table = norminv(1 - risk_thresholds);

is_single_edge = false(1, N_c);
for j = 1:N_c
    is_single_edge(j) = (size(Obstacles(j).A, 2) == 1);
end
fixed_obs = find(is_single_edge);

for k = 1:length(fixed_obs)
    j = fixed_obs(k);

    a_try = Obstacles(j).A(:,1);
    b_try = Obstacles(j).b(1);
    delta_try = risk_thresholds(j);

    [mu_try, Sigma_try,] = CCTruncation(mu_0, Sigma_0, a_try, b_try, delta_try);
    val_try = gaussKLD(mu_try, mu_0, Sigma_try, Sigma_0);

    [is_global_feasible, valid_pairs] = check_global_feasibility(...
        mu_try, Sigma_try, Obstacles, z_table, []);

    if is_global_feasible
        if verbose
            fprintf('  [Heuristic] Solution from Obstacle %d satisfies ALL constraints! Early exit.\n', j);
        end

        mu_star = mu_try;
        Sigma_star = Sigma_try;
        J_ub = val_try;
        best_constraint_pairs = valid_pairs;

        info.status = 'Solved (Single-Edge Heuristic)';
        info.iterations = 0;
        info.num_solved = k; % Count k solved subproblems.
        info.best_constraint_pairs = best_constraint_pairs;
        return;
    end
end

if verbose
    fprintf('  No single-constraint solution is globally feasible. Proceeding to BnB...\n');
end


if verbose
    fprintf('Running single-edge root setup (%d fixed constraints)...\n', numel(fixed_obs));
end

root_pairs = cell(1, numel(fixed_obs));
for k = 1:numel(fixed_obs)
    root_pairs{k} = [fixed_obs(k), 1];
end

is_fixed = false(1, N_c);
is_fixed(fixed_obs) = true;
unfathomed_polygons = find(~is_fixed);
unfathomed_edges = cell(1, N_c);
for j = 1:N_c
    if is_single_edge(j)
        unfathomed_edges{j} = [];
    else
        unfathomed_edges{j} = 1:size(Obstacles(j).A, 2);
    end
end

root_node.constraint_pairs = root_pairs;
root_node.depth = numel(root_pairs);
[sub_A0, sub_b0] = build_constraints_from_pairs(root_node.constraint_pairs, Obstacles);

if isempty(sub_A0)
    root_node.mu = mu_0;
    root_node.Sigma = Sigma_0;
    root_node.val = 0;
else
    if length(sub_b0) > 1
        [mu_root, Sigma_root, val_root] = PolCCT( ...
            mu_0, Sigma_0, sub_A0, sub_b0, z_table, mosek_data_list);
    else
        obs_idx = root_node.constraint_pairs{1}(1);
        [mu_root, Sigma_root, val_root] = CCTruncation_with_kld( ...
            mu_0, Sigma_0, sub_A0, sub_b0, z_table(obs_idx));
    end

    if val_root == -1
        [mu_star, Sigma_star, J_ub, info] = infeasible_result( ...
            unfathomed_polygons, unfathomed_edges);
        return;
    end

    root_node.mu = mu_root;
    root_node.Sigma = Sigma_root;
    root_node.val = val_root;
end

J_ub = inf;
mu_star = [];
Sigma_star = [];
best_constraint_pairs = {};
Tree = {root_node};

iter_count = 0;
pruned_count = 0;
pruned_by_unfathomed = 0;
score_lookup_count = 0;
solved_count = 0;
status = 'Solved';

KL_table = precompute_single_edge_costs(mu_0, Sigma_0, Obstacles, z_table, verbose);

while ~isempty(Tree)
    depths = cellfun(@(x) x.depth, Tree);
    max_depth = max(depths);
    deepest_indices = find(depths == max_depth);

    min_score = inf;
    selected_idx = -1;
    nodes_to_remove = [];

    for idx = deepest_indices
        node = Tree{idx};
        [score, lookups] = compute_node_score_fast(node, KL_table);
        score_lookup_count = score_lookup_count + lookups;

        if score > J_ub
            nodes_to_remove = [nodes_to_remove, idx]; %#ok<AGROW>
        elseif score < min_score
            min_score = score;
            selected_idx = idx;
        end
    end

    if ~isempty(nodes_to_remove)
        Tree(nodes_to_remove) = [];
        pruned_count = pruned_count + length(nodes_to_remove);

        if verbose
            fprintf('  [Prune by score] Removed %d nodes with score > J_ub\n', ...
                length(nodes_to_remove));
        end

        if selected_idx ~= -1
            selected_idx = selected_idx - sum(nodes_to_remove < selected_idx);
        end
    end

    if selected_idx == -1
        if verbose
            fprintf('  [Deepest layer] All nodes pruned, continue\n');
        end
        continue;
    end

    if verbose && mod(iter_count, 50) == 0
        fprintf('  [Iter %d] Min score in deepest layer: %.6f\n', iter_count, min_score);
    end

    current_node = Tree{selected_idx};
    Tree(selected_idx) = [];
    iter_count = iter_count + 1;
    depth = current_node.depth;
    force_break = false;

    [sub_A, sub_b] = build_constraints_from_pairs(current_node.constraint_pairs, Obstacles);

    if isempty(current_node.mu)
        if isempty(sub_A)
            mu = mu_0;
            Sigma = Sigma_0;
            val = 0;
        elseif length(sub_b) > 1
            [mu, Sigma, val] = PolCCT(mu_0, Sigma_0, sub_A, sub_b, ...
                z_table, mosek_data_list);
        else
            obs_idx = current_node.constraint_pairs{1}(1);
            [mu, Sigma] = CCTruncation( ...
                mu_0, Sigma_0, sub_A, sub_b, risk_thresholds(obs_idx));
            val = gaussKLD(mu, mu_0, Sigma, Sigma_0);
        end

        if val == -1
            pruned_count = pruned_count + 1;
            continue;
        end

        solved_count = solved_count + 1;
        if solved_count >= max_solves
            status = 'MaxSolvesReached';
            force_break = true;
        end
    else
        mu = current_node.mu;
        Sigma = current_node.Sigma;
        val = current_node.val;
    end

    if val >= J_ub
        pruned_count = pruned_count + 1;
        if force_break
            break;
        end
        continue;
    end

    fixed_polygons = cellfun(@(x) x(1), current_node.constraint_pairs);
    is_global_feasible = check_global_feasibility( ...
        mu, Sigma, Obstacles, z_table, fixed_polygons);

    if is_global_feasible
        if val < J_ub
            J_ub = val;
            mu_star = mu;
            Sigma_star = Sigma;
            best_constraint_pairs = current_node.constraint_pairs;

            if verbose
                fprintf('  Update Best! Cost: %.6f (Depth: %d)\n', J_ub, depth);
            end
        end

        if force_break
            break;
        end
        continue;
    end

    if force_break
        break;
    end

    is_fixed = false(1, N_c);
    is_fixed(fixed_polygons) = true;
    available_polygons = find(~is_fixed);

    if isempty(available_polygons)
        continue;
    end

    j_next = select_best_polygon_heuristic( ...
        mu, Sigma, Obstacles, available_polygons, unfathomed_edges, z_table, verbose);
    if isempty(j_next)
        continue;
    end

    available_edges = unfathomed_edges{j_next};
    if isempty(available_edges)
        continue;
    end

    for l = available_edges
        child_pairs = [current_node.constraint_pairs, {[j_next, l]}];
        a_new = Obstacles(j_next).A(:, l);
        b_new = Obstacles(j_next).b(l);
        z_val = z_table(j_next);

        std_dev = sqrt(max(a_new' * Sigma * a_new, 0));
        constraint_margin = a_new' * mu - z_val * std_dev - b_new;

        child_node.constraint_pairs = child_pairs;
        child_node.depth = depth + 1;

        if constraint_margin >= 1e-4
            child_node.mu = mu;
            child_node.Sigma = Sigma;
            child_node.val = val;
        else
            child_node.mu = [];
            child_node.Sigma = [];
            child_node.val = [];
        end

        Tree{end + 1} = child_node; %#ok<AGROW>
    end
end

info.active_constraints = [];
info.best_constraint_pairs = best_constraint_pairs;
info.iterations = iter_count;
info.num_solved = solved_count;
info.pruned = pruned_count;
info.pruned_by_unfathomed = pruned_by_unfathomed;
info.final_unfathomed_polygons = unfathomed_polygons;
info.final_unfathomed_edges = unfathomed_edges;
info.score_lookups = score_lookup_count;
info.status = status;

if isinf(J_ub)
    if verbose
        warning('No feasible solution found.');
    end
    info.status = 'Infeasible';
end

if verbose
    fprintf('\n--- Summary ---\n');
    fprintf('Total iterations: %d\n', iter_count);
    fprintf('Nodes pruned by bound: %d\n', pruned_count);
    fprintf('Nodes pruned by unfathomed set: %d\n', pruned_by_unfathomed);
    fprintf('Final objective: %.6f\n', J_ub);

    if ~isempty(best_constraint_pairs)
        fprintf('\n--- Optimal Path ---\n');
        for k = 1:length(best_constraint_pairs)
            pair = best_constraint_pairs{k};
            fprintf('  (Polygon %d, Edge %d)\n', pair(1), pair(2));
        end
    end
end
end

function [mu_star, Sigma_star, J_ub, info] = infeasible_result( ...
    unfathomed_polygons, unfathomed_edges)
mu_star = [];
Sigma_star = [];
J_ub = inf;
info.status = 'Infeasible';
info.iterations = 0;
info.num_solved = 0;
info.pruned = 0;
info.pruned_by_unfathomed = 0;
info.best_constraint_pairs = {};
info.active_constraints = [];
info.final_unfathomed_polygons = unfathomed_polygons;
info.final_unfathomed_edges = unfathomed_edges;
info.score_lookups = 0;
end

function KL_table = precompute_single_edge_costs(mu_0, Sigma_0, Obstacles, z_table, verbose)
N_c = length(Obstacles);
KL_table = cell(N_c, 1);

if verbose
    fprintf('Precomputing KL divergences for all edges...\n');
end

for j = 1:N_c
    n_edges = size(Obstacles(j).A, 2);
    KL_table{j} = inf(1, n_edges);

    for l = 1:n_edges
        a = Obstacles(j).A(:, l);
        b = Obstacles(j).b(l);
        [~, ~, KL_table{j}(l)] = CCTruncation_with_kld( ...
            mu_0, Sigma_0, a, b, z_table(j));
    end
end
end

function [A_concat, b_concat] = build_constraints_from_pairs(constraint_pairs, Obstacles)
A_concat = [];
b_concat = [];

for i = 1:length(constraint_pairs)
    pair = constraint_pairs{i};
    j = pair(1);
    l = pair(2);

    A_concat = [A_concat, Obstacles(j).A(:, l)]; %#ok<AGROW>
    b_concat = [b_concat; Obstacles(j).b(l)]; %#ok<AGROW>
end
end

function j_best = select_best_polygon_heuristic( ...
    mu_current, Sigma_current, Obstacles, available_polygons, unfathomed_edges, z_table, verbose)
% Pick the polygon with the smallest best-edge chance-constraint margin.

if nargin < 7
    verbose = false;
end

mu_current = mu_current(:);
best_score = inf;
j_best = [];

for j = available_polygons
    edges_j = unfathomed_edges{j};
    if isempty(edges_j)
        continue;
    end

    obs = Obstacles(j);
    A_all = obs.A(:, edges_j);
    b_all = obs.b(edges_j);
    b_all = b_all(:);

    if isscalar(z_table)
        z_val = z_table;
    else
        z_val = z_table(j);
    end

    mean_part = A_all' * mu_current;
    SA = Sigma_current * A_all;
    var_part = sum(A_all .* SA, 1)';
    var_part = max(var_part, 0);

    margins = mean_part - z_val * sqrt(var_part) - b_all;
    score_j = max(margins);

    if score_j < best_score
        best_score = score_j;
        j_best = j;
    end
end

if verbose && ~isempty(j_best)
    fprintf('  [Heuristic] Selected polygon %d with min(max margin) = %.6f\n', ...
        j_best, best_score);
end
end

function [score, lookups] = compute_node_score_fast(node, KL_table)
lookups = 0;

if isempty(node.constraint_pairs)
    score = 0;
    return;
end

max_kl = -inf;
for i = 1:length(node.constraint_pairs)
    pair = node.constraint_pairs{i};
    poly_idx = pair(1);
    edge_idx = pair(2);

    lookups = lookups + 1;
    kl = KL_table{poly_idx}(edge_idx);

    if isfinite(kl) && kl > max_kl
        max_kl = kl;
    end
end

if max_kl == -inf
    score = inf;
else
    score = max_kl;
end
end

function [is_feasible, valid_pairs] = check_global_feasibility( ...
    mu, Sigma, Obstacles, z_table, fixed_polygons)
% A node is globally feasible if every remaining polygon has one feasible edge.

N_c = numel(Obstacles);

if isempty(fixed_polygons)
    remaining_polygons = 1:N_c;
else
    is_fixed = false(1, N_c);
    is_fixed(fixed_polygons) = true;
    remaining_polygons = find(~is_fixed);
end

valid_pairs = cell(1, numel(remaining_polygons));
pair_count = 0;
is_feasible = true;
mu = mu(:);

for t = 1:numel(remaining_polygons)
    obs_idx = remaining_polygons(t);
    obs = Obstacles(obs_idx);

    if isscalar(z_table)
        z_val = z_table;
    else
        z_val = z_table(obs_idx);
    end

    Aobs = obs.A;
    bobs = obs.b(:);

    mean_part = Aobs' * mu;
    SA = Sigma * Aobs;
    var_part = sum(Aobs .* SA, 1)';
    var_part = max(var_part, 0);

    margin = mean_part - z_val * sqrt(var_part) - bobs;
    feasible_edge = find(margin >= -1e-6, 1, 'first');

    if isempty(feasible_edge)
        is_feasible = false;
        valid_pairs = {};
        return;
    end

    pair_count = pair_count + 1;
    valid_pairs{pair_count} = [obs_idx, feasible_edge];
end

valid_pairs = valid_pairs(1:pair_count);
end
