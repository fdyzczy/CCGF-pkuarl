% Robot-localization example with unconstrained and CC-constrained filters.

rand_seed = 32;
rng(rand_seed);

delta = risk_threshold;
is_constrained = string(method_name) == "CC";
if ~(is_constrained || string(method_name) == "UN")
    error('Unsupported method: %s', string(method_name));
end

FOV = 6;
scene_size = 60 * sqrt(num_landmarks / 16);
bias_std = 0;
lm_x = linspace(8, scene_size - 5, sqrt(num_landmarks));
lm_y = linspace(10, scene_size - 5, sqrt(num_landmarks));
[LM_X, LM_Y] = meshgrid(lm_x, lm_y);
landmarks_base = [LM_X(:), LM_Y(:)];
landmarks = landmarks_base + bias_std * randn(num_landmarks, 2);

dt = 1.0;
start_pt = [1; 1];
end_pt = [55; 55] * scene_size / 60;
K = 100;
x_path = linspace(start_pt(1), end_pt(1), K + 1);
L = end_pt(1) - start_pt(1);
amplitude = (end_pt(2) - start_pt(2)) / 2;
n_periods = 2;
y_path = start_pt(2) + amplitude * (sin(2 * pi * n_periods * (x_path - start_pt(1)) / L - pi / 2) + 1);
robot_traj = [x_path; y_path];
robot_traj = robot_traj(:, 1:K);
robot_vel = diff([robot_traj(:, 1), robot_traj], 1, 2) / dt;

sigma_r = 1;
sigma_b = deg2rad(10);
R = diag([sigma_r^2, sigma_b^2]);
Q = 0.5 * eye(2);

mixNum = 4;
x0_gsf = repmat(robot_traj(:, 1), 1, mixNum) + [[4; 4], [4; 10], [10; 4], [10; 10]];
P0_gsf = repmat(7^2 * eye(2), 1, 1, mixNum);
w0_gsf = ones(1, mixNum) / mixNum;

x0 = robot_traj(:, 1) + [6; 6];
P0 = 10^2 * eye(2);

z_all = cell(K + 1, 1);
lm_id_all = cell(K + 1, 1);

error_x1 = zeros(Nmc, K + 1);
error_x2 = zeros(Nmc, K + 1);
variance = zeros(Nmc, K + 1);
violProb = zeros(Nmc, K + 1);
totalTime = zeros(Nmc, 1);
ExpandNodeNum = zeros(Nmc, K + 1);

stateTransitionFcn = @(xk, uk) xk + uk * dt;
measurementFcn = @range_bearing_measurement;

for mc = 1:Nmc
    x0_true = robot_traj(:, 1);
    while true
        x_true = zeros(2, K + 1);
        x_true(:, 1) = x0_true;
        for k = 1:K
            x_true(:, k + 1) = stateTransitionFcn(x_true(:, k), robot_vel(:, k)) + mvnrnd([0; 0], Q)';
        end
        x_true(:, K + 1) = x_true(:, K);

        if all(x_true(1, :) >= 0 & x_true(1, :) <= scene_size & ...
               x_true(2, :) >= 0 & x_true(2, :) <= scene_size)
            break;
        end
    end

    for k = 1:K + 1
        meas = [];
        id = [];
        pos = x_true(:, k);
        for i = 1:num_landmarks
            vec = landmarks(i, :)' - pos;
            dist = norm(vec);
            if dist <= FOV && rand > p_fn
                bearing = atan2(vec(2), vec(1));
                meas = [meas; dist + sigma_r * randn(1), wrapToPi(bearing + sigma_b * randn(1))];
                id = [id; i];
            end
        end
        z_all{k} = meas;
        lm_id_all{k} = id;
    end

    switch string(filter_mode)
        case "GSF"
            subfilters = cell(1, mixNum);
            for i = 1:mixNum
                subfilters{i} = trackingUKF('StateTransitionFcn', stateTransitionFcn, ...
                    'MeasurementFcn', measurementFcn, ...
                    'State', x0_gsf(:, i), ...
                    'StateCovariance', P0_gsf(:, :, i), ...
                    'ProcessNoise', Q, ...
                    'MeasurementNoise', R, ...
                    'HasMeasurementWrapping', true);
            end
            filter = trackingGSF("TrackingFilters", subfilters, ...
                "ModelProbabilities", w0_gsf, "MeasurementNoise", R);

        case "UKF"
            filter = trackingUKF('StateTransitionFcn', stateTransitionFcn, ...
                'MeasurementFcn', measurementFcn, ...
                'State', x0, ...
                'StateCovariance', P0, ...
                'ProcessNoise', Q, ...
                'MeasurementNoise', R, ...
                'HasMeasurementWrapping', true);

        otherwise
            error('Unsupported filter mode in the CC-only release: %s', string(filter_mode));
    end

    x_est = zeros(2, K + 1);
    x_cov = zeros(2, 2, K + 1);
    x_est(:, 1) = filter.State;

    for k = 1:K + 1
        tic;
        z_k = z_all{k};
        id_k = lm_id_all{k};
        constr = {};

        for i = 1:size(landmarks, 1)
            if any(i == id_k)
                continue;
            end
            obs.type = "cur";
            obs.center = landmarks(i, :);
            obs.a = FOV;
            obs.b = obs.a;
            obs.angle = 0;
            constr{end + 1} = obs;
        end

        for j = 1:length(id_k)
            correct(filter, z_k(j, :)', landmarks(id_k(j), :));
        end

        obstacles = parse_constr_to_Obs(constr);
        if isempty(obstacles)
            x_est(:, k) = filter.State;
            x_cov(:, :, k) = filter.StateCovariance;
            totalTime(mc) = totalTime(mc) + toc;
            variance(mc, k) = trace(filter.StateCovariance);
            tic;
            if k < K + 1
                vel = robot_vel(:, k);
                predict(filter, vel);
            end
            totalTime(mc) = totalTime(mc) + toc;
            continue;
        end
        delta_set = delta * ones(1, length(obstacles)) / length(obstacles);

        if string(filter_mode) == "GSF"
            mu0 = zeros(2, mixNum);
            Sigma0 = zeros(2, 2, mixNum);
            for i = 1:length(filter.ModelProbabilities)
                mu0(:, i) = filter.TrackingFilters{i}.State;
                Sigma0(:, :, i) = filter.TrackingFilters{i}.StateCovariance;
            end
            w0 = filter.ModelProbabilities;
            mu = mu0;
            Sigma = Sigma0;
            ExpandNodeNum(mc, k) = 0;

            if is_constrained
                info = cell(1, mixNum);
                for i = 1:mixNum
                    [mu(:, i), Sigma(:, :, i), ~, info{i}] = BranchAndBound_pruned( ...
                        mu(:, i), Sigma(:, :, i), obstacles, delta_set, mosek_data_list);
                    ExpandNodeNum(mc, k) = ExpandNodeNum(mc, k) + info{i}.iterations;
                end

                ExpandNodeNum(mc, k) = ExpandNodeNum(mc, k) / mixNum;
                weight = weight_alloc(mu0, mu, Sigma0, Sigma, w0);

                for i = 1:length(filter.ModelProbabilities)
                    filter.TrackingFilters{i}.State = mu(:, i);
                    filter.TrackingFilters{i}.StateCovariance = Sigma(:, :, i);
                end
                filter.ModelProbabilities = weight;
            end
        elseif is_constrained
            [filter.State, filter.StateCovariance, ~, info] = BranchAndBound_pruned( ...
                filter.State, filter.StateCovariance, obstacles, delta_set, mosek_data_list);
            ExpandNodeNum(mc, k) = info.iterations;
        end

        x_est(:, k) = filter.State;
        x_cov(:, :, k) = filter.StateCovariance;
        totalTime(mc) = totalTime(mc) + toc;
        variance(mc, k) = trace(filter.StateCovariance);

        Nsample_viol = 10000;
        if string(filter_mode) == "GSF"
            viol_num_indiv = 0;
            Ns = zeros(1, length(filter.ModelProbabilities));
            for i = 1:length(filter.ModelProbabilities)
                mu(:, i) = filter.TrackingFilters{i}.State;
                Sigma(:, :, i) = filter.TrackingFilters{i}.StateCovariance;
            end
            w = filter.ModelProbabilities;

            for m = 1:length(w)
                Ns(m) = max(0, round(w(m) * Nsample_viol));
                if Ns(m) == 0
                    continue;
                end
                samples = mvnrnd(mu(:, m)', Sigma(:, :, m), Ns(m));
                viol_vec = zeros(Ns(m), length(constr));
                for c = 1:length(constr)
                    viol_vec(:, c) = checkConstr(samples', constr{c});
                end
                viol_num_indiv = viol_num_indiv + sum(sum(viol_vec') > 0);
            end
            cumu_mu{k} = mu;
            cumu_Sigma{k} = Sigma;
            cumu_w{k} = w;
            violProb(mc, k) = viol_num_indiv / sum(Ns);
        else
            samples = mvnrnd(filter.State', filter.StateCovariance, Nsample_viol);
            viol_vec = zeros(Nsample_viol, length(constr));
            for c = 1:length(constr)
                viol_vec(:, c) = checkConstr(samples', constr{c});
            end
            if ~isempty(constr)
                violProb(mc, k) = mean(sum(viol_vec') > 0);
            end
        end

        tic;
        if k < K + 1
            vel = robot_vel(:, k);
            predict(filter, vel);
        end
        totalTime(mc) = totalTime(mc) + toc;
    end

    error_x1(mc, :) = x_est(1, :) - x_true(1, :);
    error_x2(mc, :) = x_est(2, :) - x_true(2, :);
    x_est_1(mc, :) = x_est(1, :);
    x_est_2(mc, :) = x_est(2, :);

end

error_vec = error_x1.^2 + error_x2.^2;
rmse_avg_step = sqrt(sum(error_vec, 1) / Nmc);
rmse_single = sqrt(sum(error_vec, 2) / (K + 1));
idx = 1:Nmc;

rmse_total = mean(rmse_single(idx));
rmse_std = std(rmse_single(idx), 0);

variance_single = mean(variance, 2);
variance_total = mean(variance_single(idx));
variance_std = std(variance_single(idx), 0);

calTime_samples = totalTime / (K + 1);
calTime = mean(calTime_samples);
calTime_std = std(calTime_samples, 0);

violProb_single = max(violProb, [], 2);
vioProb_total = mean(violProb_single(idx));
vioProb_max = max(violProb_single(idx));

ExpandNodeNum_total = mean(mean(ExpandNodeNum, 2));
ExpandNodeNum_std = std(mean(ExpandNodeNum, 2));

disp('-------------------Simulation Ends-----------------------')
fprintf('landmark number: %f\n', num_landmarks);
fprintf('Filter: %s-%s\n', method_name, string(filter_mode));
fprintf('Risk Threshold: %g\n', delta);
fprintf('RMSE: %.4f\n', rmse_total);
fprintf('Computational Time: %.4f ms\n', calTime * 1000);
fprintf('Variance: %f\n', variance_total);
fprintf('Violation Probability: %f\n', vioProb_total);
