% Road-tracking example with unconstrained and CC-constrained filters.

rng(1);

delta = risk_threshold;
is_constrained = string(method_name) == "CC";
if ~(is_constrained || string(method_name) == "UN")
    error('Unsupported method: %s', string(method_name));
end

delta_set = delta * ones(1, length(constr)) / length(obstacles);
Nsample_viol = 100000;

dt = 0.5;
K = 50;
ndim = 4;

robot_traj = target_traj_list';
robot_traj(3, :) = wrapToPi(robot_traj(3, :) + pi);
robot_traj = robot_traj(:, 1:K);
robot_vel = diff([robot_traj(:, 1), robot_traj], 1, 2) / dt;
robot_vel_norm = sqrt(robot_vel(1, :).^2 + robot_vel(2, :).^2);
robot_pos = [robot_traj; robot_vel_norm];

R = 0.3^2 * eye(2);
Q = diag([0.02^2, 0.02^2, 0.25^2, 0.5^2]);

mixNum = 3;
x0_gsf = repmat(robot_pos(:, 1), 1, mixNum) + [[3; 0; pi / 2; 0], [0; 3; 0; 0.2], [-3; 0; -pi / 2; -0.2]];
P0_gsf = repmat(blkdiag(1^2 * eye(2), 0.3^2 * eye(1), 0.2^2 * eye(1)), 1, 1, mixNum);
w0_gsf = ones(1, mixNum) / mixNum;

x0 = robot_pos(:, 1) + [-3; 0; -pi / 2; 0];
P0 = blkdiag(4^2 * eye(2), 1^2 * eye(1), 0.3^2 * eye(1));

error_x1 = zeros(Nmc, K + 1);
error_x2 = zeros(Nmc, K + 1);
error_x3 = zeros(Nmc, K + 1);
error_x4 = zeros(Nmc, K + 1);
variance = zeros(Nmc, K + 1);
violProb = zeros(Nmc, K + 1);
totalTime = zeros(Nmc, 1);

stateTransitionFcn = @bicycle_model;
measurementFcn = @(x) [x(1); x(2)];

for mc = 1:Nmc
    x_true = zeros(ndim, K + 1);
    x_true(:, 1:K) = robot_pos;
    x_true(:, K + 1) = x_true(:, K);
    meas = meas_list(1:K + 1, :)';

    switch string(filter_mode)
        case "GSF"
            subfilters = cell(1, mixNum);
            for i = 1:mixNum
                subfilters{i} = trackingEKF('StateTransitionFcn', stateTransitionFcn, ...
                    'MeasurementFcn', measurementFcn, ...
                    'State', x0_gsf(:, i), ...
                    'StateCovariance', P0_gsf(:, :, i), ...
                    'HasAdditiveProcessNoise', false, ...
                    'ProcessNoise', Q, ...
                    'MeasurementNoise', R);
            end
            filter = trackingGSF("TrackingFilters", subfilters, ...
                "ModelProbabilities", w0_gsf, "MeasurementNoise", R);

        case "EKF"
            filter = trackingEKF('StateTransitionFcn', stateTransitionFcn, ...
                'MeasurementFcn', measurementFcn, ...
                'State', x0, ...
                'StateCovariance', P0, ...
                'HasAdditiveProcessNoise', false, ...
                'ProcessNoise', Q, ...
                'MeasurementNoise', R);

        otherwise
            error('Unsupported filter mode in the CC-only release: %s', string(filter_mode));
    end

    x_est = zeros(ndim, K + 1);
    x_cov = zeros(ndim, ndim, K + 1);
    x_est(:, 1) = filter.State;

    for k = 1:K + 1
        tic;
        z_k = meas(:, k);
        correct(filter, z_k');

        if string(filter_mode) == "GSF"
            mu0 = zeros(ndim, mixNum);
            Sigma0 = zeros(ndim, ndim, mixNum);
            for i = 1:length(filter.ModelProbabilities)
                mu0(:, i) = filter.TrackingFilters{i}.State;
                Sigma0(:, :, i) = filter.TrackingFilters{i}.StateCovariance;
            end
            w0 = filter.ModelProbabilities;
            mu = mu0;
            Sigma = Sigma0;

            if is_constrained
                for i = 1:mixNum
                    [mu(:, i), Sigma(:, :, i)] = BranchAndBound_pruned( ...
                        mu(:, i), Sigma(:, :, i), obstacles, delta_set, mosek_data_list);
                end
                weight = weight_alloc(mu0, mu, Sigma0, Sigma, w0);

                for i = 1:length(filter.ModelProbabilities)
                    filter.TrackingFilters{i}.State = mu(:, i);
                    filter.TrackingFilters{i}.StateCovariance = Sigma(:, :, i);
                end
                filter.ModelProbabilities = weight;
            end
        elseif is_constrained
            [filter.State, filter.StateCovariance] = BranchAndBound_pruned( ...
                filter.State, filter.StateCovariance, obstacles, delta_set, mosek_data_list);
        end

        x_est(:, k) = filter.State;
        x_est(3, k) = wrapToPi(x_est(3, k));
        x_cov(:, :, k) = filter.StateCovariance;
        totalTime(mc) = totalTime(mc) + toc;
        variance(mc, k) = trace(filter.StateCovariance);

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
            predict(filter);
        end
        totalTime(mc) = totalTime(mc) + toc;
    end

    error_x1(mc, :) = x_est(1, :) - x_true(1, :);
    error_x2(mc, :) = x_est(2, :) - x_true(2, :);
    error_x3(mc, :) = wrapToPi(x_est(3, :) - x_true(3, :));
    error_x4(mc, :) = x_est(4, :) - x_true(4, :);

    x_est_1(mc, :) = x_est(1, :);
    x_est_2(mc, :) = x_est(2, :);

end

error_vec = error_x1.^2 + error_x2.^2 + error_x3.^2 + error_x4.^2;
error_x_vec = error_x1.^2 + error_x2.^2;
error_theta_vec = error_x3.^2;
error_v_vec = error_x4.^2;

rmse_avg_step = sqrt(sum(error_vec, 1) / Nmc);
rmse_x_single = sqrt(sum(error_x_vec, 2) / (K + 1));
rmse_theta_single = sqrt(sum(error_theta_vec, 2) / (K + 1));
rmse_v_single = sqrt(sum(error_v_vec, 2) / (K + 1));
rmse_single = sqrt(sum(error_vec, 2) / (K + 1));

rmse_x_total = mean(rmse_x_single);
rmse_theta_total = mean(rmse_theta_single);
rmse_v_total = mean(rmse_v_single);
rmse_total = mean(rmse_single);

calTime = mean(totalTime) / (K + 1);
variance_single = mean(variance, 2);
variance_total = mean(variance_single);
violProb_single = max(violProb, [], 2);
vioProb_total = mean(violProb_single);
vioProb_max = max(violProb_single);

disp('-------------------Simulation Ends-----------------------')
fprintf('Filter: %s-%s\n', method_name, string(filter_mode));
fprintf('Risk Threshold: %g\n', delta);
fprintf('RMSE: %.4f\n', rmse_total);
fprintf('RMSE_position: %.4f\n', rmse_x_total);
fprintf('RMSE_angle: %.4f\n', rmse_theta_total);
fprintf('RMSE_velocity: %.4f\n', rmse_v_total);
fprintf('Computational Time: %.4f ms\n', calTime * 1000);
fprintf('Variance: %f\n', variance_total);
fprintf('Violation Probability: %f\n', vioProb_total);
