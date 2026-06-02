% Batch-reaction example using only the proposed CC truncation method.

rng(rand_seed);

constr_mode = lower(string(constr_mode));
if ~startsWith(constr_mode, "cc")
    error('code_openacc keeps only the proposed CC method. Unsupported method: %s', constr_mode);
end
constr_mode = "cc";

dt = 0.1;
kr = 0.16;
K = 20;
Q = 0.002^2 * eye(2);
R = 0.5^2;

delta = 1e-4;
Nsample_viol = 100000;

f = @(x) [ ...
    x(1) / (1 + 2 * kr * dt * x(1));
    x(2) + (kr * dt * x(1)) / (1 + 2 * kr * dt * x(1))
    ];
h = @(x) [1 1] * x;

mixNum = 3;
x0_gsf = [[0.1; 0.1], [6; 0.1], [0.1; 6]];
P0_gsf = repmat(4^2 * eye(2), 1, 1, mixNum);
w0_gsf = ones(1, mixNum) / mixNum;

x0_filter = [0.1; 6];
P0_filter = 6^2 * eye(2);
x0_true = [2; 2];

constr{1}.type = "lin";
constr{1}.a = [1; 0];
constr{1}.b = 0;
constr{2}.type = "lin";
constr{2}.a = [0; 1];
constr{2}.b = 0;
delta_set = delta * ones(1, length(constr)) / length(constr);
obstacles = parse_constr_to_Obs(constr);

error_x1 = zeros(Nmc, K + 1);
error_x2 = zeros(Nmc, K + 1);
variance = zeros(Nmc, K + 1);
violProb = zeros(Nmc, K + 1);
totalTime = zeros(Nmc, 1);

for mc = 1:Nmc
    x_true = zeros(2, K + 1);
    z = zeros(1, K + 1);
    x_true(:, 1) = x0_true;

    for k = 1:K
        x_true(:, k + 1) = f(x_true(:, k)) + mvnrnd([0; 0], Q)';
    end

    for k = 1:K + 1
        z(k) = h(x_true(:, k)) + normrnd(0, sqrt(R));
    end

    switch string(filter_mode)
        case "GSF"
            subfilters = cell(1, mixNum);
            for i = 1:mixNum
                subfilters{i} = trackingEKF( ...
                    StateTransitionFcn = f, ...
                    MeasurementFcn = h, ...
                    ProcessNoise = Q, ...
                    MeasurementNoise = R, ...
                    State = x0_gsf(:, i), ...
                    StateCovariance = P0_gsf(:, :, i));
            end
            filter = trackingGSF("TrackingFilters", subfilters, ...
                "ModelProbabilities", w0_gsf, "MeasurementNoise", R);

        case "UKF"
            filter = trackingUKF( ...
                StateTransitionFcn = f, ...
                MeasurementFcn = h, ...
                ProcessNoise = Q, ...
                MeasurementNoise = R, ...
                State = x0_filter, ...
                StateCovariance = P0_filter, ...
                kappa = 1);

        case "EKF"
            filter = trackingEKF( ...
                StateTransitionFcn = f, ...
                MeasurementFcn = h, ...
                ProcessNoise = Q, ...
                MeasurementNoise = R, ...
                State = x0_filter, ...
                StateCovariance = P0_filter);

        otherwise
            error('Unsupported filter mode in the CC-only release: %s', string(filter_mode));
    end

    x_est = zeros(2, K + 1);
    x_cov = zeros(2, 2, K + 1);

    for k = 1:K + 1
        tic;
        correct(filter, z(k));

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
            for i = 1:mixNum
                [mu(:, i), Sigma(:, :, i)] = BranchAndBound_pruned( ...
                    mu(:, i), Sigma(:, :, i), obstacles, delta_set);
            end
            weight = weight_alloc(mu0, mu, Sigma0, Sigma, w0);

            for i = 1:length(filter.ModelProbabilities)
                filter.TrackingFilters{i}.State = mu(:, i);
                filter.TrackingFilters{i}.StateCovariance = Sigma(:, :, i);
            end
            filter.ModelProbabilities = weight;
        else
            [filter.State, filter.StateCovariance] = BranchAndBound_pruned( ...
                filter.State, filter.StateCovariance, obstacles, delta_set);
        end

        x_est(:, k) = filter.State;
        x_cov(:, :, k) = filter.StateCovariance;
        totalTime(mc) = totalTime(mc) + toc;
        variance(mc, k) = trace(filter.StateCovariance);

        constr_val = @(x) max([-x(1, :); -x(2, :)], [], 1);
        if string(filter_mode) == "GSF"
            viol_num = 0;
            total_Ns = 0;
            for i = 1:length(filter.ModelProbabilities)
                Ns = max(0, round(filter.ModelProbabilities(i) * Nsample_viol));
                total_Ns = total_Ns + Ns;
                if Ns == 0
                    continue;
                end
                samples = mvnrnd(mu(:, i)', Sigma(:, :, i), Ns);
                viol_num = viol_num + sum(constr_val(samples') > 0);
            end
            violProb(mc, k) = viol_num / total_Ns;
        else
            samples = mvnrnd(filter.State', filter.StateCovariance, Nsample_viol);
            violProb(mc, k) = mean(constr_val(samples') > 0);
        end

        tic;
        if k < K + 1
            predict(filter);
        end
        totalTime(mc) = totalTime(mc) + toc;
    end

    error_x1(mc, :) = x_est(1, :) - x_true(1, :);
    error_x2(mc, :) = x_est(2, :) - x_true(2, :);
    x_est_1(mc, :) = x_est(1, :);
    x_est_2(mc, :) = x_est(2, :);
end

error_vec = error_x1.^2 + error_x2.^2;

calTime_samples = totalTime / (K + 1);
calTime = mean(calTime_samples);
calTime_std = std(calTime_samples, 0);

variance_single = mean(variance, 2);
variance_total = mean(variance_single);
variance_std = std(variance, 0);

violProb_single = max(violProb, [], 2);
vioProb_total = mean(violProb_single);
vioProb_max = max(violProb_single);
vioProb_std = std(violProb_single);

rmse_single = sqrt(sum(error_vec, 2) / (K + 1));
rmse_total = mean(rmse_single);
rmse_std = std(rmse_single);

disp('-------------------Simulation Ends-----------------------')
fprintf('Filter: %s-%s\n', constr_mode, string(filter_mode));
fprintf('Risk Threshold: %g\n', delta);
fprintf('RMSE: %.4f\n', rmse_total);
fprintf('Computational Time: %.4f ms\n', calTime * 1000);
fprintf('Variance: %f\n', variance_total);
fprintf('Violation Probability: %f\n', vioProb_total);
