clear all;
close all;
script_dir = fileparts(mfilename('fullpath'));
package_dir = fileparts(script_dir);
addpath(package_dir);
init_code_openacc(script_dir);
cd(script_dir);

mosek_data_file = fullfile(script_dir, "mosek_data_list.mat");
if ~isfile(mosek_data_file)
    error('Missing precomputed MOSEK template file: %s', char(mosek_data_file));
end

S = load(mosek_data_file);
mosek_data_list = S.mosek_data_list;
constr_data = load(fullfile(script_dir, "constr_data.mat"));
avgPoints = constr_data.old_avgPoints;
constr = constr_data.constr;
constr{end + 1}.type = 'lin';
constr{end}.a = [0; 0; 0; 1];
constr{end}.b = 0;

obstacles = parse_constr_to_Obs(constr, 4);
traj_file_list = "exp_data" + string(1:5) + ".mat";
TT = length(traj_file_list);

filter_list = ["GSF", "EKF"];
constr_list = ["CC1"];
FF = length(filter_list);
GG = length(constr_list);

Nmc = 2;

for tt = 1:length(traj_file_list)
    fprintf('*****Trajectory index: %d*****\n', tt);
    for ff = 1:length(filter_list)
        for gg = 1:length(constr_list)
            filter_mode = filter_list(ff);
            constr_mode = constr_list(gg);
            traj_data = load(fullfile(script_dir, traj_file_list(tt)));
            target_traj_list = traj_data.pos2d_interval;
            meas_list = traj_data.meas_interval;

            trackingSim;

            rmse_list{tt, ff, gg} = rmse_single;
            variance_list{tt, ff, gg} = variance_single;
            comptime_list{tt, ff, gg} = totalTime / (K + 1);
            violProb_list{tt, ff, gg} = violProb_single;

            rmse_all{tt, ff, gg} = rmse_total;
            rmse_x_all{tt, ff, gg} = rmse_x_total;
            rmse_theta_all{tt, ff, gg} = rmse_theta_total;
            rmse_v_all{tt, ff, gg} = rmse_v_total;
            variance_all{tt, ff, gg} = variance_total;
            compTime_all{tt, ff, gg} = calTime * 1000;
            violProb_all{tt, ff, gg} = vioProb_total;
            violProb_all_max{tt, ff, gg} = vioProb_max;
        end
    end
end

show_results(compTime_all, rmse_all, violProb_all, violProb_all_max, ...
    TT, FF, GG, constr_list, filter_list, ...
    rmse_x_all, rmse_theta_all, rmse_v_all);

function show_results(compTime_all, rmse_all, violProb_all, violProb_all_max, ...
    TT, FF, GG, constr_list, filter_list, ...
    rmse_x_all, rmse_theta_all, rmse_v_all, exclude_tt)

    if nargin < 14
        exclude_tt = [];
    end

    tt_idx = 1:TT;
    if ~isempty(exclude_tt)
        tt_idx(exclude_tt) = [];
    end

    compTime_mat = cell2mat(reshape(compTime_all, TT, FF, GG));
    rmse_mat = cell2mat(reshape(rmse_all, TT, FF, GG));
    violProb_mat = cell2mat(reshape(violProb_all, TT, FF, GG));
    violProbMax_mat = cell2mat(reshape(violProb_all_max, TT, FF, GG));
    rmse_x_mat = cell2mat(reshape(rmse_x_all, TT, FF, GG));
    rmse_theta_mat = cell2mat(reshape(rmse_theta_all, TT, FF, GG));
    rmse_v_mat = cell2mat(reshape(rmse_v_all, TT, FF, GG));

    compTime_mat = compTime_mat(tt_idx, :, :);
    rmse_mat = rmse_mat(tt_idx, :, :);
    violProb_mat = violProb_mat(tt_idx, :, :);
    violProbMax_mat = violProbMax_mat(tt_idx, :, :);
    rmse_x_mat = rmse_x_mat(tt_idx, :, :);
    rmse_theta_mat = rmse_theta_mat(tt_idx, :, :);
    rmse_v_mat = rmse_v_mat(tt_idx, :, :);

    avgCompTime = squeeze(mean(compTime_mat, 1));
    avgRMSE = squeeze(mean(rmse_mat, 1));
    avgViolProb = squeeze(mean(violProb_mat, 1));
    avgViolProbMax = squeeze(max(violProbMax_mat, [], 1));
    avgRMSE_x = squeeze(mean(rmse_x_mat, 1));
    avgRMSE_theta = squeeze(mean(rmse_theta_mat, 1));
    avgRMSE_v = squeeze(mean(rmse_v_mat, 1));

    if GG == 1
        avgCompTime = avgCompTime(:);
        avgRMSE = avgRMSE(:);
        avgViolProb = avgViolProb(:);
        avgViolProbMax = avgViolProbMax(:);
        avgRMSE_x = avgRMSE_x(:);
        avgRMSE_theta = avgRMSE_theta(:);
        avgRMSE_v = avgRMSE_v(:);
    end

    T_comp = array2table(avgCompTime, 'VariableNames', cellstr(constr_list), 'RowNames', cellstr(filter_list));
    T_rmse = array2table(avgRMSE, 'VariableNames', cellstr(constr_list), 'RowNames', cellstr(filter_list));
    T_viol = array2table(avgViolProb, 'VariableNames', cellstr(constr_list), 'RowNames', cellstr(filter_list));
    T_viol_max = array2table(avgViolProbMax, 'VariableNames', cellstr(constr_list), 'RowNames', cellstr(filter_list));
    T_rmse_x = array2table(avgRMSE_x, 'VariableNames', cellstr(constr_list), 'RowNames', cellstr(filter_list));
    T_rmse_theta = array2table(avgRMSE_theta, 'VariableNames', cellstr(constr_list), 'RowNames', cellstr(filter_list));
    T_rmse_v = array2table(avgRMSE_v, 'VariableNames', cellstr(constr_list), 'RowNames', cellstr(filter_list));

    disp('---- Average Computation Time Table ----');
    disp(T_comp);

    disp('---- Average RMSE Table ----');
    disp(T_rmse);

    disp('---- Average Violation Probability Table ----');
    disp(T_viol);

    disp('---- Max Violation Probability Table ----');
    disp(T_viol_max);

    disp('---- Average RMSE (Position) Table ----');
    disp(T_rmse_x);

    disp('---- Average RMSE (Angle) Table ----');
    disp(T_rmse_theta);

    disp('---- Average RMSE (Velocity) Table ----');
    disp(T_rmse_v);
end
