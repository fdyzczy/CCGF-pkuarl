clear all;
close all;
script_dir = fileparts(mfilename('fullpath'));
package_dir = fileparts(script_dir);
addpath(package_dir);
init_code_openacc(script_dir);
cd(script_dir);

mosek_data_file = fullfile(script_dir, "mosek_data_list2.mat");
if ~isfile(mosek_data_file)
    error('Missing precomputed MOSEK template file: %s', char(mosek_data_file));
end

S = load(mosek_data_file);
mosek_data_list = S.mosek_data_list;

filter_list = ["UKF", "GSF"];
num_landmarks = 16;
constr_list = ["CC1", "CC2", "CC3", "CC4"];
p_fn_list = [0];

Nmc = 5;

TT = length(p_fn_list);
FF = length(filter_list);
GG = length(constr_list);
rmse_all = zeros(TT * FF * GG, 1);
variance_all = zeros(TT * FF * GG, 1);
compTime_all = zeros(TT * FF * GG, 1);
violProb_all = zeros(TT * FF * GG, 1);
violProb_all_max = zeros(TT * FF * GG, 1);

for tt = 1:length(p_fn_list)
    for ff = 1:length(filter_list)
        for gg = 1:length(constr_list)
            p_fn = p_fn_list(tt);
            filter_mode = filter_list(ff);
            constr_mode = constr_list(gg);

            localize_sim;

            rmse_list{ff, gg} = rmse_single;
            variance_list{ff, gg} = variance_single;
            comptime_list{ff, gg} = totalTime / (K + 1);
            violProb_list{ff, gg} = violProb_single;

            idx = (tt - 1) * FF * GG + (ff - 1) * GG + gg;
            rmse_all(idx) = rmse_total;
            variance_all(idx) = variance_total;
            compTime_all(idx) = calTime * 1000;
            violProb_all(idx) = vioProb_total;
            violProb_all_max(idx) = vioProb_max;
        end
    end
end

avgCompTime = zeros(FF, GG);
avgRMSE = zeros(FF, GG);
avgViolProb = zeros(FF, GG);
avgViolProbMax = zeros(FF, GG);

for ff = 1:FF
    for gg = 1:GG
        vals_comp = zeros(TT, 1);
        vals_rmse = zeros(TT, 1);
        vals_viol = zeros(TT, 1);
        vals_viol_max = zeros(TT, 1);

        for tt = 1:TT
            idx = (tt - 1) * FF * GG + (ff - 1) * GG + gg;
            vals_comp(tt) = compTime_all(idx);
            vals_rmse(tt) = rmse_all(idx);
            vals_viol(tt) = violProb_all(idx);
            vals_viol_max(tt) = violProb_all_max(idx);
        end
        avgCompTime(ff, gg) = mean(vals_comp);
        avgRMSE(ff, gg) = mean(vals_rmse);
        avgViolProb(ff, gg) = mean(vals_viol);
        avgViolProbMax(ff, gg) = mean(vals_viol_max);
    end
end

T_comp = array2table(avgCompTime, 'VariableNames', cellstr(constr_list), 'RowNames', cellstr(filter_list));
T_rmse = array2table(avgRMSE, 'VariableNames', cellstr(constr_list), 'RowNames', cellstr(filter_list));
T_viol = array2table(avgViolProb, 'VariableNames', cellstr(constr_list), 'RowNames', cellstr(filter_list));
T_viol_max = array2table(avgViolProbMax, 'VariableNames', cellstr(constr_list), 'RowNames', cellstr(filter_list));

disp('---- Average Computation Time Table ----');
disp(T_comp);

disp('---- Average RMSE Table ----');
disp(T_rmse);

disp('---- Average Violation Probability Table ----');
disp(T_viol);

disp('---- Max Violation Probability Table ----');
disp(T_viol_max);
