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
method_list = ["UN", "CC"];
risk_threshold = 0.5;
table_col = cellstr(method_list);
p_fn_list = [0];

Nmc = 100;

TT = length(p_fn_list);
FF = length(filter_list);
MM = length(method_list);
rmse_all = zeros(TT * FF * MM, 1);
variance_all = zeros(TT * FF * MM, 1);
compTime_all = zeros(TT * FF * MM, 1);
violProb_all = zeros(TT * FF * MM, 1);
violProb_all_max = zeros(TT * FF * MM, 1);

for tt = 1:length(p_fn_list)
    for ff = 1:length(filter_list)
        for mm = 1:length(method_list)
            p_fn = p_fn_list(tt);
            filter_mode = filter_list(ff);
            method_name = method_list(mm);

            localize_sim;

            rmse_list{ff, mm} = rmse_single;
            variance_list{ff, mm} = variance_single;
            comptime_list{ff, mm} = totalTime / (K + 1);
            violProb_list{ff, mm} = violProb_single;

            idx = (tt - 1) * FF * MM + (ff - 1) * MM + mm;
            rmse_all(idx) = rmse_total;
            variance_all(idx) = variance_total;
            compTime_all(idx) = calTime * 1000;
            violProb_all(idx) = vioProb_total;
            violProb_all_max(idx) = vioProb_max;
        end
    end
end

avgCompTime = zeros(FF, MM);
avgRMSE = zeros(FF, MM);
avgViolProb = zeros(FF, MM);
avgViolProbMax = zeros(FF, MM);

for ff = 1:FF
    for mm = 1:MM
        vals_comp = zeros(TT, 1);
        vals_rmse = zeros(TT, 1);
        vals_viol = zeros(TT, 1);
        vals_viol_max = zeros(TT, 1);

        for tt = 1:TT
            idx = (tt - 1) * FF * MM + (ff - 1) * MM + mm;
            vals_comp(tt) = compTime_all(idx);
            vals_rmse(tt) = rmse_all(idx);
            vals_viol(tt) = violProb_all(idx);
            vals_viol_max(tt) = violProb_all_max(idx);
        end
        avgCompTime(ff, mm) = mean(vals_comp);
        avgRMSE(ff, mm) = mean(vals_rmse);
        avgViolProb(ff, mm) = mean(vals_viol);
        avgViolProbMax(ff, mm) = mean(vals_viol_max);
    end
end

T_comp = array2table(avgCompTime, 'VariableNames', table_col, 'RowNames', cellstr(filter_list));
T_rmse = array2table(avgRMSE, 'VariableNames', table_col, 'RowNames', cellstr(filter_list));
T_viol = array2table(avgViolProb, 'VariableNames', table_col, 'RowNames', cellstr(filter_list));
T_viol_max = array2table(avgViolProbMax, 'VariableNames', table_col, 'RowNames', cellstr(filter_list));

disp('---- Average Computation Time Table ----');
disp(T_comp);

disp('---- Average RMSE Table ----');
disp(T_rmse);

disp('---- Average Violation Probability Table ----');
disp(T_viol);

disp('---- Max Violation Probability Table ----');
disp(T_viol_max);
