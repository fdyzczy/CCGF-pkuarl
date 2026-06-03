clear all;
close all;

script_dir = fileparts(mfilename('fullpath'));
package_dir = fileparts(script_dir);
addpath(package_dir);
init_code_openacc(script_dir);
cd(script_dir);

filter_list = ["EKF", "GSF"];
method_list = ["UN", "CC"];
risk_threshold = 1e-2;
table_col = cellstr(method_list);

Nmc = 10;
rand_seed = 3;

for ff = 1:length(filter_list)
    for mm = 1:length(method_list)
        filter_mode = filter_list(ff);
        method_name = method_list(mm);

        batch_sim;

        rmse_list{ff, mm} = rmse_single;
        variance_list{ff, mm} = variance_single;
        comptime_list{ff, mm} = totalTime / (K + 1);
        violProb_list{ff, mm} = violProb_single;

        rmse_all{ff, mm} = rmse_total;
        variance_all{ff, mm} = variance_total;
        compTime_all{ff, mm} = calTime * 1000;
        violProb_all{ff, mm} = vioProb_total;
        violProb_all_max{ff, mm} = vioProb_max;
    end
end

T_rmse = array2table(cell2mat(rmse_all), ...
    'VariableNames', table_col, ...
    'RowNames', cellstr(filter_list));
disp('RMSE Table:');
disp(T_rmse);

T_time = array2table(cell2mat(compTime_all), ...
    'VariableNames', table_col, ...
    'RowNames', cellstr(filter_list));
disp('Calculation Time Table:');
disp(T_time);

T_viol = array2table(cell2mat(violProb_all), ...
    'VariableNames', table_col, ...
    'RowNames', cellstr(filter_list));
disp('Violation Probability Table:');
disp(T_viol);

T_viol_max = array2table(cell2mat(violProb_all_max), ...
    'VariableNames', table_col, ...
    'RowNames', cellstr(filter_list));
disp('Max Violation Probability Table:');
disp(T_viol_max);
