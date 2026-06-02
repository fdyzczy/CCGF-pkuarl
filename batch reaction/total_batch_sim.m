clear all;
close all;
script_dir = fileparts(mfilename('fullpath'));
package_dir = fileparts(script_dir);
addpath(package_dir);
init_code_openacc(script_dir);
cd(script_dir);

filter_list = ["EKF", "GSF"];
constr_list = ["CC"];

Nmc = 10;
rand_seed = 3;

for ff = 1:length(filter_list)
    for gg = 1:length(constr_list)
        filter_mode = filter_list(ff);
        constr_mode = lower(constr_list(gg));

        batch_sim;

        rmse_list{ff, gg} = rmse_single;
        variance_list{ff, gg} = variance_single;
        comptime_list{ff, gg} = totalTime / (K + 1);
        violProb_list{ff, gg} = violProb_single;

        rmse_all{ff, gg} = rmse_total;
        variance_all{ff, gg} = variance_total;
        compTime_all{ff, gg} = calTime * 1000;
        violProb_all{ff, gg} = vioProb_total;
        violProb_all_max{ff, gg} = vioProb_max;
    end
end

T_rmse = array2table(cell2mat(rmse_all), ...
    'VariableNames', cellstr(constr_list), ...
    'RowNames', cellstr(filter_list));
disp('RMSE Table:');
disp(T_rmse);

T_time = array2table(cell2mat(compTime_all), ...
    'VariableNames', cellstr(constr_list), ...
    'RowNames', cellstr(filter_list));
disp('Calculation Time Table:');
disp(T_time);

T_viol = array2table(cell2mat(violProb_all), ...
    'VariableNames', cellstr(constr_list), ...
    'RowNames', cellstr(filter_list));
disp('Violation Probability Table:');
disp(T_viol);

T_viol_max = array2table(cell2mat(violProb_all_max), ...
    'VariableNames', cellstr(constr_list), ...
    'RowNames', cellstr(filter_list));
disp('Max Violation Probability Table:');
disp(T_viol_max);
