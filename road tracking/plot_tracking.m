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
filter_list = ["EKF"];
constr_list = ["CC1"];
constr_legend = ["CC (0.05%)"];

Nmc = 1;

for tt = 1:1
    fprintf('*****Trajectory index: %d*****\n', tt);
    for ff = 1:length(filter_list)
        for gg = 1:length(constr_list)
            filter_mode = filter_list(ff);
            constr_mode = constr_list(gg);
            traj_data = load(fullfile(script_dir, traj_file_list(tt)));
            target_traj_list = traj_data.pos2d_interval;
            meas_list = traj_data.meas_interval;

            trackingSim;
            x_est_store{ff, gg} = x_est;
        end
    end
end

close all;
true_color = [0, 0, 0.9];
meas_color = [220, 109, 087] / 255;
cam_color = [220, 109, 087] / 255;

figure; clf; set(gcf, 'Position', [100 100 1000 600]);
hold on;
grid on;

roadhdl = plotRoad(avgPoints);
hold on; axis equal

truehdl = plot(x_true(1, 1:end), x_true(2, 1:end), '-', ...
    'color', true_color, 'LineWidth', 2, ...
    'Marker', '.', 'MarkerSize', 22, 'DisplayName', 'True Trajectory');

meashdl = scatter(meas(1, :), meas(2, :), 120, ...
    'MarkerEdgeColor', meas_color, 'marker', '+', 'LineWidth', 2);

camera_pos = [-2.933980228260293, -0.786720232065491, 0.606936347084067];
x_cam = camera_pos(1);
y_cam = camera_pos(2);
theta_cam = camera_pos(3);
arrow_length = 0.3;
dx = arrow_length * cos(theta_cam);
dy = arrow_length * sin(theta_cam);
camhdl = quiver(x_cam, y_cam, dx, dy, 0, ...
    'color', cam_color, 'LineWidth', 2, 'MaxHeadSize', 3, ...
    'Marker', 'square', 'MarkerSize', 12);

xlim([-3.2, 3.5]);
ylim([-1, 1.9]);

xlabel('X position(m)', ...
    'FontSize', 20, ...
    'HorizontalAlignment', 'right', ...
    'Units', 'normalized', ...
    'Position', [0.99, 0.12, 0], 'Interpreter', 'latex');

ylabel('Y position(m)', ...
    'FontSize', 20, ...
    'HorizontalAlignment', 'right', ...
    'Units', 'normalized', ...
    'Position', [0.22, 0.88, 0], 'Interpreter', 'latex', ...
    'Rotation', 0);

set(gca, 'FontSize', 18, 'FontName', 'Times New Roman', ...
    'XTickLabelMode', 'auto', 'YTickLabelMode', 'auto');
legend([roadhdl, truehdl, camhdl, meashdl], ...
    {'Road boundary', 'Target traj.\quad', 'Camera pose', 'Measurements'}, ...
    'FontSize', 18, 'Location', 'southoutside', 'Box', 'off', ...
    'Orientation', 'horizontal', 'Interpreter', 'latex');
box on;
drawnow;
exportgraphics(gcf, fullfile(script_dir, "tracking_scenario.pdf"));

num_methods = length(constr_list);
line_styles_rmse = {'--'};
method_colors = [0.4660 0.6740 0.1880];
true_color = [1, 1, 1] * 0.2;

figure; clf; set(gcf, 'Position', [100 100 900 600]);
font_name = 'Times New Roman';
font_size = 18;
line_width = 2;
t = tiledlayout(2, 2, 'TileSpacing', 'tight', 'Padding', 'tight');
t.TileSpacing = 'compact';
t.Padding = 'compact';

vars = {'$x$(m)', '$y$(m)', '$\theta$(rad)', '$v$(m/s)'};
var_idx = [1 2 3 4];
x_true(3, :) = unwrap(x_true(3, :));
x_est_store{1, 1}(3, :) = unwrap(x_est_store{1, 1}(3, :));

for idx = 1:4
    ax = nexttile(idx); hold on;

    for gg = 1:num_methods
        x_est = x_est_store{1, gg};
        if ~isempty(x_est)
            plot(1:(K + 1), x_est(var_idx(idx), :), line_styles_rmse{gg}, ...
                'LineWidth', line_width, 'Color', method_colors(gg, :), ...
                'DisplayName', constr_legend(gg));
        end
    end

    plot(1:(K + 1), x_true(var_idx(idx), :), 'color', true_color, ...
        'LineWidth', line_width, 'DisplayName', 'Ground Truth');
    ylabel(vars{idx}, 'Interpreter', 'latex', 'FontName', font_name, 'FontSize', font_size);
    xlim([1, K + 1]);
    grid on;
    set(gca, 'FontName', font_name, 'FontSize', font_size);
    if idx >= 3
        xlabel('Step', 'FontName', font_name, 'FontSize', font_size + 1);
    end
    if idx == 1
        leg = legend('show', 'Location', 'northoutside', 'Orientation', 'horizontal', ...
            'FontName', font_name, 'FontSize', font_size + 1);
        leg.Box = 'off';
    end
end

outerpos = t.OuterPosition;
outerpos(4) = outerpos(4) * 0.97;
t.OuterPosition = outerpos;

exportgraphics(gcf, fullfile(script_dir, "tracking_rmse.pdf"));
