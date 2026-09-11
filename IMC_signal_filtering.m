%% IMC_signal_filtering.m
% Description:
% Manual upper-threshold filtering to remove putative signal detection errors
% in marker expression values.
% Cells with expression values above a manually selected threshold
% are replaced with NaN in the output table.
%
% Input:  *.csv (raw IMC single-cell data)
% Output: *_filtered_NaN.csv
%
% Copyright (c) 2025 Momoka Hikosaka
% This code is provided for academic research purposes.

clear all; close all; warning off;
scrsz = get(0,'ScreenSize');
curr_dir = pwd;

%% ===== Settings =====
base_dir = curr_dir

data_dir = fullfile(base_dir, "raw_csv");
output_dir = fullfile(base_dir, "output");
if ~exist(output_dir, 'dir')
    mkdir(output_dir)
end

file_name = "PrecentralGyrus_01_ALS";

input_file = fullfile(data_dir, file_name + ".csv");

% Marker columns to process
marker_cols = 3:32;

%% ===== Read input table =====
data_table = readtable(input_file);
column_names = data_table.Properties.VariableNames;

% Initialize output tables
data_table_filtered_NaN = data_table;
marker_names = strings(numel(marker_cols), 1);

%% ===== Manual filtering loop =====
for k = 1:numel(marker_cols)

    col_idx = marker_cols(k);
    marker_name = string(column_names{col_idx});
    marker_names(k) = marker_name;

    % Original expression values
    marker_exp_val = data_table{:, col_idx};

    % Working variables
    marker_exp_val_NaN = marker_exp_val;
    marker_exp_val_max = marker_exp_val;

    % Initial maximum before thresholding
    current_max = max(marker_exp_val, [], "omitnan");

    % Show original plot
    fig_raw = figure;
    plot(marker_exp_val, 'o-');
    xlabel('Cell ID');
    ylabel('Expression value');
    title(marker_name, 'Interpreter', 'none');

    % Repeated manual thresholding
    while true
        disp(['Marker: ' char(marker_name)]);
        disp('Click a y-threshold on the figure, or press Enter to finish.');

        [x, threshold] = ginput(1); %#ok<ASGLU>

        if isempty(x)
            break;
        end

        hold on;
        yline(threshold, 'r--', 'Manual Threshold', ...
            'LabelHorizontalAlignment', 'left');

        % Replace values above threshold with NaN
        marker_exp_val_NaN(marker_exp_val > threshold) = NaN;

        % Compute max of remaining values
        current_max = max(marker_exp_val_NaN, [], "omitnan");

        % Replace removed values with filtered max
        marker_exp_val_max = marker_exp_val_NaN;
        marker_exp_val_max(isnan(marker_exp_val_max)) = current_max;

        % Display updated plot
        fig_updated = figure;
        plot(marker_exp_val_max, 'o-');
        xlabel('Cell ID');
        ylabel('Expression value');
        title(marker_name + " (filtered)", 'Interpreter', 'none');
    end

    close(fig_raw);

    % Save filtered values back to tables
    data_table_filtered_NaN{:, col_idx} = marker_exp_val_NaN;
    
    
    % Save final figure
    fig_final = figure;
    plot(marker_exp_val_max, 'o-');
    xlabel('Cell ID');
    ylabel('Expression value');
    title(marker_name + " (final)", 'Interpreter', 'none');

    saveas(fig_final, fullfile(output_dir, file_name + "_" + marker_name + "_signal_intensity.png"));
    close(fig_final);

    close all;
end

%% ===== Save output tables =====
output_file_NaN = fullfile(output_dir, file_name + "_filtered_NaN.csv");
writetable(data_table_filtered_NaN, output_file_NaN);

disp("Processing completed.");