%% FOV_correction
% Description:
% Spinal cord white matter filtering
% This script overlays IMC cell coordinates on a TIFF image and allows
% manual freehand selection of cells to exclude from downstream analysis.
%
% Input:
%   - TIFF image for ROI visualization
%   - CSV file containing single-cell coordinates
%
% Output:
%   - corrected CSV file
%        Sel01:
%        Flag indicating whether a row is retained or excluded from downstream analysis.
%
% Copyright (c) 2025 Momoka Hikosaka
% This code is provided for academic research purposes.

clear all; 
close all; warning('off');
scrsz = get(0,'ScreenSize');
curr_dir = pwd;
xsc_gap = scrsz(3)*20/64;

%% ===== Settings =====
base_dir = curr_dir

data_dir = fullfile(base_dir, "raw_csv");
output_dir = fullfile(base_dir, "output");
if ~exist(output_dir, 'dir')
    mkdir(output_dir)
end

file_name = "SpinalCord_02_ALS";

input_file = fullfile(data_dir, file_name + ".csv");
tiff_image = fullfile(data_dir, "Ir_image_" + file_name + ".tiff");


%% ===== Read TIFF image =====
[tiff_roi, ~] = imread(tiff_image, 1);

IMC_szY = size(tiff_roi, 1);
IMC_szX = size(tiff_roi, 2);
IMC_FOV_size = [IMC_szX, IMC_szY]; %#ok<NASGU>

%% ===== Read CSV table =====
opts = detectImportOptions(input_file);
opts.VariableNamesLine = 1;
opts.DataLines = [2 Inf];
T = readtable(input_file, opts);

required_vars = {'CellId','X_position','Y_position'};
assert(all(ismember(required_vars, T.Properties.VariableNames)), ...
    'CSV must contain CellId, X_position, and Y_position.');

% Initialize Sel01 if missing
if ~ismember('Sel01', T.Properties.VariableNames)
    T.Sel01 = ones(height(T), 1);
end

x = T.X_position;
y = T.Y_position;

%% ===== Save workspace =====
save(fullfile(output_dir, 'Cell_ident.mat'), '-v7.3');

%% ===== Display overlay and draw exclusion region =====
fig_name = "Croped";
fig = figure('Name', fig_name, ...
    'Position', [scrsz(3)*1/30, scrsz(4)*1/20, ...
                 scrsz(4)*0.5*(IMC_szX/IMC_szY), scrsz(4)*0.5]);

ax = axes(fig);
imshow(tiff_roi(:,:,1:min(3,size(tiff_roi,3))), 'Parent', ax);
hold(ax, 'on');

scatter(ax, x, y, 30, ...
    'MarkerEdgeColor', 'r', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 0.5);

xlim(ax, [0 IMC_szX]);
ylim(ax, [0 IMC_szY]);
axis(ax, 'ij');

% Draw freehand exclusion region
h = drawfreehand(ax, 'Color', [0 1 0]);
wait(h);

pos = h.Position;
xv = pos(:,1);
yv = pos(:,2);

%% ===== Exclude cells inside drawn region =====
in_region = inpolygon(T.X_position, T.Y_position, xv, yv);
T.Sel01(in_region & T.Sel01 == 1) = 0;

%% ===== Redraw final retained cells =====
cla(ax);
imshow(tiff_roi(:,:,1:min(3,size(tiff_roi,3))), 'Parent', ax);
hold(ax, 'on');

scatter(ax, T.X_position(T.Sel01 == 1), T.Y_position(T.Sel01 == 1), 30, ...
    'MarkerEdgeColor', 'r', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 0.5);

plot(ax, xv, yv, 'Color', [0 1 0], 'LineWidth', 1);

xlim(ax, [0 IMC_szX]);
ylim(ax, [0 IMC_szY]);
axis(ax, 'ij');

%% ===== Save figure and output table =====
saveas(fig, fullfile(output_dir, fig_name), 'fig');
saveas(fig, fullfile(output_dir, fig_name), 'bmp');

out_csv = fullfile(output_dir, file_name + "_crop.csv");
writetable(T, out_csv);

disp("Processing completed.");