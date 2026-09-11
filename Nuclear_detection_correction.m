%% Nuclear_detection_correction
% Description:
% This script identifies large Ir signals from a TIFF image, 
% overlays cell coordinates from a CSV file, and
% manually selects cells for correction.
%
% Selected cells can be:
%   1) replaced by an averaged row
%   2) removed without averaging
%
% Input:
%   - TIFF image for ROI visualization
%   - CSV file containing single-cell coordinates and marker values
%
% Output:
%   - corrected CSV file
%        Sel01:
%        Flag indicating whether a row is retained or excluded from downstream analysis.
%        Selection:
%        Record of the selection group indicating which manual selection round processed the cell.
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

%% ===== Detect large objects from TIFF =====
grayImg = tiff_roi;
if size(grayImg, 3) == 3
    grayImg = rgb2gray(grayImg);
end

grayImg = mat2gray(grayImg);
grayImg = imgaussfilt(grayImg, 1.2);

bw = grayImg > 0.5;
bw = bwareaopen(bw, 5);

CC = bwconncomp(bw);
stats = regionprops(CC, "Area", "Centroid", "EquivDiameter");

large_idx = find([stats.Area] > 300);

fig_name = "Large_cell_selection";
f = figure('Name', fig_name, 'Position', [300 300 800 800]);
imshow(tiff_roi(:,:,1:min(3,size(tiff_roi,3))));
hold on;

auto_centers = vertcat(stats(large_idx).Centroid);
auto_radii   = [stats(large_idx).EquivDiameter] * 1.5;

for i = 1:numel(large_idx)
    viscircles(auto_centers(i,:), auto_radii(i), 'Color', 'g', 'LineWidth', 0.5);
end

disp('Click where you want to add it and press Enter to finish.');

manual_centers = [];
manual_radii   = [];
manual_radius_default = 50;

while true
    [xclick, yclick] = ginput(1);
    if isempty(xclick)
        break;
    end
    manual_centers = [manual_centers; xclick, yclick]; %#ok<AGROW>
    manual_radii   = [manual_radii; manual_radius_default]; %#ok<AGROW>
    viscircles([xclick, yclick], manual_radius_default, 'Color', 'y', 'LineWidth', 0.5);
end

hold off;
close(f);

all_centers = [auto_centers; manual_centers];
all_radii   = [auto_radii(:); manual_radii(:)];

%% ===== Read CSV table =====
opts = detectImportOptions(input_file);
opts.VariableNamesLine = 1;
opts.DataLines = [2 Inf];
T = readtable(input_file, opts);

required_vars = {'CellId','X_position','Y_position'};
assert(all(ismember(required_vars, T.Properties.VariableNames)), ...
    'CSV must contain CellId, X_position, and Y_position.');

flag_var = 'Sel01'; % Flag indicating whether a row is retained (1) or excluded (0).
selection_var = 'Selection'; % Identifier of the manual selection group used during the correction process.

if ~ismember(flag_var, T.Properties.VariableNames)
    T.(flag_var) = ones(height(T), 1);
end

if ~ismember(selection_var, T.Properties.VariableNames)
    T.(selection_var) = zeros(height(T), 1);
end

idx_use = (T.(flag_var) == 1);
row_idx_use = find(idx_use);

T_ex = T(idx_use, {'CellId','X_position','Y_position'});
x = T_ex.X_position;
y = T_ex.Y_position;

%% ===== Save workspace =====
save(fullfile(output_dir, 'Cell_ident.mat'), '-v7.3');

%% ===== Initial overlay figure =====
fig_name = "Cell_ident";
fig1 = figure('Name', fig_name, ...
    'Position', [scrsz(3)*1/30, scrsz(4)*1/20, ...
                 scrsz(4)*0.5*(IMC_szX/IMC_szY), scrsz(4)*0.5]);

imshow(tiff_roi(:,:,1:min(3,size(tiff_roi,3))));
hold on;

scatter(x, y, 30, ...
    'MarkerEdgeColor', 'r', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 0.5);

for k = 1:numel(all_radii)
    viscircles(all_centers(k,:), all_radii(k), 'Color', [0 1 0], 'LineWidth', 0.6);
end

xlim([0 IMC_szX]);
ylim([0 IMC_szY]);
hold off;

saveas(fig1, fullfile(output_dir, fig_name), 'fig');
saveas(fig1, fullfile(output_dir, fig_name), 'bmp');

%% ===== Manual correction figure =====
fig_name = "Data_correction";
fig2 = figure('Name', fig_name, ...
    'ToolBar', 'figure', ...
    'Position', [scrsz(3)*1/30 + scrsz(4)*(IMC_szX/IMC_szY), ...
                 scrsz(4)*1/20, ...
                 scrsz(4)*0.5*(IMC_szX/IMC_szY), scrsz(4)*0.5]);

hIm = imshow(tiff_roi(:,:,1:min(3,size(tiff_roi,3))));
set(hIm, 'HitTest', 'off');
hold on;

hScatter = scatter(x, y, 30, ...
    'MarkerEdgeColor', 'flat', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 0.5);

xlim([0 IMC_szX]);
ylim([0 IMC_szY]);

C = repmat([1 0 0], numel(x), 1);
set(hScatter, 'CData', C);

findClosestPoint = @(X,Y,p) ...
    find(((X - p(1)).^2 + (Y - p(2)).^2) == min((X - p(1)).^2 + (Y - p(2)).^2), 1);

setappdata(hScatter, 'selected', false(numel(x),1));
setappdata(hScatter, 'everSelected', false(numel(x),1));

keepGoing = true;
roundIdx = 1;

%% ===== Manual selection loop =====
while keepGoing
    fprintf('--- Round %d ---\n', roundIdx);

    set(fig2, 'WindowButtonDownFcn', @(src,evt) localClick(src, evt, hScatter, findClosestPoint));
    set(fig2, 'WindowKeyPressFcn', @(src,evt) onKey(src, evt));

    uiwait(fig2);

    sel  = getappdata(hScatter, 'selected');
    ever = getappdata(hScatter, 'everSelected');

    if ~any(sel)
        disp('No selection. Skipped this round.');

        C = repmat([1 0 0], numel(ever), 1);
        C(ever,:) = repmat([0 0 1], nnz(ever), 1);
        set(hScatter, 'CData', C, 'MarkerEdgeColor', 'flat');

        nextAns = questdlg('No points selected. Continue to next round?', ...
                           'Continue', 'Yes', 'No', 'Yes');
        if strcmp(nextAns, 'Yes')
            setappdata(hScatter, 'selected', false(numel(x),1));
            roundIdx = roundIdx + 1;
            continue;
        else
            keepGoing = false;
            break;
        end
    end

    cx = mean(x(sel), 'omitnan');
    cy = mean(y(sel), 'omitnan');

    if isfinite(cx) && isfinite(cy)
        plot(cx, cy, 'x', 'MarkerSize', 8, 'LineWidth', 1.2, 'Color', 'g');
        text(cx + 8, cy, sprintf('%d', roundIdx), ...
            'Color', 'g', 'FontSize', 8, 'FontWeight', 'bold');
        drawnow;
    end

    ever = ever | sel;
    setappdata(hScatter, 'everSelected', ever);

    C = repmat([1 0 0], numel(ever), 1);
    if any(ever)
        C(ever,:) = repmat([0 0 1], nnz(ever), 1);
    end
    set(hScatter, ...
        'CData', C, ...
        'MarkerEdgeColor', 'flat', ...
        'MarkerFaceColor', 'flat', ...
        'MarkerFaceAlpha', 'flat', ...
        'SizeData', 30);

    alphaData = zeros(numel(sel), 1);
    alphaData(ever | sel) = 1;
    set(hScatter, 'AlphaData', alphaData, 'AlphaDataMapping', 'none');

    answer = questdlg('Add a new averaged row for selected cells?', ...
                      'Average', 'Yes', 'No', 'Yes');

    if strcmp(answer, 'Yes')
        mean_cols = [1, 3:53];

        dataSel = T{row_idx_use(sel), mean_cols};
        meanRowPart = mean(dataSel, 1, 'omitnan');

        newRow = array2table(roundIdx * ones(1, width(T)), ...
            'VariableNames', T.Properties.VariableNames);
        newRow{1, mean_cols} = meanRowPart;
        newRow.(flag_var) = 1;

        T = [T; newRow];

        T.(selection_var)(row_idx_use(sel)) = roundIdx;
        T.(flag_var)(row_idx_use(sel)) = 0;

        disp('Averaged row added successfully.');

    else
        T.(selection_var)(row_idx_use(sel)) = NaN;
        T.(flag_var)(row_idx_use(sel)) = 0;

        disp('Selected cells removed without averaging.');
    end

    nextAns = questdlg('Continue to next round?', 'Continue', 'Yes', 'No', 'Yes');
    if strcmp(nextAns, 'Yes')
        setappdata(hScatter, 'selected', false(numel(x),1));
        roundIdx = roundIdx + 1;
    else
        keepGoing = false;
    end

    saveas(fig2, fullfile(output_dir, fig_name), 'fig');
    saveas(fig2, fullfile(output_dir, fig_name), 'bmp');
end

%% ===== Exclude invalid rows =====
outOfBounds = isnan(T.X_position) | isnan(T.Y_position) | ...
              T.X_position < 0 | T.X_position > IMC_szX | ...
              T.Y_position < 0 | T.Y_position > IMC_szY;

T.(flag_var)(outOfBounds) = 0;

hasNaN = any(ismissing(T), 2);
T.(flag_var)(hasNaN) = 0;

%% ===== Save corrected table =====
out_csv = fullfile(output_dir, file_name + "_m.csv");
writetable(T, out_csv);

%% ===== Final overlay figure =====
fig_name = "Cell_ident_final";
fig3 = figure('Name', fig_name, ...
    'Position', [scrsz(3)*1/30, scrsz(4)*1/20, ...
                 scrsz(4)*0.5*(IMC_szX/IMC_szY), scrsz(4)*0.5]);

imshow(tiff_roi(:,:,1:min(3,size(tiff_roi,3))));
hold on;

valid_idx = (T.(flag_var) == 1);
x_valid = T.X_position(valid_idx);
y_valid = T.Y_position(valid_idx);

scatter(x_valid, y_valid, 30, ...
    'MarkerEdgeColor', 'r', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 0.5);

for k = 1:numel(all_radii)
    viscircles(all_centers(k,:), all_radii(k), 'Color', [0 1 0], 'LineWidth', 0.6);
end

xlim([0 IMC_szX]);
ylim([0 IMC_szY]);
hold off;

saveas(fig3, fullfile(output_dir, fig_name), 'fig');
saveas(fig3, fullfile(output_dir, fig_name), 'bmp');

%% ===== Local functions =====
function onKey(src, evt)
    switch evt.Key
        case {'return','enter','escape'}
            uiresume(ancestor(src, 'figure'));
        case 'z'
            fig = ancestor(src, 'figure');
            z = zoom(fig);
            if strcmp(z.Enable, 'on')
                zoom(fig, 'off');
                disp('Click Select Mode');
            else
                zoom(fig, 'on');
                disp('Zoom Mode');
            end
    end
end

function localClick(src, ~, hScatter, findClosestPoint)
    if ~strcmp(get(src, 'SelectionType'), 'normal')
        return;
    end

    axH = ancestor(hScatter, 'axes');
    obj = hittest(src);
    axOfClick = ancestor(obj, 'axes');

    if isempty(axOfClick) || axOfClick ~= axH
        return;
    end

    cp = get(axH, 'CurrentPoint');
    p = cp(1,1:2);

    xlim_ = get(axH, 'XLim');
    ylim_ = get(axH, 'YLim');

    if p(1) < xlim_(1) || p(1) > xlim_(2) || p(2) < ylim_(1) || p(2) > ylim_(2)
        return;
    end

    idx = findClosestPoint(hScatter.XData, hScatter.YData, p);

    sel = getappdata(hScatter, 'selected');
    if isempty(sel) || numel(sel) ~= numel(hScatter.XData)
        sel = false(numel(hScatter.XData), 1);
    end

    sel(idx) = ~sel(idx);
    setappdata(hScatter, 'selected', sel);

    ever = getappdata(hScatter, 'everSelected');
    tempMask = ever | sel;

    C = repmat([1 0 0], numel(sel), 1);
    C(tempMask,:) = repmat([0 0 1], nnz(tempMask), 1);

    set(hScatter, ...
        'CData', C, ...
        'MarkerEdgeColor', 'flat', ...
        'MarkerFaceColor', 'flat', ...
        'MarkerFaceAlpha', 'flat', ...
        'SizeData', 30);

    alphaData = zeros(numel(sel),1);
    alphaData(tempMask) = 1;
    set(hScatter, 'AlphaData', alphaData, 'AlphaDataMapping', 'none');

    drawnow;
    close all;
end

disp("Processing completed.");