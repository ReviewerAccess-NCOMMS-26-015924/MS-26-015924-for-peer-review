%% Spillover compensation for Imaging Mass Cytometry (Spillover_compensation_ALS_Hyperion_func.m)
% -------------------------------------------------------------------------
% Platform   : Hyperion (Standard BioTools)
% Dataset    : ALS IMC dataset
% -------------------------------------------------------------------------
% Author     : Gen Ohtsuki
% Affiliation: Kyoto University 
% Since 2026-06-01
% Copyright (c) 2026 Gen Ohtsuki
%
% Associated publication:
% Hikosaka et al. in review
% "Glial metabolic states in human sporadic amyotrophic lateral sclerosis
% uncovered by spatial single-cell proteomics."
% -------------------------------------------------------------------------
% Description:
%   This function is called from "SOcorrection_all_csv_ALS_Hyperion.m",
%   which iterates over all files matching '*_SOc.csv' in folderPath.
% -------------------------------------------------------------------------
% Input:
%   spillover_matrix  : ChannelN x ChannelN spillover matrix (%)
%   IMCmeasured       : cellN x ChannelN matrix of single-cell intensity values
% Output:
%   compensated       : cellN x ChannelN matrix of spillover-compensated values
% -------------------------------------------------------------------------
% Requirements:
%    - Signal Overlap Matrix_ALS_Hyperion.xlsx
%    - keep4.m
% -------------------------------------------------------------------------

% clear all; 
keep4 folderPath SOc_files SOCorrection_folderPath inputFile_Name outputFile_Name
close all;
% clc;
warning('off');
scrsz = get(0,'ScreenSize');
curdir = pwd;

%% Read imc table
SOfilePath = fullfile(folderPath, inputFile_Name);
T = readtable(SOfilePath);

outputFile = fullfile(folderPath, outputFile_Name);
T_out = T;

channelNames = { ...
    'Cell_141Pr_AQP4', ...
    'Cell_143Nd_SQSTM1', ...
    'Cell_145Nd_UBQLN2', ...
    'Cell_147Sm_Glutaminase', ...
    'Cell_149Sm_C9orf72', ...
    'Cell_150Nd_PDL1', ...
    'Cell_151Eu_MAP2', ...
    'Cell_152Sm_ATXN2', ...
    'Cell_153Eu_ORF1', ...
    'Cell_154Sm_pTDP43', ...
    'Cell_155Gd_FoxP3', ...
    'Cell_156Gd_CH24H', ...
    'Cell_158Gd_GFAP', ...
    'Cell_159Tb_IL6ST', ...
    'Cell_160Gd_IBA1', ...
    'Cell_161Dy_CX3CR1', ...
    'Cell_162Dy_pTBK1', ...
    'Cell_163Dy_CD4', ...
    'Cell_164Dy_TMEM119', ...
    'Cell_165Ho_PD1', ...
    'Cell_166Er_GAD67', ...
    'Cell_167Er_CD8', ...
    'Cell_168Er_FUS', ...
    'Cell_169Tm_ChAT', ...
    'Cell_170Er_APOE', ...
    'Cell_171Yb_TREM2', ...
    'Cell_172Yb_PDL2', ...
    'Cell_174Yb_MHCclassII', ...
    'Cell_175Lu_TDP43', ...
    'Cell_176Yb_TBK1'};

% Extract the channels as a table
T_channels = T(:, channelNames);

% Convert to numeric matrix (Ncells × NChannels)
IMCmeasured = T_channels{:,:};

%% Read SpillOver table of Hyperion 
SOCorrection_folderPath;
SOfilePath = fullfile(SOCorrection_folderPath, 'Signal Overlap Matrix_ALS_Hyperion.xlsx');

% Read cell E4:AH33
spillover_matrix = readmatrix(SOfilePath, 'Range', 'E4:AH33');


%% Convert percentages to fractions
SOcoeff = spillover_matrix / 100;

% Check matrix conditioning
cond_number = cond(SOcoeff);
fprintf('Condition number of spillover matrix = %.2f\n', cond_number);

if cond_number > 1e6
    warning('Spillover matrix is poorly conditioned.');
end

% -------------------------------------------------------------------------
% IMC Spillover Demixing
%
% IMCmeasured is a cellN x W x 30 IMC image stack
%   spillover_matrix   : 30 x 30 spillover matrix (%)
%
% Model:
%   Measured = True * SpilloverMatrix
%
% -------------------------------------------------------------------------

[cellN, Nchannels] = size(IMCmeasured);

assert(Nchannels == size(SOcoeff,1), ...
    'Number of channels does not match spillover matrix.');


%% ------------------------------------------------------------------------
% IMC Spillover Demixing Using NNLS
% -------------------------------------------------------------------------

% method: that solves a least squares problem using:
% (i) standard least squares (normal equations, A \ b)
% (ii) ridge regression (Tikhonov regularization)

method = 'ii';

Ytranspose = transpose(IMCmeasured);
Xt = Ytranspose;

opts = optimoptions('lsqnonlin', ...
    'Algorithm',     'levenberg-marquardt', ...
    'Display',       'final', ...
    'MaxIterations', 1000, ...
    'FunctionTolerance', 1e-10);

SOcoeff;

for p = 1:cellN

    y = double(Ytranspose(:,p));
    valid = ~isnan(y); %% NaN treatment for IMC 
    x = y;

  switch method
    case 'i' % 'i' = standard least squares
        x_lsq_valid = lsqnonneg(SOcoeff(valid,valid), y(valid), opts);
    case 'ii' % 'ii' = ridge regression
        lambda = 0.05; % Regulation parameter (0.001 - 0.1)
        A_R = [SOcoeff(valid,valid); lambda*eye(nnz(valid))];
        b_R = [y(valid); zeros(nnz(valid),1)];
        limit_regulation = 0.7; % Parameter to limit the subtraction, 70%. 
        x_lsq_valid_Tik = lsqnonneg(A_R, b_R);
        %[x x_lsq_valid x_lsq_valid_Tik]
        x_lsq_valid = max(x_lsq_valid_Tik, (1-limit_regulation)*y(valid));
  end

    if nnz(valid) < Nchannels
        x(valid) = x_lsq_valid;
        x(~valid) = NaN;
        x_lsq_valid = x;
    end

    Xt(:,p) = single(x_lsq_valid);

end

demixed_image = transpose(Xt);

fprintf('NNLS demixing completed.\n');

% % Double Check % %
% Signal values
IMCmeasured(1:5,:)
demixed_image(1:5,:)
D_demixed = IMCmeasured-demixed_image;
D_demixed(1:5,:)
% Reduction rates, > 0.3
ratio = mean(demixed_image,"omitnan")./mean(IMCmeasured,"omitnan");
table((1:30)', ratio')

%% Replace columns C:AF with demixed_image
if size(demixed_image,1) ~= height(T_out)
    error('Number of rows in demixed_image does not match table height.');
end

if size(demixed_image,2) ~= Nchannels
    error(['demixed_image must have ' num2str(Nchannels) ' columns.']);
end

T_out{:,3:(Nchannels+2)} = demixed_image;

%% Save
writetable(T_out, outputFile);
fprintf('Saved:\n%s\n', outputFile);
