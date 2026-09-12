%% Spillover compensation for Imaging Mass Cytometry (SOcorrection_all_csv_ALS_HyperionXTi.m)
% -------------------------------------------------------------------------
% Platform   : HyperionXTi (Standard BioTools)
% Dataset    : ALS IMC dataset
% -------------------------------------------------------------------------
% Author     : Gen Ohtsuki
% Affiliation: Kyoto University 
% Since 2026-06-01
% Modified 2026-07-28
% Copyright (c) 2026 Gen Ohtsuki
%
% Associated publication:
% Hikosaka et al. in review
% "Glial metabolic states in human sporadic amyotrophic lateral sclerosis
% uncovered by spatial single-cell proteomics."
% -------------------------------------------------------------------------
% Description:
%   Iterates over all '*_SOc.csv' files in folderPath and calls
%   "Spillover_compensation_IMC_GO_func.m" for each to perform spillover correction.
% -------------------------------------------------------------------------
% Input:
%   spillover_matrix  : ChannelN x ChannelN spillover matrix (%)
%   IMCmeasured       : cellN x ChannelN matrix of single-cell intensity values
% Output:
%   compensated       : cellN x ChannelN matrix of spillover-compensated values
% -------------------------------------------------------------------------
% Requirements:
%    - Signal Overlap Matrix_ALS_Hyperion.xlsx
%    - Spillover_compensation_ALS_Hyperion_func.m
%    - keep4.m
% -------------------------------------------------------------------------
%   This function is called from "SOcorrection_all_csv_ALS_HyperionXTi.m",
%   which iterates over all files matching '*_SOc.csv' in folderPath.
% -------------------------------------------------------------------------
clear all; 
close all;
warning('off');
scrsz = get(0,'ScreenSize');
curdir = pwd;

% % Cell Compartment % %
Cell_Compartment = "cytosol" %{"cytosol", "nucleus", "torus"}


FolderPath = '[path]';
addpath(FolderPath);

folderPath = fullfile(FolderPath, 'IMC_csv');
SOCorrection_folderPath = FolderPath;


%% Process "Spillover_compensation_ALS_HyperionXTi_func.m" in all .csv files.
SOc_files = dir(fullfile(folderPath, '*_SOc.csv'));

if isempty(SOc_files)
    error('No files matching "*_SOc.csv" were found in:\n%s', folderPath);
end

for k = 1:length(SOc_files)

    inputFile_Name  = SOc_files(k).name;
    outputFile_Name = SOc_files(k).name;

    Spillover_compensation_ALS_HyperionXTi_func
    fprintf('Processed: %s\n', inputFile_Name);
end

disp('All process was done.');