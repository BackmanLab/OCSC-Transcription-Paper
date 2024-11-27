%% Create LUT for sigma to D conversion
% Jane Frederick 5-20-2024

% Clear variables to prevent errors
clear variables

%% User Inputs

% ALL USERS MUST CHANGE THIS
% The file path should lead to the folder with the sigma to D conversion
functionpath = 'C:\Users\janef\OneDrive - Northwestern University\Documents - Backman Lab - Shared Folders\Subgroup Folders\Chromatin\Cancer Stem Cells\OCSC Chromatin Transcription Paper 2024\Codes\Sigma to D';
addpath(genpath(functionpath));

%% Inputs for Constants on LCPWS Systems

% Crowding volume fraction
phi = 0.37; % A549 CVC value from Li et al. 2022
% Genomic size of a packing domain in base pairs (bp)
Nf = 207e3; % A549 N_f value from Li et al. 2022
% Cell thickness in um
thickness = 2;

% Refractive index of the substance the cell is submerged in
RI = 1.337; % RI of cell culture media
% Illumination NA (determined by size of A stop)
NAi_LCPWS1 = 0.54;
NAi_LCPWS2 = 0.52; % Measured using NA calculator or hexagon size
% Collection NA (written on objective)
NAc_LCPWS1 = 1.4; % NA of 63x objective
NAc_LCPWS2 = 1.49; % NA of 100x objective
% Central imaging wavelength in nm
lambda = 585; 
% The type of objective used to image (change to false if using air
% objective)
oil_objective = true; % 60x, 63x, and 100x objectives are oil immersion
% The z location used when imaging (change to false if imaging was done
% above the glass surface)
cell_glass = true; % PWS images were taken where the cell contacts the dish

%% Determine relevant microscope parameters

% Use the Gladstone-Dale equation for refractive index
liveCellRI = S2D.RIDefinition.createFromGladstoneDale(RI, phi);
% Create configuration for LCPWS1 microscope
LCPWS1Sys = S2D.SystemConfiguration(liveCellRI, NAi_LCPWS1, NAc_LCPWS1, lambda, oil_objective, cell_glass);
% Create configuration for LCPWS2 microscope
LCPWS2Sys = S2D.SystemConfiguration(liveCellRI, NAi_LCPWS2, NAc_LCPWS2, lambda, oil_objective, cell_glass);

%% Convert the RMS values in the CSV files

% Generate a range of RMS values
rms = linspace(0,0.497);
% Convert RMS values to D values for LCPWS1 system
[dOut_LCPWS1,dCorrected_LCPWS1, Nf_expected_LCPWS1,lmax_corrected_LCPWS1] = SigmaToD_AllInputs(rms, LCPWS1Sys, Nf, thickness);
% Convert RMS values to D values for LCPWS2 system
[dOut_LCPWS2,dCorrected_LCPWS2, Nf_expected_LCPWS2,lmax_corrected_LCPWS2] = SigmaToD_AllInputs(rms, LCPWS2Sys, Nf, thickness);
% Create lookup table for LCPWS1 system
convlut_LCPWS1 = [rms;dOut_LCPWS1;dCorrected_LCPWS1;Nf_expected_LCPWS1;lmax_corrected_LCPWS1];
% Create lookup table for LCPWS2 system
convlut_LCPWS2 = [rms;dOut_LCPWS2;dCorrected_LCPWS2;Nf_expected_LCPWS2;lmax_corrected_LCPWS2];

% Convert lookup tables to tables and write to CSV files
conversion_table_LCPWS1 = array2table(convlut_LCPWS1','VariableNames',{'RMS','Db','D','Nf','lmax'});
conversion_table_LCPWS2 = array2table(convlut_LCPWS2','VariableNames',{'RMS','Db','D','Nf','lmax'});
writetable(conversion_table_LCPWS1,'SigmaToDLUT_LCPWS1.csv');
writetable(conversion_table_LCPWS2,'SigmaToDLUT_LCPWS2.csv');