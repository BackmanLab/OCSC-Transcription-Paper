%% Convert all RMS values from several folders to D
% Jane Frederick 5-20-2024

% Clear variables to prevent errors
clear variables

%% User Inputs

% ALL USERS MUST CHANGE THIS
% The file path should lead to the folder with the sigma to D conversion
functionpath = 'C:\Users\janef\OneDrive - Northwestern University\Documents - Backman Lab - Shared Folders\Subgroup Folders\Chromatin\Cancer Stem Cells\OCSC Chromatin Transcription Paper 2024\Codes\Sigma to D\PWS-SigmaToD';
addpath(genpath(functionpath));

% ALL USERS MUST CHANGE THIS
% Ask user for main folder with all experiments to be converted
%folderdir = uigetdir('C:\','Select Experiment Folder');
% Use this line if you don't want to select the folder manually
folderdir = 'C:\Users\janef\OneDrive - Northwestern University\Documents - Backman Lab - Shared Folders\Subgroup Folders\Chromatin\Cancer Stem Cells\OCSC Chromatin Transcription Paper 2024\Data\PWS RMS and D CSV Files';

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

% Find all CSV files to convert
fileinfo = dir([folderdir,'\**\*.csv']);

% Count number of files to convert
numfiles = size(struct2table(fileinfo),1);

% Create a folder for converted values in the parent folder
outdir = [folderdir,'\..\D Files'];
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% Loop through all folders and convert RMS values to D
for i=1:numfiles
     filepath = [fileinfo(i).folder,'\',fileinfo(i).name];
     opts = detectImportOptions(filepath, 'VariableNamingRule', 'preserve');
     myTable = readtable(filepath, opts);

     % Check if the 'RMS' column exists, if not rename the first column to 'RMS'
     if(isequal(any(strcmp('RMS', myTable.Properties.VariableNames)),0))
          myTable.Properties.VariableNames{1} = 'RMS';
     end
     rms = table2array(myTable(:, "RMS"));

     % Extract the year from the filename to determine the system used
     newstr = extract(fileinfo(i).name,"20" + digitsPattern(2));
     if(str2double(newstr{1})<2019)
          % Extra reflection correction factor for older data
          rms = rms.*2.43; 
          nuSys = LCPWS1Sys;
          % The RMS measurement in background regions where true RMS should be 0
          noiseRms = 0.009;
     else
          nuSys = LCPWS2Sys;
          % The RMS measurement in background regions where true RMS should be 0
          noiseRms = 0.1;
     end

     % Noise subtraction
     rms = sqrt(rms.^2 - noiseRms^2);

     % Convert RMS to D using the SigmaToD_AllInputs function
     [dOut,dCorrected, Nf_expected,lmax_corrected] = SigmaToD_AllInputs(rms, nuSys, Nf, thickness);
     myTable = [myTable table(dCorrected, 'VariableNames', {'D'}) table(Nf_expected, 'VariableNames', {'Nf'}) table(lmax_corrected, 'VariableNames', {'lmax'})];

     % Create output folder if it doesn't exist
     outputfold = [outdir,extractAfter(fileinfo(i).folder,folderdir)];
     if ~exist(outputfold, 'dir')
         mkdir(outputfold);
     end

     % Write the converted table to a new CSV file
     [~,name,ext] = fileparts(filepath);
     writetable(myTable, strcat(outputfold,'\',name,'_D.csv'));

     % Print progress to the console
     fprintf('%1.0f out of %1.0f files processed.\n',i,numfiles);
end
% Print completion message to the console
fprintf('RMS to D conversion for all %1.0f files completed.\n',numfiles);