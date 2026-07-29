clear all

% Define user folder path
userfolder = 'C:\Users\tkx1376';
% Define paper folder path
paperfolder = '\OneDrive - Northwestern University\Documents - Backman Lab - Shared Folders\Lab Paper Drafts\HP1a Paper';
% Concatenate user folder and paper folder paths
parentfolder = strcat(userfolder,paperfolder);

% Define CSV folder and file names
csvfolder = '\OneDrive - Northwestern University\Documents - Backman Lab - Shared Folders\Lab Paper Drafts\HP1a Paper\Data\Auxin PWS Data\PWS Data\PWS Dynamics Measurements';
csvfile = "090924_HP1a_AID_PWS.csv";
outfile = "090924_HP1a_AID_PWS_Dynamics.csv";

% Concatenate paths to create full CSV file paths
csvName = strcat(parentfolder,csvfolder,csvfile);
outName = strcat(parentfolder,csvfolder,outfile);

% The RMS measurement in background regions where true RMS should be 0
noiseRms = 0.007; 

% Add folder with functions to convert RMS to D to MATLAB path
functionpath = '\OneDrive - Northwestern University\Documents - Backman Lab - Shared Folders\Lab Paper Drafts\HP1a Paper\Sigma to D\PWS-SigmaToD';
addpath(genpath(strcat(parentfolder,functionpath)));

% Crowding volume fraction
phi = 0.37; % A549 CVC value from Li et al. 2022
% Genomic size of a packing domain in base pairs (bp)
Nf = 207e3; % A549 N_f value from Li et al. 2022
% Cell thickness in um
thickness = 2;

% Refractive index of the substance the cell is submerged in
RI = 1.337; % RI of cell culture media
% Illumination NA (determined by size of A stop)
NAi = 0.52; % Measured using NA calculator or hexagon size
% Collection NA (written on objective)
NAc = 1.49; % NA of 100x objective
% Central imaging wavelength in nm
lambda = 585; 
% The type of objective used to image (change to false if using air objective)
oil_objective = true; % 60x, 63x, and 100x objectives are oil immersion
% The z location used when imaging (change to false if imaging was done above the glass surface)
cell_glass = true; % PWS images were taken where the cell contacts the dish

% Create refractive index definition for live cell
liveCellRI = S2D.RIDefinition.createFromGladstoneDale(RI, phi);
% Create system configuration for SigmaToD conversion
nuSys = S2D.SystemConfiguration(liveCellRI, NAi, NAc, lambda, oil_objective, cell_glass);

% Detect import options for the CSV file
opts = detectImportOptions(csvName, 'VariableNamingRule', 'preserve');
% Read the CSV file into a table
myTable = readtable(csvName, opts);
% Extract RMS values from the table
rms = table2array(myTable(:, "RMS"));

% Correct RMS values by subtracting noise RMS
rms = sqrt(rms.^2 - noiseRms^2);

% Convert RMS to diffusion coefficient D
[dOut,dCorrected, Nf_expected,lmax_corrected] = SigmaToD_AllInputs(rms, nuSys, Nf, thickness);

% Add corrected D values to the table
myTable = [myTable table(dCorrected, 'VariableNames', {'D'})];
% Write the updated table to a new CSV file
writetable(myTable, outName);