%% function [cubeRms, BW] = load_pws(folder, analysis, roi)
%
% DESCRIPTION
%   This function loads RMS files and nuclear ROI files created by the PWS
%   analysis codes.
%
% INPUT ARGUMENTS
%   folder:
%       Optional string input to define the path to the cell folder that
%       the user wishes to load. If no input is given, a dialog box will
%       open asking the user to choose a folder.
%   analysis:
%       Optional input to define the analysis name if multiple analyses
%       were run on the same data. If no input is given and only 1 analysis
%       file is available, the code will default to that analysis file. If
%       multiple analyses files are in the folder, a dialog box will open
%       asking the user which analysis file they wish to load.
%   roi:
%       Optional string input to define the ROI name if multiple ROIs were
%       drawn on the same data. If no input is given and only 1 roi file is
%       available, the code will default to that ROI file. If multiple roi
%       files are in the fodler, a dialog box will open asking the user
%       which roi file they would like to load. If no ROI files are found,
%       then no ROIs will be loaded (the code will still run).
%
% OUTPUT ARGUMENTS
%   cubeRms:
%       The RMS file defined by the user, calculated from the raw PWS data
%       in the folder. cubeRms should be a 2D array.
%   BW:
%       A stack of binary images defining ROIs loaded from the file defined
%       by the user. The first 2 dimensions are equivalent to X and Y
%       dimanesion of cubeRms, while the third dimension is the ROI number.
%
% EXAMPLES
%   [cubeRms BW] = load_pws('/Users/anneshim/Desktop/WT/Cell1', 'p0', 'nuc')
%   [cubeRms BW] = load_pws([], [], 'nuc_2')
%   [cubeRms BW] = load_pws
%   cubeRms = load_pws
%
% REFERENCES
%   https://doi.org/10.1101/2020.01.26.920363 
%
% Author: Adam Eshein (aeshein@u.northwestern.edu) 5.10.2020
function [cubeRms, BW] = load_pws(folder, analysis, roi)
%% Check to make sure inouts are correct and get the cell folder
if nargin < 1 || isempty(folder)
    folder = uigetdir([],'Choose a Cell folder to load.');
elseif nargin > 3
    error('NNET:Arguments', 'Wrong number of inputs. This function requires 0, 1, 2 or 3 inputs.')
end

%% load RMS
if isdir([folder '/PWS/analyses'])
    if nargin >= 2 && ~isempty(analysis)
        cubeRms = h5read([folder '/PWS/analyses/analysisResults_' analysis '.h5'],'/rms');
    elseif size(ls([folder '/PWS/analyses/analysisResults_*']),1) == 1
        cubeRms = h5read([folder '/PWS/analyses/' ls([folder '/PWS/analyses/analysisResults_*'])],'/rms');
    elseif size(ls([folder '/PWS/analyses/analysisResults_*']),1) == 0
        error('No analysis files found.')
    else
        file_rms = uigetfile([folder '/PWS/analyses/.h5'],'Choose an analysis file to load.');
        cubeRms = h5read([folder '/PWS/analyses/' file_rms],'/rms');
    end
else
    error('No "analyses" folder found.')
end

%% load BWs

if nargout > 1 % don't bother loading BW if there's no output to assign it to
    
    if nargin ==3 && ~isempty(roi)
        hinfo = hdf5info([folder '/ROI_' roi '.h5']);
    elseif size(ls([folder '/ROI*']),1) == 1
        hinfo = hdf5info([folder '/' ls([folder '/ROI*'])]);
    elseif size(ls([folder '/ROI*']),1) == 0
        fprintf('\n No ROI file was found. \n\n ')
        hinfo = [];
        BW = [];
    else
        file_bw = uigetfile([folder '/.h5'],'Choose an ROI file to load.');
        hinfo = hdf5info([folder '/' file_bw]);
    end
    
    if ~isempty(hinfo)
        for ii = 1:length(hinfo.GroupHierarchy.Groups)
            BW(:,:,ii) = hdf5read(hinfo.GroupHierarchy.Groups(ii).Datasets(1));
        end
    end
else
    fprintf('No ROI file was loaded, since the user did not define an output for it.')
end