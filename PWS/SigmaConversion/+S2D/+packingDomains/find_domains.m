%% function [L, pws_BW, params] = find_domains(pws, BW, params)
%
% DESCRIPTION
%   This funtion takes an input PWS image and does binarization and
%   segmentation to find domains.
%
% INPUT ARGUMENTS
%   pws:
%       An image of cells collected using Partial Wave Spectroscopic
%       Microscopy (PWS). The image may be in any units (D or Sigma), but
%       should be created using codes from the Backman Lab Github page.
%   BW:
%       A binarized region of interest(ROI) selecting the nucleus of a cell
%       in the PWS image. Note that a single nucleus (I.e., a single 
%       connected image) should be used.
%   params:
%       An optional input containing parameters used for the binarization
%       and segmentation process. params should be a structure array with 
%       two optional fields of:
%
%       "threshold" which may be a single number in the same units as the
%       input PWS image and will be used as a threshold to identify
%       domains, or may be the string "global" which asks the code to
%       identify the best global threshold value (using Otsu's method) to 
%       binarize the image, or may be the string "adaptive" which uses a 
%       spatially varying threshold to identify lo domains based on local 
%       values in the image. The default value is "global".
%
%       "minimum_domain_size" which sets a lower threshold on the smallest
%       domain size the algorithm will segment. The number here should be
%       in units of pixels and is the maximum radius. The defailt value is
%       0.
%
% OUTPUT ARGUMENTS
%   L:
%       A 3D binary matrix where the first two dimensions are in the units
%       of x and y of the 2D input "pws" and the third dimension is a new
%       domains. Essentially, each layer of L is a mask for one domain.
%   pws_BW:
%       A 2D binary matrix that is a mask showing all domains. This image
%       is segmented to create L.
%   params:
%       A structure array containing fields that define various params used
%       in the analaysis.
%
% EXAMPLES
%   params.threshold = 0.15
%   params.minimum_domain_size = 5
%   [L, pws_BW, params] = find_domains(cubeRms, BW, params)
%   [L, pws_BW, params] = find_domains(cubeRms, n1)
%   L = find_domains(cubeRms, BW)
%   
% REFERENCES
%   https://doi.org/10.1101/2020.01.26.920363 
%
%
% Author: Adam Eshein (aeshein@u.northwestern.edu) 5.8.2020
function [L, pws_BW, params] = find_domains(pws, BW, params)
%% Check the input arguments and assign parameters for the code

pws = double(pws); % make this a double to avoid issues down the road...

if nargin == 2
    thresh = 'global'; %Let the algorithm pick a threshold based on Otsu's method.
    l_min = 0; % threshold for minimum size of a domain. Set the default to zero.
    params.minimum_domain_size = 0;
elseif nargin == 3
    if isstruct(params)
        if isfield(params, 'threshold')
            if ischar(params.threshold)
                thresh = params.threshold; % the user can write 'adaptive' in the input field.
            else
                thresh = rescale(params.threshold, 'InputMin', min(pws(BW==1)), 'InputMax', max(pws(BW==1))); % Note that the threshold is in the units of the original image and then rescaled.
            end
        else
            thresh = 'global'; %Let the algorithm pick a threshold based on Otsu's method.
        end
        
        if isfield(params, 'minimum_domain_size')
            l_min = params.minimum_domain_size; % l_min should be a number.
        else
            l_min = 0; % set the default to 0. This might have the effect of over segmenting the domains.
            params.minimum_domain_size = 0;
        end
        
    else
        error('NNET:Arguments', 'Third input must be a structure array.')
    end
else
    error('NNET:Arguments', 'Wrong number of inputs. This function requires 2 or 3 inputs.')
end

%% Binarize the pws image to find domains

if ischar(thresh)
    if strcmp(thresh, 'global')
        thresh = graythresh(rescale(pws(BW==1))); % graythresh requires input image with values between 0 and 1.
        params.threshold = rescale(thresh, min(pws(BW==1)), max(pws(BW==1)), 'InputMin', 0, 'InputMax', 1); % rescale the threshold back to units of the input image
    end
end

pws_BW = double(BW) .* imbinarize(rescale(pws, 'InputMin', min(pws(BW==1)), 'InputMax', max(pws(BW==1))), thresh);
% Note that the entire PWS image is rescaled, but only based on the min and max values within the selected ROI.
% This should be compatible with adaptive thresholding and global thresholding.

%% find euclidean distances
dist = -bwdist(~pws_BW);
dist(~pws_BW) = -Inf;

%% Watershed segmentation
if l_min > 0
    dist = imhmin(dist, l_min); % the height threshold for suppressing shallow minima
end
L0 = double(watershed(dist));

for ii = 2:max(L0(:))
    L(:,:,ii-1) = L0==ii;
end

