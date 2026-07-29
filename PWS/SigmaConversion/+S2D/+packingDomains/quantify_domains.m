%% function [stats, Lmetrics] = quantify_domains(D, L, pws_BW, BW)
%
% DESCRIPTION
%   This function uses the outputs of find_domains.m and calculates useful
%   statistics from the domains found in PWS microscopy images.
%
% INPUT ARGUMENTS
%   D:
%       Chromatin packing scaling which can be calculated using SigmaToD.m.
%       Sigma images may also be used, however note that outputs referring
%       to "mean D" will be in Sigma rather than D.
%   L:
%       A 3D stack of segmented domain maps. This is an output from
%       find_domains.m.
%   pws_BW:
%       A binary mask showing regions which were thresholded and found to
%       be domains. This is an output of find_domains.m
%   BW:
%       ROI map selecting the nucleus of a cell. All analysis will only be 
%       performed inside the ROI.
%
% OUTPUT ARGUMENTS
%   stats:
%       Structure array with fields defining various statistics of domains.
%       The fields are:
%
%       "mean_D_domains" is the average D value from all domains. Each
%       domain is given equal weight regardless of size.
%
%       "mean_D_domain_space" is the spatially averaged D valuse of the
%       space of the nucleus made of domains.
%
%       "mean_D_nondomain_space" is the spatially averaged D of all space
%       within the nucleus not part of domains.
%
%       "domain_projection_fraction" is the fraction of space in the
%       nucleus made of domains. 
%
%       "number_domains" is the total number of domains in the nucleus.
%       Note that if the resolution of the PWS system is not adequate to
%       resolve individual domains, this number may be biased too small.
%
%       "mean_domain_size" is the average size (measured in pixels) of a
%       domains. Note that if the resolution of the PWS system 
%       is not adequate to resolve individual domains, then this number may
%       be biased too large.
%
%       "Lmetrics_labels" labels for each column in "Lmetrics"
%   Lmetrics:
%       Statistics on individual domains. As defined in
%       stats.Lmetrics_labels, the first column defines the average D value
%       of each domain, while the secomd column defines the size of that
%       domain measured in pixels. 
%
% EXAMPLES
%   [stats, Lmetrics] = quantify_domains(D, L, pws_BW, BW)
%   [stats, Lmetrics] = quantify_domains(SigmaToD(cubeRms), L, pws_BW, n2)
%   stats = quantify_domains(SigmaToD(cubeRms), L, pws_BW, BW)
%
% REFERENCES
%   https://doi.org/10.1101/2020.01.26.920363 
%
% Author: Adam Eshein (aeshein@u.northwestern.edu) 5.9.2020      
function [stats, Lmetrics] = quantify_domains(D, L, pws_BW, BW)

Lmetrics(:,1) = squeeze(sum(sum(double(L) .* double(repmat(D, [1 1 size(L,3)])) .* repmat(double(BW), [1 1 size(L,3)]))) ./ sum(sum(L)));
Lmetrics(:,2) = squeeze(sum(sum(L)));

stats.Lmetrics_labels = [{'Average D'}, {'Size in pixels'}];

stats.mean_D_domains = mean(Lmetrics(:,1));
stats.mean_D_nondomain_space = sum(sum(double(D) .* double(~pws_BW) .* double(BW))) / sum(sum(double(~pws_BW) .* double(BW)));
stats.mean_D_domain_space = sum(sum(double(D) .* double(pws_BW) .* double(BW))) / sum(sum(double(pws_BW) .* double(BW)));
stats.domain_projection_fraction = sum(sum(double(pws_BW) .* double(BW))) / sum(sum(BW));
stats.number_domains = size(Lmetrics,1);
stats.mean_domain_size = mean(Lmetrics(:,2));