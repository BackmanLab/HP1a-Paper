function [dOut, dCorrected, Nf_expected, lmax_corrected] = SigmaToD_AllInputs(sigmaIn, system_config, Nf, thickIn)
% This is the most complete calculation of the Sigma to D conversion. Check
% Aya's paper "Characterizing chromatin packing scaling in whole nuclei
% using interferometric microscopy" for more information.
%
% Parameters:
%   sigmaIn: The sigma values you want to convert.
%   system_config: The SystemConfiguration object to use.
%   Nf: The genomic length of a packing domain. Often values of 5e5 or 1e6
%   are used.
%   thickIn: The expected thickness of the sample. If this is greater than
%   the calculated DOF of the microscope then it has no effect. The DOF of
%   LCPWS2 is 1.1um.
%
% Returns:
%   dOut: This is analogous to `D_b` in Aya's paper. The `model parameter`.
%   dCorrected: This is the true `D`. This is usually what we care about.
%   Nf_expected: The genomic length we expect based on D and the calculated lMax.
%   lmax_corrected: Calculation of LMax from Nf and Db based on eqn. 2.
arguments
    sigmaIn
    system_config S2D.SystemConfiguration
    Nf double
    thickIn double
end

lmin = 1; % Minimum fractal length that is possible % 
ri_glass = system_config.ri_definition.ri_glass;
ri_chromatin = system_config.ri_definition.ri_chromatin;
ri_media = system_config.ri_definition.ri_media;
ri_immersion = system_config.immersion_ri;

% Fresnel Explanation:
%
% n0 = RI of material the light is traveling from. n1 = RI of the
% material the light is entering.
%
% Forward amplitude transmission coeff:
% t_f = (2 * n0) / (n0 + n1)
%
% Backward amplitude transmission coeff t_b is the same as forward but with RI
% values swapped.
%
% Reflectance amplitude transmission coeff:
% r = (n0 - n1) / (n0 + n1)
%
% Power transmission is defined as:
% T_f = n1 / n0 * t_f^2 and T_b = n0 / n1 * t_b^2
%
% Power reflection is:
% R = r^2 = ((n0 - n1) / (n0 + n1))^2
%
% The "Fresnel Coefficient" (equation 8 from the paper, confusingly named
% as `R` but named `Fcoeff` here) is:
% Fcoeff = t_f * t_b * r / R_reference
%
% Combining all the above equations we get:
% Fcoeff = (2 * n0 / (n0 + n1)) * (2 * n1 / (n0 + n1)) * ((n0 - n1) / (n1 + n0)) / R_reference
%
% Terms cancel out to produce the following simpler expression:
% Fcoeff = 4 * n0 * n1 * (n0 - n1) / (n0 + n1)^3 / R_reference
R_reference = ((ri_glass - ri_media) / (ri_glass + ri_media))^2;
if system_config.cell_glass_interface
    n0 = system_config.ri_definition.ri_glass;
    n1 = system_config.ri_definition.ri_chromatin;
else
    n0 = system_config.ri_definition.ri_media;
    n1 = system_config.ri_definition.ri_chromatin;
end
Fcoeff = 4 * n0 * n1 * (n0 - n1) / (n0 + n1)^3 / R_reference;  % Eqn 8 from paper. Combination of fresnel coefficients normalized by reflectance of reference.


% defining the cell effective thickness
dof = system_config.center_lambda * ri_immersion / (2 * system_config.na_i^2) / 1e3;  % Depth of focus in microns. Berek Formula.
thick = min(dof, thickIn)*1e3;

% TODO the paper says "k is the center wavelength in vacuum". First, it should say
% wavenumber, not wavelength. More importantly. In the code it is defined as the 
% wavenumber in chromatin. I decided to go with the paper.
k = 2 * pi / system_config.center_lambda;

% Closed form solution for sigma^2 for an exponential ACF. Eqn. 7
s2_to_int = @(lc, L, lambda, ric, na, lmin, lmax, d) ...
    Fcoeff ^ 2 * system_config.ri_definition.sigma_n^2 * (...
    (2/pi * (L .* (k^4 .* lc.^3) * na^2) ...  
    ./ (1 + (k.*lc).^2*(4 + na^2)) ...
    ./ (1 + (k.*lc).^2 * 4)) ...
     + 1/4 * (1 - 1./sqrt(1 + (k.*lc * na).^2)));
       
% Probability distribution function. Eqn. 5 from the paper.
Pr = @(lc, lmin,lmax,D) lc.^(D-4) .* (D-3)./(lmax.^(D-3) - lmin.^(D-3)); 


dfs = linspace(2.1,2.9,20); % Range of possible Df values. Same as Db in the paper
sigmaLUT = zeros(size(dfs));   % This will be a lookup table to sigma.
lmax_r = 10:20:10000;  % Array of lmax values. In units of nm.

for idf = 1:length(dfs) % range of Db's
    df = dfs(idf); % Db from the paper

    eqn2 = @(lmax_r) 6.*(df-3)./df .* (1 - (lmin./lmax_r).^df)./((lmin./lmax_r).^3 - (lmin./lmax_r).^df); %Eqn 2. from the Paper. Array of Nfs for corresponding values of lmax_r % This is equation 2 in the paper.
    lmax = S2D.utility.numericalInversion(eqn2, lmax_r, Nf);  % to get lmax corresponding to Nf. Nf = (lmax/lmin)^D cannot be inverted since we know Db and not D
    sigmaLUT(idf) = sqrt(integral(@(lc) (Pr(lc,lmin,lmax,df).* s2_to_int(lc, thick, system_config.center_lambda, ri_chromatin, system_config.na_c, lmin, lmax, df)), lmin, lmax));
end

dVSigma = polyfit(sigmaLUT,dfs, 3); % Fit a cubic to the sigma - Db relationship
dOut = polyval(dVSigma,sigmaIn);  % get Db from our input sigma value.

%% Correcting for lmax instead of nf
lmax_corrected = S2D.LMaxCorrection(dOut, lmin, Nf);  % Calculate LMax from Nf and Db based on eqn. 2.

dCorrected  = S2D.acfd(dOut, 1, lmax_corrected); % Calculate true D from Db
Nf_expected = (lmax_corrected/lmin).^dCorrected;  % The genomic length we expect based on D and the calculated lMax.
% %}