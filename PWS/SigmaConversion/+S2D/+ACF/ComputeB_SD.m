function bnr    = ComputeB_SD(d,lmin,lmax,r)
% Compute the ACF used in the "Size Distribution" model. See
% "Characterizing chromatin packing scaling in whole nuclei using
% interferometric microscopy".

% Evaluates our model of the ACF of chromatin for the input parameters.
% d: The uncorrected fractal dimension (Db). Must be scalar.
% lmin: The minimum size in nanometers that we expect this ACF to be valid
% for. Must be scalar.
% lmax: The maximum size in nanometers that we expect this ACF to be
% fractal for. Must be scalar.
% r: The size or sizes to evaluate the function at. May be a 1d array or scalar.
bnr = (3-d).* r.^(d-3) ./ (lmin.^(d-3) - lmax.^(d-3)) .*  (S2D.utility.gamma_incomplete(r./lmax,3-d) - S2D.utility.gamma_incomplete(r./lmin, 3-d));
end