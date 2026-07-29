function bnr_theory    = ComputeB_WM(An,Ln,D,r)
% Compute the Whittle-Matern ACF used in "Review of interferometric
% spectroscopy of scattered light for the quantification of
% subdiffractional structure of biomaterials" 

bnr_theory = An.*(r./Ln).^((D-3)/2).*besselk((D-3)/2,r./Ln);