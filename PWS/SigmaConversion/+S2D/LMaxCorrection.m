function [lmax_corrected] = LMaxCorrection(Db, lmin, Nf)
%LMAXCORRECTION Performs Aya's correction fo lMax rather than the previous
%approach of correcting for Nf. Corresponds to Eqn. 2 from Aya's paper.
% Db: Can be a scalar or a 2d array of uncorrected D values (Db)
% lmin: The minimum size (in nanometers) over which we expect the structure
% to be fractal. Usually set to 1nm.
% Nf: The genomic length of the packing domain.

lmax_r = 10:30:25000;   % An array of possible lmax values
[n, m] = size(Db);
lmax_corrected = zeros(n, m);
for i=1:n
   for j=1:m
       d = Db(i, j);
       eqn2 = @(lmax_r) 6.*(d-3)./d .* (1 - (lmin./lmax_r).^d)./((lmin./lmax_r).^3 - (lmin./lmax_r).^d);  % We invert this function to get lmax as a function of Nf and D.
       lmax_corrected(i, j) = S2D.utility.numericalInversion(eqn2, lmax_r, Nf);
   end
end

end

