function d_corrected = acfd(d, lmin, lmax)
% This function is used to convert uncorrected D (Db) to D_corrected.
% Corresponds to Eqn. 3 from the paper.
% d: may be scalar or a 2d array.
% lmin: Must be a scalar
% lmax: The shape must match that of d.

if ~isequal(size(lmax), size(d))
   error('lmax and d must be scalar or 2d arrays of the same size') 
end

delta = 0.1;
mid = 60; % Aya: Unfortunately, this fit location is a very sensitive component.   Technically, you want it somewhere between lmin and lmax (which are on the order of 1 and 200), and more specifically, you want the mid point to be in the center of the linear region of the ACF.  So given all the parameters, what does your Bn(r) look like (on a loglog scale).  Where is the linear region (eg. where r^(3-D) can be appropriately fit)?  What is the midpoint of that linear regime?  
point = (lmax+lmin)/mid;

[n, m] = size(d);
d_corrected = zeros(n, m);
for i=1:n
    for j=1:m
        d_corrected(i, j) = 3 + (log(S2D.ACF.ComputeB_SD(d(i, j), lmin, lmax(i, j), point(i, j) + delta)) - ...
            log(S2D.ACF.ComputeB_SD(d(i, j), lmin, lmax(i, j), point(i, j))))./ ...
            (log(point(i, j) + delta) - log(point(i, j)));
    end
end

