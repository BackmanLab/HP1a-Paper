function [x] = numericalInversion(func, xVals, y)
%NUMERICALINVERSION For a function y = f(x), use interpolation and a lookup
%table to approximate x = g(y) where g is the inverse of f.
%   func: A `function_handle` that accepts an array of x values as the
%       first argument and returns a corresponding array of y values.
%   xVals: The array of x values to pass to `func` for the purpose of generating a lookup table
%   y: The input to the inverted function for which you want to get an x
%       value.
% Author: Nick Anthony
yVals = func(xVals); % Generate y values corresponding to each x value.
x = interp1(yVals, xVals, y, 'pchip', 'extrap');  % The original code used linear interpolation with no exptrapolation. Probably would be faster but could also have problems.
if (x < min(xVals)) || (x > max(xVals))
    warning(['Numerical Inversion resorting to extrapolation. x: ', num2str(x), ' y: ', num2str(y)]);
end

end

