function [polyFit, mse] = SigmaToD_polyApprox(system_config, Nf, thickIn, sigmaMin, sigmaMax, polynomialOrder, debug)
%SigmaToD_polynomialApprox This function uses the exact solution to generate a lookup table
%and then fits a polynomial to it.
%
% The return value `polyFit` is the polynomial that has been fit tot he exact solution.
% It can be passed to the `polyval` function to quickly approximate D
% values for a large number of sigma values.
arguments
    system_config S2D.SystemConfiguration
    Nf double
    thickIn double
    sigmaMin double  % The minimum value of the sigma range to fit to.
    sigmaMax double  % The maximum value of the sigma range to fit to.
    polynomialOrder int32; % Use `debug` to help select the best polynomial order for the fit.
    debug logical = 0  % If true then a figure will be opened comparing the exact and approximate solutions. fitError ouptut argument will be the mean squared error (MSE) between the two methods.
end
mse = nan; % This is only evaluated if debug is true.
sigmaLUT = linspace(sigmaMin, sigmaMax, 1000);
[Db, D_LUT] = SigmaToD_AllInputs(sigmaLUT, system_config, Nf, thickIn);
polyFit = polyfit(sigmaLUT, D_LUT, polynomialOrder);
if debug
    D_approx = polyval(polyFit, sigmaLUT);  % D values approximated from our recent fit.
    mse = immse(D_LUT, D_approx);  % The error between exact and approximated.
    figure;
    hold on;
    plot(sigmaLUT, D_LUT, 'DisplayName', 'Exact');
    plot(sigmaLUT, D_approx, 'DisplayName', 'Approx.');
    legend;
    xlabel('Sigma');
    ylabel('D');
    title(['MSE: ', num2str(mse)]);
end
end

