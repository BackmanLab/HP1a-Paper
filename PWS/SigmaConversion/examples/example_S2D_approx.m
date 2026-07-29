close all
clear all
phi = 0.35;
Nf = 1e6;
ncCenterLambda = 2 * pi / ((2 * pi / 450 + 2 * pi / 700) / 2);
thickness = 2; % um
sigma = linspace(0.02, 0.25);
sMin = min(sigma);
sMax = max(sigma);

liveCellRI = S2D.RIDefinition.createFromGladstoneDale(1.337, phi);  

nuSys = S2D.SystemConfiguration(liveCellRI, 0.52, 1.49, 585, true, true);

fig = figure;
hold on;
polyVals = SigmaToD_polyApprox(nuSys, Nf, thickness, sMin, sMax, 5, true);
figure(fig);
plot(sigma, polyval(polyVals, sigma), 'DisplayName', '5th order');
polyVals2 = SigmaToD_polyApprox(nuSys, Nf, thickness, sMin, sMax, 2, true);
figure(fig);
plot(sigma, polyval(polyVals2, sigma), 'DisplayName', '2nd order');
legend;
xlabel('Sigma');
ylabel('D');

tic;
sigma = randi(200, 10000, 10000) / 1000.0;
dOut = polyval(polyVals, sigma);
['Approximated a hundred million sigma values in ', num2str(toc), ' seconds.']

