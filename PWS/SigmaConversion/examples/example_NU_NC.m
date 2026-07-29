phi = 0.35;
Nf = 1e6;
ncCenterLambda = 2 * pi / ((2 * pi / 450 + 2 * pi / 700) / 2);
thickness = 2; % um
sigma = linspace(0.02, 0.3);

liveCellRI = S2D.RIDefinition.createFromGladstoneDale(1.337, phi);

nuSys = S2D.SystemConfiguration(liveCellRI, 0.52, 1.49, 585, true, true);
ncSys = S2D.SystemConfiguration(liveCellRI, 0.6, 0.8, ncCenterLambda, false, true);

figure;
hold on
grid on;
[dOut,dCorrected, Nf_expected,lmax_corrected] = SigmaToD_AllInputs(sigma, nuSys, Nf, thickness);
plot(dCorrected, sigma, 'DisplayName', 'NU NA_i=0.52, NA_c=1.49');
[dOut,dCorrected, Nf_expected,lmax_corrected] = SigmaToD_AllInputs(sigma, ncSys, Nf, thickness);
plot(dCorrected, sigma, 'DisplayName', 'NanoCytomics NA_i=0.6, NA_c=0.8');
legend('Location', 'NorthWest');
ylabel('Sigma');
xlabel('D');