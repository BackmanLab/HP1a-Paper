%%% This example plots how changes to the illumination numerical aperture
%%% will affect the relationship between Sigma and D.

close all

phi = 0.35;  % The CVC of the cells in question. Measured by ChromSTEM
Nf = 1e6;  % The genomic length of a packing domain. Measured by Hi-C
thickness = 2; % in microns. The sample thickness. This only has an effect if it is less thane the depth of focus. In the case of a 0.52 NA_i the dof is 1.07 microns
sigma = linspace(0.02, 0.3);  % We are just generating a range of sigma values to use.

liveCellRI = S2D.RIDefinition.createFromGladstoneDale(1.337, phi);  % The Refractive index of sigma_n of chromatin can be calculated based on the CVC and the immersion RI. In this case the immersion is PBS (RI=1.337)

figure;
hold on;
for NA_i = 0.1:.1:1  % Iterate over a range of illumination numerical apertures.
    nuSys = S2D.SystemConfiguration(liveCellRI, NA_i, 1.49, 585, true, true);  % System config for the Northwestern LCPWS2, STORM, LCPWS3 systems.
    [dOut, dCorrected, Nf_expected, lmax_corrected] = SigmaToD_AllInputs(sigma, nuSys, Nf, thickness); % Compute the D values.
    plot(sigma, dCorrected, 'DisplayName', num2str(NA_i)); % Plot sigma vs D.
end
legend;
xlabel('Sigma');
ylabel('D');

figure;
hold on;
for NA_c = 0.6:.1:1.5  % Iterate over a range of NA collection values.
    nuSys = S2D.SystemConfiguration(liveCellRI, 0.52, NA_c, 585, true, true);
    [dOut,dCorrected, Nf_expected,lmax_corrected] = SigmaToD_AllInputs(sigma, nuSys, Nf, thickness);
    plot(sigma, dCorrected, 'DisplayName', num2str(NA_c));
end
legend;
xlabel('Sigma');
ylabel('D');

figure;
hold on;
for oilImmersion = [false, true]  % Iterate over using oil immersion or media immersion.
    nuSys = S2D.SystemConfiguration(liveCellRI, 0.6, .8, 585, oilImmersion, true);
    [dOut,dCorrected, Nf_expected,lmax_corrected] = SigmaToD_AllInputs(sigma, nuSys, Nf, thickness);
    plot(sigma, dCorrected, 'DisplayName', num2str(oilImmersion));
end
legend;
xlabel('Sigma');
ylabel('D');

figure;
hold on;
for glassRef = [false, true]  % Iterate over using the bottom or top surface of the cell as the "reference" reflection.
    nuSys = S2D.SystemConfiguration(liveCellRI, 0.6, .8, 585, true, glassRef);
    [dOut,dCorrected, Nf_expected,lmax_corrected] = SigmaToD_AllInputs(sigma, nuSys, Nf, thickness);
    plot(sigma, dCorrected, 'DisplayName', num2str(glassRef));
end
legend;
xlabel('Sigma');
ylabel('D');

figure;
hold on;
for cellRi = 1.37:.02:1.5  % Iterate over various Chromatin RI values. Assume that sigma_n stays as 0.04
    liveCellRI = S2D.RIDefinition(1.361, cellRi, 0.04);  
    ncSys = S2D.SystemConfiguration(liveCellRI, 0.6, .8, 585, true, true);  % Standard NanoCytomics config.
    [dOut,dCorrected, Nf_expected,lmax_corrected] = SigmaToD_AllInputs(sigma, ncSys, Nf, thickness);
    plot(sigma, dCorrected, 'DisplayName', num2str(cellRi));
end
legend;
xlabel('Sigma');
ylabel('D');
title(['Sigma Vs. D as a function of RI_{chromatin}']);