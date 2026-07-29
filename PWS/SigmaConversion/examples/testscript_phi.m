clear all;
close all;
data = linspace(0.02, .3, 500);
phis = [0.32, .5, .6, .7, .8];
Nfs = [5e5, 1e6, 1.5e6];
nai = 0.52;
nac_input = 1.49;
thickIn = 2;

ri_media = 1.337;
riDef = S2D.RIDefinition(ri_media, nan, nan);
sysConfig = S2D.SystemConfiguration(riDef, nai, nac_input, 585, true, true);

for Nf = Nfs
    figure;
    hold on;
    for phi = phis
        sysConfig.ri_definition = S2D.RIDefinition.createFromGladstoneDale(ri_media, phi);
        [dOut, dCorrected, Nf_expected, lmax_corrected] = SigmaToD_AllInputs(data, sysConfig, Nf, thickIn);
        plot(dCorrected, data, 'DisplayName', ['Phi=', num2str(phi)]);
    end
    legend;
    ylabel('Sigma');
    xlabel('D');
    xlim([2 2.8]);
end
