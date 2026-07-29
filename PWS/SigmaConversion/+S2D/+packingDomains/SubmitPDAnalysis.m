%% About this code
%Author: ARS 3/17/21 
%Needs: load_pws.m, find_domains.m, quantify_domains.m
%Needs: SigmaToD.m, acf_1.m, acfd.m, SigmaToD_coefs.mat
%% Run Code
mean_D_domains={};
number_domains={};
mean_domain_size={};

directory='/Users/anneshim/Desktop/WT/';
cells={'1','2'};

for cellnum=1:length(cells)
    %f = fullfile(directory,'Cell',cells{cellnum})
    f = [directory 'Cell' cells{cellnum}]
    [cubeRms BW] = S2D.packingDomains.load_pws(f, 'p0', 'nuc'));
    sz=size(BW);
    indexes=sz(3);
    
    for index=1:indexes
        n1=BW(:,:,index);
        [L, pws_BW, params] = S2D.packingDomains.find_domains(cubeRms, n1);
        [stats, Lmetrics] = S2D.packingDomains.quantify_domains(S2D.SigmaToD(cubeRms,2.43, 0.52, 0.036,[]), L, pws_BW, n1);
        mean_D_domains=[mean_D_domains,stats.mean_D_domains];
        number_domains=[number_domains,stats.number_domains];
        mean_domain_size=[mean_domain_size,stats.mean_domain_size];
    end
end

disp(vertcat(mean_D_domains{:}))
disp(vertcat(number_domains{:}))
disp(vertcat(mean_domain_size{:}))