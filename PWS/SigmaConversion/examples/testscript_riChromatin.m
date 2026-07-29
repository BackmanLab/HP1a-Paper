clear all
close all
figure
phi = linspace(.3, 1, 100);

ri_media = 1.363; % Ethanol  1.332 Cell Media
lmin = 1; % Minimum fractal length that is possible % 
phi_mc = 0.05; % 0.05; % Phi is CVC - mobile crowders - 5% of remaning volume that is unoccupied
rho_protein = 1.35; %1.35;  % 1.35 in literautre - Data from Sciences paper 2017 on ChromEM - Density of protein
rho_chromatin = 0.56; % Data from Sciences paper 2017 on ChromEM - Calculated by VB
% effective refractive index of each component, RI_increment (alpha),
% RIcell, RIglass, RImedia
lambda  = 585; % nm % Change to NC peak wavelength?? - center wavelength based on k
riinc   = 0.1799; %0.1899*0 + (0.1899*1/2 + 0.17*1/2 + 0.185*0); %?? - Gladstone equation's alpha
densityProtein = (rho_chromatin*phi + rho_protein*phi_mc*(1-phi));
ric     = ri_media + riinc*densityProtein; %RI chromatin    
sigma_n = sqrt(phi.*(1 - phi))* (riinc*rho_chromatin*(1 - phi_mc) - riinc * rho_protein * phi_mc); % Eqn10 from the paper.    

plot(phi, ric)
xlabel('CVC');
ylabel('RI chromatin')
title('Gladstone-Dale Equation (Ethanol)')

figure
plot(phi, sigma_n)
xlabel('CVC');
ylabel('sigma_n')
title('Gladstone-Dale Equation (Ethanol)')

figure
plot(phi, densityProtein)
xlabel('CVC');
ylabel('Density of Protein')
title('Gladstone-Dale Equation (Ethanol)')