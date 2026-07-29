classdef RIDefinition
    %RIDEFINITION This class is used to bundle together RI values used in
    %the sigma to D conversion.
    
    properties
        ri_glass = 1.517;
        ri_media;
        ri_chromatin;  % average RI in chromatin
        sigma_n; % sigma_n is the standard deviation of RI fluctuations in chromatin
    end
    
    methods
        function obj = RIDefinition(ri_media, ri_chromatin, sigma_n)
            %RIDEFINITION Construct an instance of this class by passing
            %the RI values in directly. More commonly you will use
            %`createFromGladstoneDale` to create an instance from `phi`.
            obj.ri_chromatin = ri_chromatin;
            obj.ri_media = ri_media;
            obj.sigma_n = sigma_n;
        end
    end
    
    methods (Static)
       function obj = createFromGladstoneDale(ri_media, phi)
            %Create a new instance of this class where `ri_chromatin` and `sigma_n` based on `phi` using the gladstone-
            % dale equation. WARNING: MATLAB doesn't modify the existing
            % instance. It returns a modified copy.
            %   phi: The CVC of the cell nucleus.
            phi_mc = 0.05; % 0.05; % Phi is CVC - mobile crowders - 5% of remaning volume that is unoccupied
            rho_protein = 1.35; % 1.35 in literature - Data from Sciences paper 2017 on ChromEM - Density of protein
            rho_chromatin = 0.56; % Data from Sciences paper 2017 on ChromEM - Calculated by VB. Density of chromatin.
            riinc   = 0.1799; %0.1899*0 + (0.1899*1/2 + 0.17*1/2 + 0.185*0); %?? - Gladstone equation's alpha
            ri_chromatin = ri_media + riinc*(rho_chromatin*phi + rho_protein*phi_mc*(1-phi)); %RI chromatin   
            sigma_n = sqrt(phi*(1 - phi))* (riinc*rho_chromatin*(1 - phi_mc) - riinc * rho_protein * phi_mc); % Eqn10 from the paper.    
            obj = S2D.RIDefinition(ri_media, ri_chromatin, sigma_n);
        end
    end
end

