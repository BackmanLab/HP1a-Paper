classdef SystemConfiguration
    %SYSTEMCONFIGURATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        cell_glass_interface = true; % If true then the fresnel coefficients will reflect that the cell_glass interface is the reference plane. Otherwise the cell/media interface will be used.
        ri_definition; % The colleciton of refractive indices used.
        immersion_ri;  % The refractive index of the immersion media for the objective.
        center_lambda = 585; % nm % center wavelength based on k
        na_c; % Collection NA
        na_i; % Illumination NA
    end
    
    methods
        function obj = SystemConfiguration(ri_def, na_i, na_c, center_lambda, oil_immersion, cell_glass_interface)
            %SYSTEMCONFIGURATION Construct an instance of this class
            %   oil_immersion: True/False. If true then immersion_ri will
            %   be set to the RI of glass found in `ri_def`. Otherwise the
            %   RI of media will be used. Air objectives are not currently
            %   supported.
            arguments
                ri_def S2D.RIDefinition
                na_i double
                na_c double
                center_lambda double  % The lambda associated with the center wavenumber of the system: lambda = 2 * pi / (k_max + k_min / 2)
                oil_immersion logical
                cell_glass_interface logical
            end
            obj.ri_definition = ri_def;
            obj.na_i = na_i;
            obj.na_c = na_c;
            obj.center_lambda = center_lambda;
            if oil_immersion
                obj.immersion_ri = ri_def.ri_glass;
            else
                obj.immersion_ri = ri_def.ri_media;
            end
            obj.cell_glass_interface = cell_glass_interface;
        end
        
    end
end

