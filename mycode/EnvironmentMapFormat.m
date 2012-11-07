classdef EnvironmentMapFormat
    % Represents the various kinds of environment map formats
    %   
    %   Currently:
    %   - LatLong (latitude-longitude)
    %   - Angular
    %   - Cube
    %   - SkyAngular 
    %       (just the top half of the hemisphere, looking towards zenith)
    
    enumeration
        LatLong, Angular, Cube, SkyAngular
    end
    
    methods (Static)
        
        function f = str2format(str)
            % Useful method to convert from string to class format
            %   
            %   f = str2format(str)
            %
            switch lower(str)
                case 'latlong'
                    f = EnvironmentMapFormat.LatLong;
                    
                case 'angular'
                    f = EnvironmentMapFormat.Angular;
                    
                case 'cube'
                    f = EnvironmentMapFormat.Cube;
                    
                case 'skyangular'
                    f = EnvironmentMapFormat.SkyAngular;
                    
                otherwise
                    error('Unknown string: %s\n', str);
            end
        end
    end
end

