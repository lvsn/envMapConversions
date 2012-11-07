classdef EnvironmentMapFormat
    %UNTITLED13 Summary of this class goes here
    %   Detailed explanation goes here
    
    enumeration
        LatLong, Angular, Cube, SkyAngular
    end
    
    methods (Static)
        % useful method to convert from string to class format
        function f = str2format(str)
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

