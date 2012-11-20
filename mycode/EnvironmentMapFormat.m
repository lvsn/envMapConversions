classdef EnvironmentMapFormat
    % Represents the various kinds of environment map formats
    %   
    %   Currently:
    %   - LatLong (latitude-longitude)
    %   - Angular
    %   - SkyAngular 
    %       (just the top half of the hemisphere, looking towards zenith)
    %   - Cube
    %   - Octahedral
    %   - Sphere
    
    enumeration
        LatLong, Angular, SkyAngular, Cube, Octahedral, Sphere
    end
    
    methods (Static)
        
        function f = format(input)
            % Useful method to parse inputs for environment map formats
            %   
            %   f = format('formatName')
            %
            % Creates a format from a string. This will find the 'nearest'
            % string. Note: might do crazy things if the input string is
            % not "close" enough to the set of available options.
            %
            %   f = format(EnvironmentMapFormat.formatName)
            %
            % Simply returns the format directly. 
            %
            
            if isa(input, 'EnvironmentMapFormat')
                % we're good, just return the input
                f = input;
                
            elseif ischar(input)
                % find nearest string. 
                % 
                [members, names] = enumeration('EnvironmentMapFormat');
                
                dists = cellfun(@(s) strdist(lower(input), s), lower(names));
                [~,mind] = min(dists);
                
                f = members(mind);
            else
                error('Unsupported input for format conversion.');
            end
        end
    end
end

