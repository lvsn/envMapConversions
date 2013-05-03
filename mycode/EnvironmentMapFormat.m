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
    %   - SkyOctahedral
    %   - Sphere
    %   - SkySphere
    %   - Fisheye (r = f*theta)
    %   - Stereographic (specific type of fisheye, r = 2r tan(theta/2))
    
    enumeration
        LatLong, Angular, SkyAngular, Cube, Octahedral, SkyOctahedral, ...
            Sphere, SkySphere, Fisheye, Stereographic
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
                [~, names] = enumeration('EnvironmentMapFormat');
                
                % make sure that format exists
                if nnz(strcmpi(names, input)) ~= 1
                    fprintf('Available formats: \n');
                    fprintf('%s\n', names{:});
                    
                    error('envmapFormat:badformat', ...
                        'Unsupported format %s.', input);
                end
                
            else
                error('envmapFormat:badformat', ...
                    'Unsupported input for format conversion.');
            end
        end
    end
end

