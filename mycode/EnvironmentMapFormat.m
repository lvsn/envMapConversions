classdef EnvironmentMapFormat
    % Enumeration class to represent the various kinds of formats supported by
    % the EnvironmentMap class.
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
    %   - Omnidirectional (could be fisheye, catadioptric, etc. See
    %       documentation for the OCamCalib toolbox by Davide Scaramuzza at 
    %       https://sites.google.com/site/scarabotix/ocamcalib-toolbox
    
    enumeration
        LatLong, SkyLatLong, Angular, SkyAngular, Cube, Octahedral, SkyOctahedral, ...
            Sphere, SkySphere, Fisheye, Stereographic, Omnidirectional
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
                
                % make sure that format exists
                formatInd = strcmpi(names, input);
                
                if nnz(formatInd) ~= 1
                    fprintf('Available formats: \n');
                    fprintf('%s\n', names{:});
                    
                    error('EnvironmentMapFormat:badformat', ...
                        'Unsupported format %s.', input);
                end
                
                % return desired format
                f = members(formatInd);
                
            else
                error('EnvironmentMapFormat:badformat', ...
                    'Unsupported input for format conversion.');
            end
        end
    end
end

