classdef EnvironmentMap
    % Stores an environment map and its format.
    %   
    %   Very easy to convert to and from different formats. For example:
    %
    %   % Create environment map directly from image file, and specify
    %   % which format it's stored in
    %   e = EnvironmentMap(pathToFile, EnvironmentMapFormat.LatLong);
    %
    %   % Convert to the 'Angular' format
    %   e = e.convertTo(EnvironmentMapFormat.Angular);
    %
    %   % Resize 
    %   e = imresize(e, [e.nrows/2, NaN]);
    %
    %   % Display
    %   imshow(e);
    
    properties (GetAccess = public, SetAccess = private)
        data = [];
        format = 0;
    end
    
    properties (Dependent, SetAccess = private)
        % for convenience
        nrows;
        ncols;
        nbands;
    end
    
    methods
        function e = EnvironmentMap(varargin)
            % creates an EnvironmentMap object
            %
            %   e = EnvironmentMap(envmap, format)
            %
            % Creates an environment map from an image, and its associated
            % format.
            %
            %   e = EnvironmentMap(filename, format)
            % 
            % Creates an environment map from an image _file_, and its
            % associated format.
            
            assert(length(varargin)==2, 'Must have two inputs');
            assert(isa(varargin{2}, 'EnvironmentMapFormat'), ...
                'Second input must be of type ''EnvironmentMapFormat''');
            
            if ischar(varargin{1})
                % we're given the filename
                filename = varargin{1};
                [~,~,ext] = fileparts(filename);
                
                switch (ext)
                    case '.hdr'
                        % hdrread is from matlab's image processing toolbox
                        e.data = im2double(hdrread(filename));
                        
                    case '.exr'
                        % exrread requires the MatlabEXR package
                        % from K. Johnson
                        e.data = exrread(filename);
                        
                    otherwise
                        % we'll just rely on imread
                        e.data = im2double(imread(filename));
                        
                end
                
                
            else
                % we're given the image directly
                assert(isdouble(varargin{1}), ...
                    'First input must be an image in double format');
                assert(ndims(varargin{1}) <= 3, ...
                    'Can''t handle multi-dimensional environment maps');
                
                e.data = im2double(varargin{1});
            end
            
            e.format = lower(varargin{2});
        end
        
        % getters for the dependent properties
        function nrows = get.nrows(e)
            [nrows, ~, ~] = size(e.data);
        end
        
        function ncols = get.ncols(e)
            [~, ncols, ~] = size(e.data);
        end
        
        function nbands = get.nbands(e)
            [~, ~, nbands] = size(e.data);
        end
        
        % additional overloads
        function imshow(e)
            imshow(e.data);
        end
        
        function varargout = size(e, varargin)
            varargout = size(e.data, varargin{:}); 
        end
        
        function e = imresize(e, varargin)
            e.data = imresize(e.data, varargin{:});
        end
        
        function e = convertTo(e, tgtFormat)
            % Converts the environment map to the specified format
            %
            %   e = convertTo(e, tgtFormat)
            %
            
            assert(isa(tgtFormat, 'EnvironmentMapFormat'), ...
                'Format must be of type ''EnvironmentMapFormat''');
            
            tgtDim = e.nrows;
            
            switch (e.format)
                case EnvironmentMapFormat.LatLong
                    switch tgtFormat
                        case EnvironmentMapFormat.Angular
                            e.data = envmapLatLong2Angular(e.data, tgtDim);
                            
                        case EnvironmentMapFormat.Cube
                            e.data = envmapLatLong2Cube(e.data, tgtDim);
                            
                        case EnvironmentMapFormat.SkyAngular
                            e.data = envmapLatLong2SkyAngular(e.data, tgtDim);
                            
                        otherwise
                            error('Cannot convert from %s to %s\n', ...
                                e.format.char, tgtFormat.char);
                    end
                    
                case EnvironmentMapFormat.Angular
                    switch tgtFormat
                        case EnvironmentMapFormat.LatLong
                            e.data = envmapAngular2LatLong(e.data, tgtDim);
                            
%                         case EnvironmentMapFormat.Cube
%                             e.data = envmapAngular2Cube(e.data, tgtDim);
                            
%                         case EnvironmentMapFormat.SkyAngular
%                             e.data = envmapAngular2SkyAngular(e.data, tgtDim);

                        otherwise
                            error('Cannot convert from %s to %s\n', ...
                                e.format.char, tgtFormat.char);
                    end
                    
                case EnvironmentMapFormat.Cube
                    switch tgtFormat
%                         case EnvironmentMapFormat.LatLong
%                             e.data = envmapCube2LatLong(e.data, tgtDim);
                            
                        case EnvironmentMapFormat.Angular
                            e.data = envmapCube2Angular(e.data, tgtDim);
                            
%                         case EnvironmentMapFormat.SkyAngular
%                             e.data = envmapCube2SkyAngular(e.data, tgtDim);

                        otherwise
                            error('Cannot convert from %s to %s\n', ...
                                e.format.char, tgtFormat.char);
                    end
                    
                case EnvironmentMapFormat.SkyAngular
                    switch tgtFormat
                        case EnvironmentMapFormat.LatLong
                            e.data = envmapSkyAngular2LatLong(e.data, tgtDim);
                            
%                         case EnvironmentMapFormat.Angular
%                             e.data = envmapSkyAngular2Angular(e.data, tgtDim);
                            
%                         case EnvironmentMapFormat.Cube
%                             e.data = envmapSkyAngular2Cube(e.data, tgtDim);

                        otherwise
                            error('Cannot convert from %s to %s\n', ...
                                e.format.char, tgtFormat.char);

                    end
            end
            
            e.format = tgtFormat;
            
        end
    end
end

