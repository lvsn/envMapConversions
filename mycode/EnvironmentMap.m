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
                assert(isnumeric(varargin{1}), ...
                    'First input must be an image');
                assert(ndims(varargin{1}) <= 3, ...
                    'Can''t handle multi-dimensional environment maps');
                
                e.data = im2double(varargin{1});
            end
            
            e.format = EnvironmentMapFormat.format(varargin{2});
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
            varargout{:} = size(e.data, varargin{:}); 
        end
        
        function e = imresize(e, varargin)
            e.data = imresize(e.data, varargin{:});
        end
        
        function e = imfilter(e, varargin)
            e.data = imfilter(e.data, varargin{:});
        end
        
        function display(e)
            fprintf('EnvironmentMap, [%dx%dx%d], %s\n', ...
                e.nrows, e.ncols, e.nbands, e.format.char);
        end
        
        function e = intensity(e)
            % Returns intensity-version of the environment map
            if e.nbands == 1
                warning('EnvironmentMap:intensity', ...
                    'Environment map already is intensity-only');
            else
                e.data = rgb2gray(e.data);
            end
        end
            
        
        function e = convertTo(e, tgtFormat, tgtDim)
            % Converts the environment map to the specified format
            %
            %   e = convertTo(e, tgtFormat, <tgtDim>)
            %
            
            tgtFormat = EnvironmentMapFormat.format(tgtFormat);
            
            if ~exist('tgtDim', 'var')
                tgtDim = e.nrows;
            end
            
            if tgtFormat == e.format
                % we already have the right format! 
                % Let's make sure we have the right size
                e = imresize(e, [tgtDim NaN]);
                return;
            end
            
            % Get the world coordinates
            [dx, dy, dz, valid] = EnvironmentMap.worldCoordinatesStatic(tgtFormat, tgtDim);
            
            % Get the environment map representation
            switch (e.format)
                case EnvironmentMapFormat.LatLong
                    e.data = envmapWorld2LatLong(e.data, dx, dy, dz);
                    
                case EnvironmentMapFormat.Angular
                    e.data = envmapWorld2Angular(e.data, dx, dy, dz);
                    
                case EnvironmentMapFormat.Cube
                    e.data = envmapWorld2Cube(e.data, dx, dy, dz);
                    
                case EnvironmentMapFormat.SkyAngular
                    e.data = envmapWorld2SkyAngular(e.data, dx, dy, dz);
                    
                case EnvironmentMapFormat.Octahedral
                    e.data = envmapWorld2Octahedral(e.data, dx, dy, dz);
                    
                case EnvironmentMapFormat.Sphere
                    e.data = envmapWorld2Sphere(e.data, dx, dy, dz);
                    
                otherwise
                    error('Conversion from %s to %s unsupported!', ...
                        e.format.char, tgtFormat.char);
            end
            
            e.data(~valid(:,:,ones(1,e.nbands))) = 0;
            e.format = tgtFormat;            
        end
        
        function [x, y, z, valid] = worldCoordinates(e)
            % Returns the [x,y,z] world coordinates
            [x, y, z, valid] = EnvironmentMap.worldCoordinatesStatic(e.format, e.nrows);
        end

    end
    
    methods (Static)
                
        function [x, y, z, valid] = worldCoordinatesStatic(format, dims)
            % Returns the [x,y,z] world coordinates
            
            assert(isa(format, 'EnvironmentMapFormat'), ...
                'Input format is expected to be of type ''EnvironmentMapFormat''.');

            switch(format)
                case EnvironmentMapFormat.LatLong
                    [x, y, z, valid] = envmapLatLong2World(dims);
                                        
                case EnvironmentMapFormat.Angular
                    [x, y, z, valid] = envmapAngular2World(dims);
                    
                case EnvironmentMapFormat.Cube
                    % we need to divide the dimensions by 4 since they
                    % assume we're talking about only one side of the cube
                    % and not the cube "image" 
                    [x, y, z, valid] = envmapCube2World(dims./4);
                    
                case EnvironmentMapFormat.SkyAngular
                    [x, y, z, valid] = envmapSkyAngular2World(dims);
                    
                case EnvironmentMapFormat.Octahedral
                    [x, y, z, valid] = envmapOctahedral2World(dims);
                    
                case EnvironmentMapFormat.Sphere
                    [x, y, z, valid] = envmapSphere2World(dims);

                otherwise
                    error('IlluminationModel:getWorldCoordinates', ...
                        'Unsupported format: %s', format.char);
            end
        end
    end
end

