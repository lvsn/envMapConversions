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
        focalLength = [];
        hFov = [];
        vFov = [];
        
        fisheyeMappingFcn = [];
        fisheyeInverseFcn = [];
    end
    
    properties (Dependent, SetAccess = private)
        % for convenience
        nrows;
        ncols;
        nbands;
    end
    
    methods
        function e = EnvironmentMap(input, format, varargin)
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
            %
            %   e = EnvironmentMap(dims, format)
            %
            % Creates an all-zeros environment map with vertical
            % dimensions 'dims' (horizontal dimensions automatically
            % determined based on format).
            %
            %   e = EnvironmentMap(..., <par1, val1>, ...)
            %
            % Creates environment map with additional name/value pairs:
            % 
            %   - 'hfov': horizontal field of view (degrees)
            %   - 'vfov': vertical field of view (degrees)
            %   - 'dfov': diagonal field of view (degrees)
            
            
            hfov = [];
            vfov = [];
            dfov = [];
            
            focal = [];
            
            parseVarargin(varargin{:});
                        
            if ischar(input)
                % we're given the filename
                filename = input;
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
                
            elseif isscalar(input)
                % we're given the vertical dimension only --> create
                % all-zeros
                switch EnvironmentMapFormat.format(format)
                    case EnvironmentMapFormat.LatLong
                        e.data = zeros(input, input*2);
                        
                    case EnvironmentMapFormat.Cube
                        e.data = zeros(4*dims, 3*dims);
                        
                    otherwise
                        % square data
                        e.data = zeros(input, input);
                end
                
            else
                % we're given the image directly
                assert(isnumeric(input), ...
                    'First input must be an image');
                assert(ndims(input) <= 3, ...
                    'Can''t handle multi-dimensional environment maps');
                
                e.data = im2double(input);
            end
            
            e.format = EnvironmentMapFormat.format(format);
            
            switch e.format 
                case EnvironmentMapFormat.Stereographic
                    assert(~isempty(hfov) | ~isempty(vfov) | ~isempty(dfov), ...
                        'Stereographic format requires field of view.');
                    % re-warp the environment map such that the horizontal fov
                    % maps the input (adjust the vertical fov accordingly assuming square pixels)
                    
                    if ~isempty(hfov)
                        %                     e.hFov = hfov*pi/180;
                        %                     e.vFov = e.hFov*size(e.data,1)/size(e.data,2);
                        hfov = hfov*pi/180;
                        theta = hfov/2;
                        r = size(e.data,2)/2;
                        
                    elseif ~isempty(vfov)
                        vfov = vfov*pi/180;
                        theta = vfov/2;
                        r = size(e.data,1)/2;
                        
                    elseif ~isempty(dfov)
                        r = sqrt(size(e.data,1)^2+size(e.data,2)^2)/2;
                        dfov = dfov*pi/180;
                        theta = dfov/2;
                        
                    end
                    e.focalLength = r/(2*tan(theta/2));
                    
                    % compute the vertical/horizontal fovs
                    e.hFov = 2*atan2(size(e.data,2)/2, 2*e.focalLength);
                    e.vFov = 2*atan2(size(e.data,1)/2, 2*e.focalLength);
                    
                    % the mapping function is r = 2*f*tan(theta/2)
                    e.fisheyeMappingFcn = @(theta) 2.*e.focalLength.*tan(theta/2);
                    % the inverse function is theta = 2*atan(r/2f);
                    e.fisheyeInverseFcn = @(r) 2*atan2(r, 2*e.focalLength);
                    
                case EnvironmentMapFormat.Fisheye
                    assert(~isempty(focal), ...
                        'Fisheye format requires focal length');
                    
                    e.focalLength = focal;
                    
                    % compute the vertical/horizontal fovs
                    e.hFov = (size(e.data,2)/2)/e.focalLength;
                    e.vFov = (size(e.data,1)/2)/e.focalLength;
                    
                    % the mapping function is r = f * theta;
                    e.fisheyeMappingFcn = @(theta) e.focalLength.*theta;
                    % the inverse function is theta = r/f;
                    e.fisheyeInverseFcn = @(r) r./e.focalLength;
                    
                    
                otherwise
                assert(isempty(hfov), ...
                    'Use the "Stereographic" format when specifying a field of view');
            end
            
            % check for input validity
            if e.format == EnvironmentMapFormat.Sphere || ...
                e.format == EnvironmentMapFormat.Angular || ...
                e.format == EnvironmentMapFormat.SkySphere || ...
                e.format == EnvironmentMapFormat.SkyAngular
                assert(size(e.data,1) == size(e.data,2), ...
                    'When the format is sphere/angular, both dimensions must be equal.');
            end
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
            % we'll use the 'hdr' version if available
            if exist('imshowHDR', 'file')
                imshowFcn = @imshowHDR;
            else
                imshowFcn = @imshow;
            end
            imshowFcn(e.data);
        end
        
        function varargout = size(e, varargin)
            varargout{:} = size(e.data, varargin{:}); 
        end
        
        function e = imresize(e, varargin)
            origSize = [e.nrows, e.ncols];
            e.data = imresize(e.data, varargin{:});
            
            if e.format == EnvironmentMapFormat.Stereographic || ...
                    e.format == EnvironmentMapFormat.Fisheye
                
                % must adapt the focal length when resizing!
                ratio = [e.nrows e.ncols]./origSize;
                if ratio(1) ~= ratio(2)
                warning('EnvironmentMap:resize', ...
                    'Aspect ratio changed by %.2f%%', abs(ratio(2)-ratio(1))/ratio(1)*100);
                end
                e.focalLength = e.focalLength * ratio(1);
            end
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
                % rgb2gray does't work for HDR data. 
                e.data = 0.299 * e.data(:,:,1) + ...
                    0.587 * e.data(:,:,2) + ...
                    0.114 * e.data(:,:,3);
            end
        end
        
        function [e, R] = rotate(e, conversionFormat, input, varargin)
            % Rotates the environment map based on various types of input
            % rotation formats
            %
            %   e = e.rotate(conversionFormat, input, <tolerance>, <checks>)
            %
            %   'conversionFormat' can be either:
            %       - 'DCM'  : direct cosine matrix (input = 3x3 matrix)
            %       - 'EA###': euler angles (input = [psi, theta, phi])
            %       - 'EV'   : euler rotation vector & angle (input = [vec angle])
            %       - 'Q'    : quaternion (input = [q1 q2 q3 q4])
            % 
            % See also:
            %   SpinCalc.m
            
            % Get the world coordinates
            [dx, dy, dz, valid] = e.worldCoordinates();
            
            % Get rotation matrix from input
            if ~strcmp(conversionFormat, 'DCM')
                conversion = sprintf('%stoDCM', conversionFormat);
                R = SpinCalc(conversion, input, varargin{:});
            else
                R = input;
            end
            
            % Rotate the data
            ptR = R*[row(dx); row(dy); row(dz)];
            dx = reshape(ptR(1,:), size(dx));
            dy = reshape(ptR(2,:), size(dy));
            dz = reshape(ptR(3,:), size(dz));
            
            % make sure we're still in the [-1,1] interval (might
            % _slightly_ overflow, thus causing problems later on)
            dx = max(min(dx, 1), -1);
            dy = max(min(dy, 1), -1);
            dz = max(min(dz, 1), -1);
            
            % Create new environment map
            [u, v] = e.world2image(dx, dy, dz);
            e = e.interpolate(u, v, valid);
        end
            
        
        function e = convertTo(e, tgtFormat, tgtDim, varargin)
            % Converts the environment map to the specified format
            %
            %   e = convertTo(e, tgtFormat, <tgtDim>, ...)
            %
            %   e = convertTo(e, ..., <par1, val1>, ...)
            %
            % Useful to convert to the Stereographic format, which requires
            % the specification of the field of view (either horizontal,
            % vertical, or diagonal)
            
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
            
            % Create temporary environment map object, and ask it for its
            % image coordinates
            eTmp = EnvironmentMap(tgtDim, tgtFormat, varargin{:});
            [dx, dy, dz, valid] = eTmp.worldCoordinates();
            
            % Put in image coordinates
            [u, v] = e.world2image(dx, dy, dz);

            % Interpolate
            e = e.interpolate(u, v, valid);
            
            % Change format
            e.format = tgtFormat;            
        end
        
        function e = interpolate(e, u, v, valid)
            % Interpolate to get the desired pixel values
            envMap = zeros(size(u, 1), size(u, 2), e.nbands);
            for c=1:size(envMap,3)
                envMap(:, :, c) = reshape(...
                    interp2(linspace(0,1,e.ncols), linspace(0,1,e.nrows), ...
                    e.data(:, :, c), ...
                    u(:), v(:)), size(envMap,1), size(envMap,2));
            end
            valid = valid & any(~isnan(envMap), 3);
            envMap(~valid(:,:,ones(1,e.nbands))) = 0;
            e.data = envMap;
        end
        
        function [x, y, z, valid] = worldCoordinates(e)
            % Returns the [x,y,z] world coordinates for each pixel in the
            % environment map
            [u, v] = meshgrid(linspace(0, 1, e.ncols), linspace(0, 1, e.nrows));
            [x, y, z, valid] = e.image2world(u, v);
        end
        
        function [u, v] = world2image(e, x, y, z)
            % Returns the [u,v] coordinates (in the [0,1] interval)
            %
            %   [u, v] = world2image(e, x, y, z);
            %
            %   [u, v] = world2image(e, pt);
                
            if nargin == 2
                % shortcut for a single point
                assert(length(x) == 3);
                y = x(2); 
                z = x(3); 
                x = x(1);
            end
            
            assert(all(size(x) == size(y)) && all(size(x) == size(z)), ...
                'All of x, y, and z must have the same size');
            
            % Returns the [u, v] image coordinates
            switch (e.format)
                case EnvironmentMapFormat.LatLong
                    [u, v] = e.world2latlong(x, y, z);
                    
                case EnvironmentMapFormat.Angular
                    [u, v] = e.world2angular(x, y, z);
                    
                case EnvironmentMapFormat.Cube
                    [u, v] = e.world2cube(x, y, z);
                    
                case EnvironmentMapFormat.SkyAngular
                    [u, v] = e.world2skyangular(x, y, z);
                    
                case EnvironmentMapFormat.Octahedral
                    [u, v] = e.world2octahedral(x, y, z);
                    
                case EnvironmentMapFormat.Sphere
                    [u, v] = e.world2sphere(x, y, z);
                    
                case EnvironmentMapFormat.SkySphere
                    [u, v] = e.world2skysphere(x, y, z);
                    
                case {EnvironmentMapFormat.Fisheye, ...
                        EnvironmentMapFormat.Stereographic}
                    [u, v] = e.world2fisheye(x, y, z);
                    
                otherwise
                    error('EnvironmentMap:world2image', ...
                        'Unsupported format: %s', e.format.char);
            end
        end
        
        function [x, y, z, valid] = image2world(e, u, v)
            % Returns the [x,y,z] coordinates in the [-1,1] interval
            %
            %   [x, y, z, valid] = world2image(e, u, v);
            %
            
            assert(all(size(u) == size(v)), ...
                'u and v must have the same size');
            
            % Returns the [x, y, z] world coordinates
            switch (e.format)
                case EnvironmentMapFormat.LatLong
                    [x, y, z, valid] = e.latlong2world(u, v);
                    
                case EnvironmentMapFormat.Angular
                    [x, y, z, valid] = e.angular2world(u, v);
                    
                case EnvironmentMapFormat.Cube
                    [x, y, z, valid] = e.cube2world(u, v);
                    
                case EnvironmentMapFormat.SkyAngular
                    [x, y, z, valid] = e.skyangular2world(u, v);
                    
                case EnvironmentMapFormat.Octahedral
                    [x, y, z, valid] = e.octahedral2world(u, v);
                    
                case EnvironmentMapFormat.Sphere
                    [x, y, z, valid] = e.sphere2world(u, v);
                    
                case EnvironmentMapFormat.SkySphere
                    [x, y, z, valid] = e.skysphere2world(u, v);
                    
                case EnvironmentMapFormat.Stereographic
                    [x, y, z, valid] = e.stereographic2world(u, v);
                    
                otherwise
                    error('EnvironmentMap:world2image', ...
                        'Unsupported format: %s', e.format.char);
            end
        end

    end
    
    methods (Access=private)
        % Below are all the conversion functions. 
        % In all of the following: x,y,z,u,v in [0,1]
        
        function [u, v] = world2latlong(~, x, y, z)
            % world -> lat-long
            u = 1 + (1/pi) .* atan2(x, -z);
            v = (1/pi) .* acos(y);
            
            u = u./2; % because we want [0,1] interval
        end
        
        function [u, v] = world2angular(~, x, y, z)
            % world -> angular
            rAngular = acos(-z) ./ (2.*pi.*sqrt(x.^2 + y.^2));
            v = 1/2-rAngular.*y;
            u = 1/2+rAngular.*x;
        end
        
        function [u, v] = world2cube(~, x, y, z)
            % world -> cube
            error('world2cube: fixme');
        end
        
        function [u, v] = world2octahedral(~, x, y, z)
            % world -> octahedral
            error('world2octahedral: fixme');
        end
        
        function [u, v] = world2skyangular(~, x, y, z) 
            % world -> skyangular
            thetaAngular = atan2(x, z); % azimuth
            phiAngular = atan2(sqrt(x.^2+z.^2), y); % zenith
            
            r = phiAngular./(pi/2);
            
            u = r.*sin(thetaAngular)./2+1/2;
            v = 1/2-r.*cos(thetaAngular)./2;
        end
        
        function [u, v] = world2skysphere(~, x, y, z)
            % world -> skysphere
            r = sin(.5.*acos(y)) ./ (sqrt(x.^2+z.^2)) * 2/sqrt(2);
            
            u = .5*r.*x+.5;
            v = 1-(.5*r.*z+.5);
        end
        
        function [u, v] = world2sphere(~, x, y, z)
            % world -> sphere
            r = sin(.5.*acos(-z)) ./ (2.*sqrt(x.^2+y.^2));
            
            u = .5+r.*x;
            v = .5-r.*y;
        end
        
        function [u, v] = world2fisheye(e, x, y, z)
            % world -> stereographic
            assert(~isempty(e.hFov));
            assert(~isempty(e.vFov));
            assert(~isempty(e.focalLength));
            

            theta = pi/2-atan2(-z,sqrt(x.^2 + y.^2));
            phi = atan2(y, x);
            r = e.fisheyeMappingFcn(theta);
%             r = 2.*e.focalLength.*tan(theta/2);
            
            % (u,v) are in pixels
            u = r.*cos(phi);
            v = -r.*sin(phi);
            
            % normalize in the [-1,1] interval
            u = u./e.fisheyeMappingFcn(e.hFov);
            v = v./e.fisheyeMappingFcn(e.vFov);
%             u = u./(2*e.focalLength*tan(e.hFov/2));
%             v = v./(2*e.focalLength*tan(e.vFov/2));
            
            % u and v are now in the [-1,1] interval. 
            u = (u+1)/2;
            v = (v+1)/2;
            
            % camera is pointing in the negative z direction
            u(z > 0) = -1;
            v(z > 0) = -1;
        end
        
        function [x, y, z, valid] = latlong2world(~, u, v)
            u = u*2;
            
            % lat-long -> world
            thetaLatLong = pi.*(u-1);
            phiLatLong = pi.*v;
            
            x = sin(phiLatLong).*sin(thetaLatLong);
            y = cos(phiLatLong);
            z = -sin(phiLatLong).*cos(thetaLatLong);
            
            valid = true(size(x));
        end
        
        function [x, y, z, valid] = angular2world(~, u, v)
            % angular -> world
            thetaAngular = atan2(-2.*v+1, 2.*u-1);
            phiAngular = pi.*sqrt((2.*u-1).^2 + (2.*v-1).^2);
            
            x = sin(phiAngular).*cos(thetaAngular);
            y = sin(phiAngular).*sin(thetaAngular);
            z = -cos(phiAngular);
            
            r = (u-0.5).^2 + (v-0.5).^2;
            valid = r <= .25; % .5^2
        end
        
        function [x, y, z, valid] = cube2world(~, u, v)
            error('cube2world: fixme!');
        end
        
        function [x, y, z, valid] = octahedral2world(~, u, v)
            error('octahedral2world: fixme!');
        end
        
        function [x, y, z, valid] = sphere2world(~, u, v)
            u = u*2-1;
            v = v*2-1;
            
            % sphere -> world
            r = sqrt(u.^2 + v.^2);
            theta = atan2(u, -v);
            
            phi = zeros(size(theta));
            valid = r<=1;
            phi(valid) = 2.*asin(r(valid));
            
            x = sin(phi).*sin(theta);
            y = sin(phi).*cos(theta);
            z = -cos(phi);
        end

        function [x, y, z, valid] = skysphere2world(~, u, v)
            % skysphere -> world
            u = u*2-1;
            v = v*2-1;
            
            valid = u.^2 + v.^2 <= 1;
            
            uM = u*sqrt(2)/2;
            vM = v*sqrt(2)/2;
            
            r = sqrt(uM.^2 + vM.^2);
            
            phi = zeros(size(valid));
            theta = zeros(size(valid));
            
            theta(valid) = atan2(uM(valid), -vM(valid));
            phi(valid) = 2.*asin(r(valid));
            
            x = sin(phi).*sin(theta);
            y = cos(phi);
            z = sin(phi).*cos(theta);
        end
        
        function [x, y, z, valid] = skyangular2world(~, u, v)
            % skyangular -> world
            thetaAngular = atan2(-2.*v+1, 2.*u-1); % azimuth
            phiAngular = pi/2.*sqrt((2.*u-1).^2 + (2.*v-1).^2); % zenith
            
            x = sin(phiAngular).*cos(thetaAngular);
            z = sin(phiAngular).*sin(thetaAngular);
            y = cos(phiAngular);
            
            r = (u-0.5).^2 + (v-0.5).^2;
            valid = r <= .25; % .5^2
        end
        
        function [x, y, z, valid] = stereographic2world(e, u, v)
            % stereographic -> world
            assert(~isempty(e.hFov));
            assert(~isempty(e.vFov));
            assert(~isempty(e.focalLength));
            
            % put back in the [-1,1] interval
            u = u*2-1;
            v = v*2-1;
            
            % scale by image coordinates
            u = u*e.fisheyeMappingFcn(e.hFov);
            v = v*e.fisheyeMappingFcn(e.vFov);
%             u = u*(2*e.focalLength*tan(e.hFov/2));
%             v = v*(2*e.focalLength*tan(e.vFov/2));
            
            r = sqrt(u.^2+v.^2);
            
            % recover spherical angles
            theta = e.fisheyeInverseFcn(r) - pi/2;
%             theta = 2*atan2(r, 2*e.focalLength) - pi/2;
            phi = -atan2(v, u);
            
            % convert back to x, y, z from spherical angles.
            x = cos(theta).*cos(phi);
            y = cos(theta).*sin(phi);
            z = sin(theta);
            
            valid = true(size(u)); % all valid
        end

    end
end

