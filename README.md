This package contains MATLAB functions that convert to and from different types 
of environment map representations. Most of the code was implemented 
following the book:

> E. Reinhard, G. Ward, S. Pattanaik, and P. Debevec. High dynamic range 
> imaging. Morgan Kaufman, 2005.

Please see the book for more details. 

Features
========

The main features of this package are:
1. Support for many different environment map formats
2. Easy integration of additional formats
3. Automatic computation of solid angles

Formats
=======

The supported formats are:

- *Sphere*: spherical format, useful for mapping from a photo of a metallic
sphere. Rays are very distorted close to the edge of the circle. 

- *Angular*: typical format for storing environment maps, because it minimizes
distortions and is relatively compact;

- *LatLong*: latitude-longitude (equirectangular).

- *Cube*: 6-faced representation. Commonly used. 

- *SkyAngular*: convenient way to visualize the top hemisphere of the 
environment map only, with minimal distorsions. _Warning_: converting to this 
format will entirely drop the bottom half of the environment map!

- *Octahedral*: projection of a sphere onto an octahedron. 

Coordinates
===========

The world coordinates `(dx,dy,dz)` have the following reference frame:

<pre><code>
     ^ y
     |
     |
     .----> x
    /
z  v
</code></pre>
  
The camera would be looking in the negative z direction.

The EnvironmentMap class
========================

The `EnvironmentMap` class provides unified functionality and renders the 
use of the individual functions obsolete. It makes the creation and manipulation
of environment maps very easy! Here's an example:

1. Load an environment map from file (this assumes the file stores the 
environment map in the latitude-longitude format):

        $ envmap = EnvironmentMap('path/to/file', 'latlong');

2. Resize the environment map:

        $ envmap = imresize(envmap, [500 NaN]);

3. Display the environment map:

        $ imshow(envmap);

4. Convert to the 'angular' representation:

        $ envmap = envmap.convertTo('sphere');
        $ imshow(envmap);

5. Rotate the environment map (now in angular representation) by 90 degrees 
around the x-axis (see `SpinCalc.m` in the `3rd_party` directory for more on
the nomenclature for representing rotations):

        $ envmap = envmap.rotate('EA213', [90 0 0]);
        $ imshow(envmap);

The EnvironmentMapFormat class
------------------------------

The `EnvironmentMapFormat` class is used as an `enum` type to store all supported
formats (see above). For convenience, it contains a `format()` function which
allows conversion from strings. For example:

    $ fmt = EnvironmentMapFormat.format('sphere');

Strings don't have to match the format perfectly, the closest string (according
to the Levenshtein distance) is used. 

Dependencies
============

This software package depends on the following (included) 3rd-party libraries
(see the `3rd_party` directory):

- [hdrutils](https://github.com/lvsn/hdrutils)

- `SpinCalc` for rotation representations, by John Fuller 
([http://www.mathworks.com/matlabcentral/fileexchange/20696-function-to-convert-between-dcm-euler-angles-quaternions-and-euler-vectors](http://www.mathworks.com/matlabcentral/fileexchange/20696-function-to-convert-between-dcm-euler-angles-quaternions-and-euler-vectors))

Notes
-----

All functions assume that their inputs is in double format in the [0, 1] 
range. `im2double` should do the trick.


History
=======

- 10/30/12: Moved to github. Refer to commit messages for updates.

- 04/19/11: Fixed small but in envmapWorld2SkyAngular.m which prevented 
  correct conversion from skyAngular to latlong format.
