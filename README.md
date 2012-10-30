Readme
======

This package contains functions that convert to and from different types 
of environment map representations. Most of the code was implemented 
following the book:

> E. Reinhard, G. Ward, S. Pattanaik, and P. Debevec. High dynamic range 
> imaging. Morgan Kaufman, 2005.

Please see the book for more details. 

The "sky angular" format is a convenient way to visualize the top 
hemisphere of the environment map only, with minimal distorsions. 

Coordinates
-----------

The world coordinates `(dx,dy,dz)` have the following reference frame:

       ^ y
       |
       |
       .----> x
      /
  z  v
  
The camera would be looking in the negative z direction.

Note
----

All functions assume that their inputs is in double format in the [0, 1] 
range. `im2double` should do the trick.

History
-------

10/30/12: Moved to github. Refer to commit messages for updates.

04/19/11: Fixed small but in envmapWorld2SkyAngular.m which prevented 
correct conversion from skyAngular to latlong format.