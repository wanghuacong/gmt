.. index:: ! grdslice
.. include:: ../module_supplements_purpose.rst_

********
grdslice
********

|grdslice_purpose|

Synopsis
--------

.. include:: ../../common_SYN_OPTs.rst_

**gmt grdslice**
*grid*
|-C|\ *interval*
[ |-A|\ [*area*] ]
[ |-D|\ [*distance*] ]
[ |-E|\ [*slicefile*] ]
[ |-F|\ [*foundationfile*] ]
[ |-I|\ [*indexfile*] ]
[ |-L|\ *low/high* ]
[ |-S|\ *smoothfactor* ]
[ |-T|\ *foundation* ]
[ |SYN_OPT-V| ]
[ |-Z|\ [**+s**\ *factor*][**+o**\ *shift*] ]
[ |SYN_OPT-f| ]
[ |SYN_OPT--| ]

|No-spaces|

Description
-----------

**grdslice** is used to determine local *peaks* in a 2-D grid.  To do so, it reads a binary 2-D
grid and contour-slices the surface using the selected constant contour interval,
starting from the max grid value and moving down to lower and lower values.  As new contours are
first encountered, their centers are considered the location of a new peak, and we trace the hierarchy of how
each new slice belongs to groups of previously detected peaks by virtue of enclosing them.  Slice shapes
are equated to ellipses and their best-fitting major and major axes and orientations are computed.
Information are written to separate data files, as requested.  The primary output (to stdout) is a table
with the vertical traces of peak centers.

Required Arguments
------------------

*grid*
    2-D gridded data set to be sliced. (See GRID FILE FORMATS below). We automatically detect if this
    grid is geographic or in spherical Mercator units (such as created by :doc:`img2grd </supplements/img/img2grd>`.
    If plain Cartesian you may need to use modifiers to convert horizontal distance units to km.  E.g.,
    if your Cartesian grid has *x* and *y* in meters then append **+Uk** to auto-convert to km.

.. _-C:

**-C**\ *interval*
    Constant contour *interval* for slicing.  The smaller the interval, the more detailed the results
    and the longer the execution time.

Optional Arguments
------------------

.. _-A:

**-A**\ *area*
    Ignore contours whose areas are less than *area* in km^2 (or user units
    for Cartesian grids) [retain all contours].

.. _-D:

**-D**\ *distance*
    Examine the peak locations found and eliminate those that fall within *distance* km of another
    peak with a larger amplitude [retain all peaks].

.. _-E:

**-E**\ *slicefile*
    Set the filename for the output of all slice polygons [no slice output].

.. _-F:

**-F**\ *foundationfile*
    Set the filename for the output of all peak foundation polygons [no foundation output].

.. _-I:

**-I**\ *indexfile*
    Set the filename for the output of all peak indices [no index output].

.. _-L:

**-L**\ *low/high*
    Limit contour range: Only do the contour slicing within the specified range of contours [use he entire range].

.. _-S:

**-S**\ *smoothfactor*
    Resample the contour lines at roughly every (gridbox_size/*smoothfactor*) interval [no resampling].

.. _-T:

**-T**\ *foundation*
    Select the *foundation* contour level for detecting a peak [Default uses all, but see **-L**].
    **Note**: If **-L** is used, the *foundation* must be equal or larger than the *low* value.

.. _-Z:

**-Z**\ [**+s**\ *factor*][**+o**\ *shift*]
    Use to subtract *shift* from the data and multiply the results by
    *factor* before contouring starts [1/0]. (Numbers in **-C**, **-L**
    and **-T** refer to values *after* this scaling has occurred.)

.. |Add_-f| unicode:: 0x20 .. just an invisible code
.. include:: ../../explain_-f.rst_

.. |Add_-V| unicode:: 0x20 .. just an invisible code
..  include:: ../../explain_-V.rst_

.. include:: ../../explain_colon.rst_

.. include:: ../../explain_help.rst_

Examples
--------

To slice the grid gravity_smts.grd every 0.1 mGal and write all results, try::

    grdslice gravity_smts.grd -C0.1 -Dindeces.txt -Eslices.txt -Ffoundations.txt -V > centers.txt

References
----------

Kim, S.-S., and P. Wessel (2011), New global seamount census from altimetry-derived gravity data,
*Geophys. J. Int., 186*, 615-631,
`http://dx.doi.org/10.1111/j.1365-246X.2011.05076.x <http://dx.doi.org/10.1111/j.1365-246X.2011.05076.x>`_.

Wessel, P., D. T. Sandwell, and S.-S. Kim (2010), The global seamount census, *Oceanography, 23(1)*, 24-33,
`http://dx.doi.org/10.5670/oceanog.2010.60 <http://dx.doi.org/10.5670/oceanog.2010.60>`_.

See Also
--------

:doc:`gmt.conf </gmt.conf>`, :doc:`gmt </gmt>`,
:doc:`grdcontour </grdcontour>`,
:doc:`grdmath </grdmath>`,
:doc:`grdseamount </supplements/potential/grdseamount>`,
:doc:`img2grd </supplements/img/img2grd>`
