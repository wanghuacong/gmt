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
[ |-Q|\ *divisor* ]
[ |-S|\ *smoothfactor* ]
[ |-T|\ *foundation* ]
[ |SYN_OPT-V| ]
[ |-Z|\ [**+s**\ *factor*][**+o**\ *shift*] ]
[ |SYN_OPT-f| ]
[ |SYN_OPT--| ]

|No-spaces|

Description
-----------

**grdslice** is used to determine local *peaks* in a 2-D surface.  To do so, it reads a binary 2-D
grid and contour-slices the surface using the selected constant contour interval,
starting from the max grid value and moving down to lower and lower values.  As new closed contours are
first encountered, their centers are considered the location of a new peak, and we trace the hierarchy of how
each new slice belongs to groups of previously detected peaks by virtue of enclosing them.  Slice areal shapes
are equated to ellipses and their best-fitting major and major axes and orientations are computed.
A rich set of information is written to separate data files, as requested via options.  The primary
output (written to *stdout*) is a table with the vertical traces of peak centers per contour. Options
are available to set limits on how small contour areas to consider and if we should eliminate closely-spaced
peaks, for instance.

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
    Skip contours whose areas are less than *area* in km^2 (or squared user units
    for Cartesian grids) [retain all closed contours].

.. _-D:

**-D**\ *distance*
    Examine the peak locations found and eliminate those that fall within *distance* km (or Cartesian units)
    of another peak with a larger amplitude [retain all peaks].

.. _-E:

**-E**\ *slicefile*
    Set the filename for the output of all slice polygons found [no slice output].  Each segment header
    will list the contour area (via **-Z**), the contour level (via **-L**), the number of identified peaks
    inside this contour (via **-N**), the parameters of the best-fitting ellipsoid (via **-S**\ *x/y/azim/major/minor/fit*),
    and the *ID* of the foundation contour it belongs to (via **-I**). The *fit* is a measure between 0-100 of now
    close the contour is to the best-fit ellipsoidal shale of the same area. Each segment header is followed by the
    coordinates of the contour via *x y z* records, with *z* being the contour value.

.. _-F:

**-F**\ *foundationfile*
    Set the filename for the output of all foundation polygons [no foundation output]. Each segment header
    will list the contour area (via **-Z**), the contour level (via **-L**), the number of identified peaks
    inside this contour (via **-N**), the parameters of the best-fitting ellipsoid (via **-S**\ *x/y/azim/major/minor/fit*),
    and the *ID* of the foundation contour it belongs to (via **-I**). The *fit* is a measure between 0-100 of now
    close the foundation contour is to the best-fit ellipsoidal shale of the same area. Each segment header is followed by the
    coordinates of the contour via *x y z* records, with *z* being the contour value.

.. _-I:

**-I**\ *indexfile*
    Set the filename for the output of all peak indices [no index output].  The format of this file is a set of
    *x y *index* records for the peaks.

.. _-L:

**-L**\ *low/high*
    Limit contour slicing range: Only do the contour slicing within the specified range of contours [use the entire range].

.. _-Q:

**-Q**\ *divisor*
    The location and amplitude of each peak is refined by a grid search on a sub-pixel lattice inside the closed
    uppermost contour.  The *divisor* is used to subdivide the grid intervals by that factor.  A higher value yields
    a more precise result.  **Note**: The interpolated grid value at the sub-pixel locations depends on the interpolation
    selected via **-n** [8].

.. _-S:

**-S**\ *smoothfactor*
    Resample all contour lines at roughly every (*gridbox_size*/*smoothfactor*) interval [no along-contour resampling].

.. _-T:

**-T**\ *foundation*
    Select the *foundation* contour level for detecting a peak [Default uses all, but see **-L**].
    **Note**: If **-L** is used, the *foundation* must be equal or larger than the *low* value.

.. _-Z:

**-Z**\ [**+s**\ *factor*][**+o**\ *shift*]
    Subtract *shift* from the data grid and multiply the results by
    *factor* before contouring starts [1/0]. (Values in **-C**, **-L**
    and **-T** refer to values *after* any scaling has occurred.)

.. |Add_-f| unicode:: 0x20 .. just an invisible code
.. include:: ../../explain_-f.rst_

.. |Add_-V| unicode:: 0x20 .. just an invisible code
..  include:: ../../explain_-V.rst_

.. include:: ../../explain_colon.rst_

.. include:: ../../explain_help.rst_

Examples
--------

We illustrate a typical workflow by using a synthetic data set with five peaks reflecting
elliptical Gaussian shapes created by::

    gmt grdseamount -R-0:30/0:30/-0:30/0:32 -I1m -Gsynth.grd -Cg -E -fg <- EOF
    0       0       30  80  40  131
    -0.3    -0.3    70  25  16  60
    0.3     0.3     50  25  25  60
    -0.35   0.35    60  25  20  45
    0.35    -0.35   -30 25  15  47
    EOF

The synthetic grid is visualized via :doc:`grdview </grdview>` below:

.. figure:: /_images/GMT_grdslice_view.*
   :width: 500 px
   :align: center

   Three-dimensional view of a cluster of five peaks, two which are isolated and three which
   are part of an overlapping cluster.

We can run :doc:`grdcontour </grdcontour>` to visualize what the data looks like in map-view:

.. figure:: /_images/GMT_grdslice_contours.*
   :width: 400 px
   :align: center

   The contour map reveals in 10-unit increments what the shapes look like.  For this quick
   analysis we decide to use a 10-unit interval and use *z = 10* as the foundation contour.

To determine the centers, slices, foundations and index files for this data we try::

    gmt grdslice synth.grd -C10 -Isynth_index.txt -Esynth_slice.txt -Fsynth_foundation.txt -T10 -L10/130 -A5 -nl -Q8 -S8 > synth_center.txt

The foundations, indices, slices and best-fitting ellipses are then presented graphically below:

.. figure:: /_images/GMT_grdslice_products.*
   :width: 750 px
   :align: center

   (left) We paint the three foundations and label them with the index numbers. (middle)
   All the contour slices are reproduced and we plot the five peaks that were identified.
   (right) Each contour is approximated by an ellipse whose parameters are used to make
   this plot.

Finally, one can use the contour centers to visualize how these peaks might be related
hierarchically:

.. figure:: /_images/GMT_grdslice_connect.*
   :width: 500 px
   :align: center

   We attempt to show the 3-D lines connecting the slice centers as the start at each
   peak but then merge as peaks start to share common contours, all the way down to
   the foundation contour.

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
