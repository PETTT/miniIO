# miniIO

miniIO is a series of mini-apps that provide the data structures necessary 
to reproduce the I/O behavior of representative types of HPC simulation 
applications.  The computational code for each mini-app is only enough to fill
the data structures with reasonable values and scales trivially in MPI. This
simulates the expected growth in scaling of codes without having to wait for
such codes to be optimized or refactored, allowing one to perform I/O
scaling runs to maximum system or queue sizes.

The mini-apps include stubbed interface methods and make files to add I/O
libraries or routines in modular fashion. Several I/O libraries/methods are 
already supported, ADIOS, HDF5, and MPI-IO (not in all mini-apps).  This allows
users of the mini-apps to:

- add new I/O types to the mini-apps easily, 
- switch between I/O types at will for benchmarks, 
- disable each I/O library component if it is not available for a given system. 

The mini-apps also share a simplectic noise method to generate 
pseudo-realistic synthetic data to fill data structure arrays so that one can 
produce reliable performance results with I/O compression algorithms.

There are currently 4 mini-apps:

### 1. struct

The struct mini-app produces Cartesian, rectilinear, and curvilinear grids.
Structured grids have regular, implicit connectivity, resulting in a constant
number of points along each dimension with two neighboring points on each
dimension (except edges).  However, these structured grids can include a
_missing_ data mask that marks grid points as somehow unavailable for
computation.  A real-world example would be land points in an ocean model. Masked points are important because they create load imbalance in a simulation. Simulation codes often load balance their domain decomposition for better scalability. Unfortunately this unbalances the I/O load since masked data points still get written out. This mini-app helps characterize the effects of this I/O load imbalance on I/O performance.

### 2. unstruct

The unstruct mini-app produces unstructured grids.	Unstructured grids do not
have implicit connectivity.  Any given data point could be connected to any
nearby point and a special data structure of connections must be provided (or
inferred from base geometry in the case of strand grids).  The grids themselves
can be steady or time-variant.  The grid generation is implemented using a
superquadric equation whose _roundness_ parameters can be varied over time to
emulate a time-variant grid.  The grid is also constructed as a series of
volume-inclusive shells so that both traditional unstructured grids and newer
strand grids can be represented.

### 3. amr

Adaptive Mesh Refinement (AMR) grids adapt their resolution to fit a solution
space as appropriate.  AMR often starts with a coarse Cartesian base grid. Then,
when some parameter of the model at a given location indicates that higher
resolution is required, that region or single cell of the grid is subdivided
into smaller Cartesian sub-grids. 

### 4. cartiso

Derived output occurs when analysts that use a simulation code wish to pull out
subsets of the simulation data for post-hoc analysis.  Since it is difficult to
predict where these derivations lie within the model data space, often the
derivations are very load imbalanced.  For highly scalable models, the load
imbalance across tens of thousands of compute cores reduces performance
significantly.  The cartiso mini-app is named as such because it generates a Cartesian volumetric grid as the basis for computing an isosurface. As such, it is two mini-apps in one, as one can write out either the Cartesian arrays or the isosurface geometry or both.

## Contents Overview

- LICENSE - BSD 3-clause license text
- Makefile.inc - common included makefile
- pdirs.h, timer.h - common definitions, functions, etc.
- osn - OpenSimplexNoise library; used by all mini-apps
- struct, unstruct, amr, cartiso - the mini-apps!
    - README - all mini-apps provide their own README
    - Makefile - all mini-apps use a similar Make approach
    - {miniapp}.c - mini-app main is in .c file of its own name
    - adios*.c|.h - ADIOS I/O methods, as supported
    - hdf5*.c|.h - HDF5 I/O methods, as supported
    - Other I/O methods may or may not be the same across mini-apps.

## Building

For now, there is no single Makefile for all mini-apps.  One enters the
directory of each desired mini-app to build.  In general, just follow the README
in each mini-app directory.  The general approach is:

1. Edit Makefile.inc and set the correct compilers and flags for the system.
2. Enter desired mini-app directory.
3. Edit Makefile to enable/disable I/O modules
4. ``make``

Repeat as desired for each mini-app.  The mini-apps are aware of their shared
dependencies (e.g., osn), so they will build it if necessary.

### Adding HDF5 support

In order to build the mini-apps with HDF5 support, simply set the environment variable
*HDF5_DIR* to the directory where HDF5 is installed. Additionally, the HDF5 output
can be generated by setting the command line option --hdf5 for most of the
mini-apps (the exception being the _cartiso_ app).


## Running

The mini-apps use a common command-line argument format where possible.  In
general, a run will look something like this:
```sh
mpiexec -np N ./miniapp --required_arg option [--optional_arg option]
```
where:
- mpiexec = actual executable depends on your MPI implementation
- N = number of MPI tasks
- miniapp = name of actual mini-app
- required_arg = some required command-line argument; varies per mini-app
    - option = some arguments require options after them
- optional_arg = optional command-line argument

See the README and USAGE statement for each mini-app for more details.

## Adding an I/O Module to a Mini-app

Follow the README for each mini-app for specifics.  Each mini-app carefully and
clearly comments where new I/O methods get added.  The general overview of
adding an I/O module for each mini-app is as follows:

1. Add to the I/O Modules sections in the Makefile
2. Add includes to main .c file
3. Recommended: add usage string to main .c file
4. Recommended: Add state variables; you'll need some for command line options
5. Recommended: add your command line options
6. Optional: Add initialization, if your code will need it
7. Add output routine or function call in marked section of time step loop
8. Optional: Add cleanup/finalization code
9. Update Make dependencies (``make depend``) if new files were added

If any of the above steps, especially initialization (#6), time step output
(#7), or cleanup (#8), are longer than 2-3 lines, please put them in a function
in a separate .h and .c file (and update dependencies).  This keeps the main
mini-app code much cleaner.  We won't accept pull requests or submitted code
that doesn't follow this standard.

The "recommended" steps aren't necessary to get working code.  They are
recommended to provide proper options and documentation to users of your code.
We won't accept pull requests or submitted code that don't do these.

