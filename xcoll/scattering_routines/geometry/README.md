# Developer Notes on Xcoll Geometry

The geometry code is meant to be as independent as possible, as to allow a flexible way to define new geometries.
It provides all the necessary toolery to do a 2.5D geometry: a 2D geometry in the $(s, x')$ plane (where $x'$ is already rotated over the collimator angle) with optional vertical bounds.
The main purpose of the code is to find the crossings of a given particle `trajectory` with a given `Shape2D`: a very flexible way to construct any shape from a set of `Segment`s.


## Trajectories
A trajectory is a particle moving along the $s$ direction. As such it should be a well-defined function of $s$ which returns $x$ (perpendicular trajectories, that would move along $x$ at constant $s$, are hence not allowed).

The textbook (and currently only implemented) trajectory is the `drift`, defined as a line with slope $m_x$ through a point $(s_0, x_0)$

## Segments
are defined parametrically (e.g. a perpendicular line segment can exist), e

## Shapes

Pretty straightforward to define a new shape
To define a new type of object, one has to:
1. Define the `create_` and `destroy_` functions in objects.h
2. It might be useful to define an API to link it to Xsuite as done in collimator_geometry.h
3. Add tests for this new object (and potentially for its API)

## Crossing Functions in C

## Defining New Segment Types
To define a new type of segment, one has to:
1. Create the struct definition and ID macro - make sure the ID is unique (in segments.h)
2. Define a constructor function (in segments.h)
3. For each trajectory type, define a crossing function for that segment and update the general function `crossing_...` (in crossing_....h). Roots with multiplicity should be added a number of times equal to the multiplicity (this is because the algorithm will count IN/OUT trajectories based on the ordered roots)
4. For each trajectory type, adapt `max_array_size_..` to give the maximum number of crossings that can occur (in crossing_....h)
5. Add tests for this new segment

add to xcoll init

code auto-generation in Shape2D

## Defining New Trajectories
To define a new type of trajectory (e.g. multiple coulomb scattering), one has to:
1. Create a new file crossing_trajectoryname.h and add it to geometry.h
2. For all segments, define the crossing functions and the general `crossing_trajectoryname`,`crossing_trajectoryname_vlimit`, and `max_array_trajectoryname_drift` (as done in crossing_drift.h)
3. Add the four `s` functions (`crossing_trajectoryname_first`, `crossing_trajectoryname_vlimit_first`, `crossing_trajectoryname_after_s`, `crossing_trajectoryname_vlimit_after_s`) in get_s.h
4. Add tests for this new trajectory


## Outlook

One could think of using the `Segment` definitions to define a `FrontShape` class to define collimators with a particular transversal profile, but this is currently not yet implemented.



## Files

**sort.h:**
A basic low-level fast implementation to sort arrays of ints or floats.

**methods.h:**
A collection of useful functions. Currently there is only one: to calculate the overlap between two intervals.

**find_root.h:**
First steps towards a smart implementation of a root finder.
Not yet used.

**segments.h:**
The segment definitions, as polymorphous variants of the `Segment` parent struct.
Currently there are a `LineSegment`, a `HalfOpenLineSegment`, and a `CircularArcSegment`.

**objects.h:**
Objects are defined as arrays of `Segments`.
For each object there are a constructor function that allocates the memory and returns the array (`create_objectname`) and a destructor that frees the memory (`destroy_objectname`).
Currently there are a `jaw`, a `polygon`, an `open_polygon`, and a `crystal`.

**crossing_drift.h:**
This file contains, for each segment, the function that defines how to find the crossing of a drift trajectory with that segment.
These need a very specific syntax:
`crossing_drift_segmentname(void* segment, int8_t* n_hit, double* s, _driftargs_)`
The void pointer points to the segment, `n_hit` and `s` are accumulation variables that collect the crossings, and `_driftargs_` are the arguments that define the trajectory (`double part_s, double part_x, double part_tan` for a drift trajectory).
Finally, two general crossing functions are provided:
- `crossing_drift(segments, n_segments, n_hit, s, _driftargs_)` which gives all crossings of a drift trajectory with an object defined by an array of `Segments`
- `crossing_drift_vlimit(segments, n_segments, n_hit, s, _driftargs_, _driftargs_V_, y_min, y_max)` which gives all crossings of a drift trajectory with an object, with an additional vertical constraint given by `y_min` and `y_max` (which is why we need `__driftargs_V__`)

**get_s.h:**
For each trajectory type (currently only drift), there are four functions:
- `crossing_drift_first(segments, n_segments, _driftargs_)` which gives the s-value of the first crossing
- `crossing_drift_vlimit_first(segments, n_segments, _driftargs_, _driftargs_V_, y_min, y_max)` as above but with a vertical constraint
- `crossing_drift_after_s(segments, n_segments, _driftargs_, current_s)` which gives the s-value of the first crossing after a given `current_s`
- `crossing_drift_vlimit_after_s(segments, n_segments, _driftargs_, _driftargs_V_, y_min, y_max, current_s)` as above but with a vertical constraint


### The following files are dependent on the Xsuite API:

**rotation.h**:
A `YRotation` without the drift.

**collimator_geometry.h**:
API to simplify the integration with collimators in Xsuite: a struct to pass the geometry information (as we cannot call the Xsuite functions at this level), a function to check the impact locations and move to the collimator jaw frame, and a function to return from the jaw frame.

**crystal_geometry.h**:
The same but for crystals.

**geometry.py**:
An Xobject that contains all C dependencies, such that `BeamElements` can be dependent on them and automatically get the C code.


## Defining new segments



## Defining new objects
