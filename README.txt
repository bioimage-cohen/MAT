The multitemporal association tracker is a multitarget tracking system, this
version has a cost function which was developed for tracking cells in 2D,
though it would not be difficult to modify the cost function for 3D.

This directory includes Visual Studio 2008 solution files, but the underlying
tracker should build on other platforms if necessary.


Directory Structure
-------------------
- root: README and project files
|
|-- src: All source files for the MAT cell tracker
|
|-- samples: Some sample MATLAB read/write code and sample data


Input File Format
-----------------
The tracker reads detection (segmentation) information using a simple text
format. If a more complex or faster format is required, or new detection info
must be integrated into the tracker then all the related codecan be found in
detection.cpp. A sample MATLAB writer for the current format, along with some
sample data is included in the samples subdirectory.

Here is an abstract example of the file format:
Note: <> indicate the description of a single text number in the file format

<imW> - The width of the source image or movie (in pixels)
<imH> - The height of the souce image or movie (in pixels)
<numT> - The total number of frames in the source movie
<numD> - The total number of detections found over all frames
<numD_#> - The total number of detections found in frame #
<comX_$> - The centroid X coordinate of detection $ in frame #
<comY_$> - The centroid X coordinate of detection $ in frame #
<numP_$> - The number of pixels in connected-component of detection $ in frame #
<pixX_@> - The X coordinate of pixel @ in detection $ in frame #
<pixY_@> - The X coordinate of pixel @ in detection $ in frame #

---------------------
<imW> <imH>
<numT> <numD>

<numD_1>
<comX_1> <comY_1> <numP_1>: (<pixX_1>,<pixY_1>) (<pixX_2>,<pixY_2>) (<pixX_2>,<pixY_2>) ...
<comX_2> <comY_2> <numP_2>: (<pixX_1>,<pixY_1>) (<pixX_2>,<pixY_2>) (<pixX_2>,<pixY_2>) ...
.
.
.
<numD_N>
<comX_1> <comY_1> <numP_1>: (<pixX_1>,<pixY_1>) (<pixX_2>,<pixY_2>) (<pixX_2>,<pixY_2>) ...
.
.
.


Output File Format
------------------
The MAT tracker outputs two important structures to a text file, first the
edges which were assigned during tracking and a global ID which indicates to
which "track" they were assigned.

We also output the tracking graph which contains the costs the multitemporal
tracking assigned to each feasible edge (connection) between detections in 
nearby frames.

Here is an abstract example of the file format:
Note 1: <> indicate the description of a single text number in the file format
Note 2: ID fields are global detection IDs, these begin at 0 and are assigned
        to the detections in the order they are read from the input file.

<inID> - Incoming vertex ID of a directed edge (ie. destination vertex)

<outID> - Outgoing vertex ID of a directed edge (ie. source vertex)

<labelT> - Global track label, numbers are arbitrary but all edges with the same
           label are on the same track

<cost> - Edge connection cost in the tracking graph

------------
<labelT>,<outID>,<inID>
<labelT>,<outID>,<inID>
<labelT>,<outID>,<inID>
.
.
.
-1,-1,-1
<outID>,<inID>,<cost>
<outID>,<inID>,<cost>
<outID>,<inID>,<cost>
.
.
.

Thus, reading until the -1,-1,-1 line will give all the necessary tracking info
but reading further gives costs and information for other edges which were not
assigned to tracks.