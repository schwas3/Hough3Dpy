This is a python implementation of the source code associated with (https://doi.org/10.5201/ipol.2017.208):

Christoph Dalitz, Tilman Schramke, and Manuel Jeltsch, Iterative Hough Transform for Line Detection in 3D Point Clouds, Image Processing On Line, 7 (2017), pp. 184–196. https://doi.org/10.5201/ipol.2017.208


The exact performance of the code will eventually be modifed to be application specific and NOT exactly replicate the behavior of the source code.

-----

Code is currently written to require manual changing of input parameters in `hough3dlines.py` (and in other files??). The output is textbased only, `python hough3dlines.py > output_file.out` will produce output to a file.

To-do:
- Argparse
- Visualization (+boolean arg, output file arg, more parameter args)
- Output points in lines (either ID's or lists of points or ? maybe some intergration with flow is possible?)
- Output confidence values? i.e. eigenvalues or LSQ-fit value to identify line strengths? (Note: this is obtainable with a description of a line and list of points, resource scope exceeded?)
- Optimization of sphere.py - don't allow 3 forbidden vertices to produce a vertex (this is an important nuance, forbidden vertices and allowable vertices may produce allowable vertices)
  - Not sure if this is actually possible... Triangles can be forbidden but can vertices? (EDIT: Vertex trees extend outwards, meaning this optimization is forbidden without substantial consideration - i.e. not worth doing considering limited benefit)
  - sphere.py is likely not as impactful on performance as the granularity itself (i.e. making the list is cheap, applying it is expensive and not impacted by sphere.py)
- Add energy dependent functionalities? Not likely, but could have minEnergy? Again, this likely exceeds the scope of this resource.
- Optimize code by using list indices instead of lists.
  - Rather than initialize `X` and remove candidate points in `Y` from `X`, include a boolean flag in `X` and `Y` (`X` init to all `1` and `Y` to all `0`) and perhaps make the flag a line index in `Y` (would instead init `Y` to all `-1` and assign line index as they are found)?
  - ?
- **Remove sphere.py in favor of sphere.dat - a set of pregenerated vertices, will be useful when Hough3D-py is implemented as a step in production**
  - `sphere.py` benchmark testing
  - `sphere.dat` file?
  - EDIT: should be able to remove useless triangles

----
Overview

`sphere.py`
Generates vectors in a manner to maximize the their uniform distribution about the surface of the unit sphere.
- Starts with 12 vertices and 20 triangles of a unit icosahedron
- All triangles are tesselated (1 triangle => 4 smaller triangles, defined by reference to 3 vertex indices), then vertices are normalized to put them on the unit sphere, repeat `SphereGranularity` times
- Degenerate vectors (`v=-v`) are removed

`pointcloud.py`
Defines `PointCloud` class:
- `PointCloud.shift` stores total displacement of origin from true `<0,0,0>`
- `PointCloud.points` a `Nx3` numpy array of points
- `getMinMax3D` finds extreme `x`, `y`, `z` coordinates
- `shiftToOrigin` shifts `PointCloud.points` s.t. the Min/Max points found by `getMinMax3D` are equidistant from the origin. Updates `PointCloud.shift` accordingly.
- `meanValue` returns CoM (or `<0,0,0>` if empty) of all unweighted points
- `readFromFile` sets `PointCloud.points` from a file
- `pointsCloseToLine(a,b,opt_dx)` gets `PointCloud.points` < `opt_dx` away from line `ab`
- `removePoints(Y)` removes `Y.points` from `PointCloud.points`

`hough.py`
Stores accumulator array (`VotingSpace`) and performs geometric algorithms described in reference paper to extract most favored line (or end code execution)

`hough3Dlines.py`
Main run code. Also includes orthogonal least-square fit procedure definition (`orthogonal_LSQ`). Details can be found in code/IPOL-article. (This section to be updated?)
