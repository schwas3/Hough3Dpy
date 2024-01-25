This is a python implementation of the source code associated with (https://doi.org/10.5201/ipol.2017.208):

Christoph Dalitz, Tilman Schramke, and Manuel Jeltsch, Iterative Hough Transform for Line Detection in 3D Point Clouds, Image Processing On Line, 7 (2017), pp. 184–196. https://doi.org/10.5201/ipol.2017.208


The exact performance of the code has been modifed to be application specific and NOT exactly replicate the behavior of the source code.

-----
To-do?:
- Output points in lines (either ID's or lists of points or ? maybe some intergration with flow is possible?)
- Output confidence values? i.e. eigenvalues or LSQ-fit value to identify line strengths? (Note: this is obtainable with a description of a line and list of points, resource scope exceeded?)
- Add energy dependent functionalities? Not likely, but could have minEnergy? Again, this likely exceeds the scope of this resource.
- Optimize code by using list indices instead of lists.
  - Rather than initialize `X` and remove candidate points in `Y` from `X`, include a boolean flag in `X` and `Y` (`X` init to all `1` and `Y` to all `0`) and perhaps make the flag a line index in `Y` (would instead init `Y` to all `-1` and assign line index as they are found)?
  - ?

----
Overview

`sphere.py`
Generates vectors in a manner to maximize the their uniform distribution about the surface of the unit sphere.
- Starts with 12 vertices and 20 triangles of a unit icosahedron
- All triangles are tesselated (1 triangle => 4 smaller triangles, defined by reference to 3 vertex indices), then vertices are normalized to put them on the unit sphere, repeat `SphereGranularity` times
- Degenerate vectors (`v=-v`) are removed
- A `sphere.pkl` file is used to store `uniqueVectors`, `vectors`, and `triangles` in a manner that automatically pulls from the best source. `sphereVertices(...)` can also take the argument `saveToFile=False` and will recompute sphere direction vectors procedurally instead. If encountering any issues, consider deleting `sphere.pkl` and rerunning to refactor this cache file.

`hough3Dlines.py`
Main run code. Also includes orthogonal least-square fit procedure definition (`orthogonal_LSQ`). Details can be found in code/IPOL-article. (This section to be updated?)
