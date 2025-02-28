----------------------------------------------------------------
Package name: LpCVT
----------------------------------------------------------------
The essential code of our LpCVT is provided in this package. It
contains Restricted Voronoi diagram computation in a bounded
domain, L_p Lloyd energy computation on surface/volume.
----------------------------------------------------------------

Usage: LpCVT data/three_holes.obj data/three_holes.pts

>Combinatorial tests:
  The program outputs the restricted Voronoi diagram (rvd.obj), the restricted Delaunay triangulation (rdt.obj) and 
the clipped Voronoi diagram (cvd.obj).
>Algebraic tests:
  The program outputs the value and gradient of F_Lp in surface and volume mode.

----------------------------------------------------------------
Installation:

Prerequisites:
    CGAL library (http://www.cgal.org);
Build Tool:
    CMake (http://www.cmake.org)
Tested Compilers:
    Microsoft Visual C++ 9.0
    GNU G++

Generating Makefile:

For both Windows and Linux users:
    you can generate Visual C++ project file or Makefile by running the following
command on LpCVT directory
    cmake-gui .

For Linux user:
    A stand-alone Makefile is provided also.

For Windows user:
    If you have not any experience on CGAL installation, please
    follow the instruction in INSTALL_WIN.txt.

----------------------------------------------------------------
