# LpCVT
![MacOS](https://github.com/Bigger-and-Stronger/LpCVT/actions/workflows/macos-build.yml/badge.svg)

This repo is a reorganization of the source code of the paper "**Lp Centroidal Voronoi Tessellation and its applications**". [[Paper]](https://dl.acm.org/doi/abs/10.1145/1778765.1778856)

The original source code is from https://xueyuhanlang.github.io, and the code has been rearranged by [secantalpha](https://github.com/secantalpha), with some modifications made to ensure compatibility with CGAL 5.0 or a later version. The code of this repo is based on the version [secantalpha/LpCVT](https://github.com/secantalpha/LpCVT).

If you use this code, you may need to cite the paper:

```
@article{10.1145/1778765.1778856,
    author = {L\'{e}vy, Bruno and Liu, Yang},
    title = {Lp Centroidal Voronoi Tessellation and its applications},
    year = {2010},
    issue_date = {July 2010},
    publisher = {Association for Computing Machinery},
    address = {New York, NY, USA},
    volume = {29},
    number = {4},
    issn = {0730-0301},
    url = {https://doi.org/10.1145/1778765.1778856},
    doi = {10.1145/1778765.1778856},
    month = jul,
    articleno = {119},
    numpages = {11},
    keywords = {Centroidal Voronoi Tessellation, anisotropic meshing, hex-dominant meshing, quad-dominant meshing}
}
```

# :bulb: What did I do?
- The data files have been reorganized to make this project lighter.
- Parts of the code have been rearranged and reorganized.

# :link: Requirements / dependencies

- [Boost](https://www.boost.org)
- [CGAL](https://www.cgal.org/index.html) (version >= 5.0)
- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)
- [LBFGSpp](https://github.com/yixuan/LBFGSpp) (already included)

# :checkered_flag: Getting started

```
git clone https://github.com/Bigger-and-Stronger/LpCVT.git
```

# :hammer: Build

In Windows

In macOS
```
mkdir build
cd build
cmake ..
make
```

# :computer: Usage

# :bar_chart: Output

# :children_crossing: Add to your code