# edgeDrawing (An implementation of the ED algorithm)

Version 0.9, 07/02/2026.

Future versions: <https://github.com/pmonasse/Edge-Drawing/tree/standalone>

The ED algorithm[^1] is a vectorized edge detector. It extends the Canny edge detector by replacing the hysteresis thresholding with edge tracking along level lines. 

This code is linked to an IPOL submission[^2].

## Build
*Requirements:*

  - CMake >= 3.13 <https://cmake.org/download/>
  - C++ compiler
  
*Build instructions:*
- Unix, MacOS:
  ```
  $ cd /path_to_this_file/
  $ cmake -DCMAKE_BUILD_TYPE:bool=Release -B build
  $ cmake --build build
  ```
- Windows with MinGW:
  ```
  $ cd /path_to_this_file/
  $ cmake -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE:bool=Release -B build
  $ cmake --build build
  ```

## Usage
```
Usage: ./build/edgeDrawing [options] in.png [out.png]
-g, --grad-min=ARG Min gradient (6)
-a, --anchor-gap=ARG Min gap of gradient for anchor (2)
-l, --length-min=ARG Min length of edge segment (10)
-s, --sigma=ARG Sigma of Gaussian blur (1)
-t, --t-junctions=ARG T-junctions output image
-e, --epsNFA=ARG log10(NFA) for validation, normally 0 or negative (0)
-S, --segLevel=ARG Segmentation level for valid edges (0)
out.png may be omitted if -t or --t-junctions was used
NFA validation only with -e and/or -S
segLevel:
 0: validate only full edges
 1: validate full edge if a segment is valid
 other: validate only segments
```

Typical settings for a contrario validation of edges:

- a smaller value of -g could be used (2).
- a smaller value of -a could be used (0).
- the NFA threshold -e could be 0.

Example:
```
build/edgeDrawing data/shapes.png shapes_out.png
```
Compare `shapes_out.png` and reference `data/shapes_out_png`.

## Files
Reviewed for IPOL publication:
chain_tree.h chain_tree.cpp edge_drawing.h edge_drawing.cpp main.cpp
maxTree.h validation.cpp

Support files:
CMakeLists.txt cmdLine.h image.h image.cpp io_png.h io_png.c README.md

Other:
third_party/zlib-xxx third_party/libpng-yyy data/zzz.png

[^1]: Original implementation by the authors of ED, Cihan Topal and Cuneyt Akinlar (https://github.com/CihanTopal/ED_Lib)

[^2] Edge Drawing: A Fast Edge Segment Detector, Adle Ben Salem and Pascal Monasse, IPOL, 2026
