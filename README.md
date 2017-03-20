# Blended barycentric coordinates Ver. 1.0.0

### Description

This set of C++ classes provides you with an implementation of blended barycentric coordinates from the paper: D. Anisimov, D. Panozzo, and K. Hormann. Blended barycentric coordinates. Computer Aided Geometric Design, to appear, 2017. The implementation is based on the pseudocodes from Appendix A of the paper that can be found [here](http://www.inf.usi.ch/hormann/papers/Anisimov.2017.BBC.pdf).

##### NOTE: This code has been tested only on Mac OS!

### Run the code

In order to run the code, do the following:

1. Install [macports](https://www.macports.org/install.php)
2. Open terminal and type the following:

```bash
  sudo port install cmake
```
```bash
  cd path_to_the_folder/blmv/2d/
```
```bash
  mkdir bin
```
```bash
  cd bin
```
```bash
  cmake -DCMAKE_BUILD_TYPE=Debug ..
```
```bash
  make
```
```bash
  ./blmv
```

For the release version use instead: 

```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
```

### Example

```C++
// Polygon.
std::vector<VertexR2> poly(7);

poly[0] = VertexR2(-0.542, -0.740);
poly[1] = VertexR2(-0.066, -0.740);
poly[2] = VertexR2(   0.0, -0.086);
poly[3] = VertexR2( 0.066, -0.740);
poly[4] = VertexR2( 0.542, -0.740);
poly[5] = VertexR2( 0.406,  0.444);
poly[6] = VertexR2(-0.406,  0.444);

// Evaluation point.
VertexR2 query(0.0, 0.2);

// Storage for the computed blended coordinates.
std::vector<double> b;

// Compute blended coordinates.
BlendedR2 blc(poly);

blc.setContinuity(2);
blc.compute(query, b);

// Output the resulting coordinates.
std::cout << "\nResult: ";
for (size_t i = 0; i < b.size(); ++i) std::cout << b[i] << " ";
std::cout << "\n\n";
```

##### NOTE: For the complete example see main.cpp!

### Bugs

If you find any bugs, please report them to me, and I will try to fix them as soon as possible! Please also note that this code does not provide the same timings as in the paper for the two reasons.
First, the triangle search here is implemented brute-force and second, the O(n) linear boundary behaviour is a part of the computation, where n is the number of the polygon's vertices (see comments in the code for more details).
