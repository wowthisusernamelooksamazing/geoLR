## geoLR (Version 1.0.0)
geoLR is a package and framework to perform fast geometric low-rank approximation for dense possibly rectangular matrices associated with a smooth kernel function. geoLR takes two sets of points as input and does not require forming the kernel matrix. Though designed for general rectangular matrices corresponding to two sets of points, it can also be used for symmetric kernel matrices associated with one set of points.


## Main Features
* geoLR is generic, allowing user-defined kernel functions, both uniform and non-uniform data. 
* AnchorNet works for datasets in arbitrary dimensions.
* AnchorNet does not require forming the kernel matrix and does not assume the matrix to be symmetric or square. The algorithm computes low-rank factors given data and kernel function only.
* The total complexity of AnchorNet is O(rN) for computing a rank-r approximaiton to an N-by-N kernel matrix.

## Use
geoLR is in active development (currently 1.0.0) and its interface may change.

## Drivers
main.m

## Reference
 -  [Data-Driven Linear Complexity Low-Rank Approximation of General Kernel Matrices: A Geometric Approach](https://arxiv.org/abs/2212.12674) Difeng Cai, Edmond Chow, Yuanzhe Xi


## Versioning
geoLR attempts to follow [semantic versioning](https://www.semver.org). Do note, that in it's current (1.0.0) development, such versioning may not be strictly
followed.