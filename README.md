## geoLR (Version 1.0.0)
geoLR is a package and framework to perform fast geometric low-rank approximation for dense possibly rectangular matrices associated with a smooth kernel function. geoLR takes two sets of points as input and does not require forming the kernel matrix. Though designed for general rectangular matrices corresponding to two sets of points, it can also be used for symmetric kernel matrices associated with one set of points.


## Main Features
* geoLR is generic, allowing user-defined kernel functions, both uniform and non-uniform data. 
* AnchorNet works for datasets in arbitrary dimensions.
* AnchorNet does not require forming the kernel matrix and does not assume the matrix to be symmetric or square. The algorithm computes low-rank factors given data and kernel function only.
* The total complexity of AnchorNet is O(rN) for computing a rank-r approximaiton to an N-by-N kernel matrix.

## Use
geoLR is in active development (currently 1.0.0) and its interface may change.

## Files
main.m: Driver.    
geoLR: Geometric Low-Rank Compression for Kernel Matrices with I/O below.

* Input: data sets 'X' and 'Y' (can be identical), low-rank approximation level 'pk', kernel function handle 'ff'.
* Output: low-rank factors U,V for approximating A with U*V.

## Reference
 -  [Data-Driven Linear Complexity Low-Rank Approximation of General Kernel Matrices: A Geometric Approach](https://arxiv.org/abs/2212.12674) Difeng Cai, Edmond Chow, Yuanzhe Xi
