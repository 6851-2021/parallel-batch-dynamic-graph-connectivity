# Parallel Batch Dynamic Graph Connectivity

This project shows a fast parallel implementation of batch-dynamic graph connectivity algorithm first described by [Acar, Anderson, Blelloch and Dhulipala](https://arxiv.org/abs/1903.08794). 
To simplify implementation, we make some modifications which increases the
depth bound from *O*(log<sup>3</sup>&nbsp;*n*) to
*O*(log<sup>4</sup>&nbsp;*n*).

The main item of interest in this project is `connectivity.h` which can be found in `include/batch_dynamic_connectivity/connectivity.h`. 

--------------------------------

## Backends
The codebase uses a couple of important data structures and libraries under the hood:
- `parallel_euler_tour_trees` for each level of our batch-dynamic HDT structure. They come from Tom T. Tseng's implementation of [batch-dynamic ETTs](https://github.com/tomtseng/parallel-euler-tour-tree). 
- `include/sequence` comes from the [pbbs](https://github.com/tomtseng/parallel-euler-tour-tree) implementation of sequences and is used heavily in `parallel_euler_tour_trees`.
- `include/utilities` is derived from [pbbs](https://github.com/tomtseng/parallel-euler-tour-tree) utilities
- `include/catch.hpp` is based on [Catch](https://github.com/catchorg/Catch2) which is a fast framework for testing
- `external/parlay` is [parlaylib](https://github.com/cmuparlay/parlaylib) and is used for doing essential parallel primitives throughout our code.

All of the parallelism is supported through [OpenCilk](https://cilk.mit.edu). We use the OpenCilk v1 compiler for most of this C++ compilation. The compiler can be changed easily within the `Makefile`.

---------------------------------
## Compilation

To compile, we recommend you install the OpenCilk compiler and substitute that into the `Makefile`:
```make
CXX = /path/to/opencilk/clang++
```
After that run
```sh
make CILK=1
```
Run ``make clean`` to delete the compiled binaries.

---------------------------------
## Resources
- Thomas Tseng, Laxman Dhulipala, and Guy Blelloch. [Batch-parallel Euler tour
trees](https://arxiv.org/abs/1810.10738). In _Proceedings of the Twenty-First
Workshop on Algorithm Engineering and Experiments_, page to appear.  Society for
Industrial and Applied Mathematics, 2019.
- Acar, Umut A., Daniel Anderson, Guy E. Blelloch, and Laxman Dhulipala. "[Parallel batch-dynamic graph connectivity](https://dl.acm.org/doi/abs/10.1145/3323165.3323196?casa_token=7Vwyr7HhOz0AAAAA:5_wxIRJTV9zc7qpo0AfXPcxfngngeR-6Tw6F3Ex4dz0L6yxOurmp1fcd_4LHyz7H989sTnXfLca5)." In _The 31st ACM Symposium on Parallelism in Algorithms and Architectures_, pp. 381-392. 2019.
