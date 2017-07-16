var usermanual =
[
    [ "Introduction", "manual-introduction.html", null ],
    [ "Installation", "manual-installation.html", [
      [ "Dependencies", "manual-installation.html#manual-installation-dependencies", null ],
      [ "Generic Installation of ViennaCL", "manual-installation.html#manual-installation-generic", null ],
      [ "Get the OpenCL Library", "manual-installation.html#manual-installation-opencl", [
        [ "NVIDIA Driver", "manual-installation.html#manual-installation-opencl-nvidia", null ],
        [ "AMD Accelerated Parallel Processing SDK (formerly Stream SDK)", "manual-installation.html#manual-installation-opencl-amd", null ],
        [ "INTEL OpenCL SDK", "manual-installation.html#manual-installation-opencl-intel", null ]
      ] ],
      [ "Enabling OpenMP, OpenCL, or CUDA Backends", "manual-installation.html#manual-installation-backends", null ],
      [ "Building the Examples and Tutorials", "manual-installation.html#manual-installation-examples", [
        [ "Manual Builds (no CMake)", "manual-installation.html#manual-installation-examples-manual", null ],
        [ "CMake-Assisted Builds", "manual-installation.html#manual-installation-examples-cmake", null ],
        [ "Linux", "manual-installation.html#manual-installation-examples-linux", null ],
        [ "Mac OS X", "manual-installation.html#manual-installation-examples-macos", null ],
        [ "Windows", "manual-installation.html#manual-installation-examples-windows", null ]
      ] ]
    ] ],
    [ "Basic Types", "manual-types.html", [
      [ "Scalar Type", "manual-types.html#manual-types-scalar", [
        [ "Example Usage", "manual-types.html#manual-types-scalar-usage", null ],
        [ "Members", "manual-types.html#manual-types-scalar-members", null ]
      ] ],
      [ "Vector Type", "manual-types.html#manual-types-vector", [
        [ "Example Usage", "manual-types.html#manual-types-vector-usage", null ],
        [ "Members", "manual-types.html#manual-types-vector-members", null ]
      ] ],
      [ "Dense Matrix Type", "manual-types.html#manual-types-matrix", [
        [ "Example Usage", "manual-types.html#manual-types-matrix-usage", null ],
        [ "Members", "manual-types.html#manual-types-matrix-members", null ]
      ] ],
      [ "Initializer Types", "manual-types.html#manual-operations-initializers", null ],
      [ "Sparse Matrix Types", "manual-types.html#manual-types-sparse", [
        [ "Compressed Matrix", "manual-types.html#manual-types-sparse-compressed", [
          [ "Example Usage", "manual-types.html#manual-types-sparse-compressed-usage", null ],
          [ "Members", "manual-types.html#manual-types-sparse-compressed-members", null ]
        ] ],
        [ "Coordinate Matrix", "manual-types.html#manual-types-sparse-coordinate", [
          [ "Example Usage", "manual-types.html#manual-types-sparse-coordinate-usage", null ],
          [ "Members", "manual-types.html#manual-types-sparse-coordinate-members", null ]
        ] ],
        [ "ELL Matrix", "manual-types.html#manual-types-sparse-ell", null ],
        [ "Sliced ELL Matrix", "manual-types.html#manual-types-sparse-sliced-ell", null ],
        [ "Hybrid Matrix", "manual-types.html#manual-types-sparse-hyb", null ],
        [ "Compressed Compressed Matrix", "manual-types.html#manual-types-sparse-compressed-compressed", null ]
      ] ],
      [ "Proxies", "manual-types.html#manual-types-proxies", null ]
    ] ],
    [ "Basic Operations", "manual-operations.html", [
      [ "Elementary Vector and Matrix Operations (BLAS Level 1)", "manual-operations.html#manual-operations-blas1", null ],
      [ "Matrix-Vector Operations (BLAS Level 2)", "manual-operations.html#manual-operations-blas2", null ],
      [ "Matrix-Matrix Operations (BLAS Level 3)", "manual-operations.html#manual-operations-blas3", null ],
      [ "Sparse Matrix Operations", "manual-operations.html#manual-operations-sparse", [
        [ "Sparse Matrix-Vector Products", "manual-operations.html#manual-operations-sparse-spmv", null ],
        [ "Sparse Matrix-Matrix Products", "manual-operations.html#manual-operations-sparse-spgemm", null ]
      ] ],
      [ "Row, Column, and Diagonal Extraction", "manual-operations.html#manual-operations-row-column-diagonal", null ],
      [ "Conversion Operations", "manual-operations.html#manual-operations-conversion", null ]
    ] ],
    [ "Algorithms", "manual-algorithms.html", [
      [ "Direct Solvers", "manual-algorithms.html#manual-algorithms-direct-solvers", null ],
      [ "Iterative Solvers", "manual-algorithms.html#manual-algorithms-iterative-solvers", [
        [ "Conjugate Gradients (CG)", "manual-algorithms.html#manual-algorithms-iterative-solvers-cg", null ],
        [ "Mixed-Precision Conjugate Gradients", "manual-algorithms.html#manual-algorithms-iterative-solvers-mixed-cg", null ],
        [ "Stabilized Bi-CG (BiCGStab)", "manual-algorithms.html#manual-algorithms-iterative-solvers-bicgstab", null ],
        [ "Generalized Minimum Residual (GMRES)", "manual-algorithms.html#manual-algorithms-iterative-solvers-gmres", null ]
      ] ],
      [ "Preconditioners", "manual-algorithms.html#manual-algorithms-preconditioners", [
        [ "Parallel Incomplete LU Factorization with Static Pattern (Chow-Patel-ILU0)", "manual-algorithms.html#manual-algorithms-preconditioners-parallel-ilu0", null ],
        [ "Parallel Incomplete Cholesky Factorization with Static Pattern (Chow-Patel-IChol0)", "manual-algorithms.html#manual-algorithms-preconditioners-parallel-icc0", null ],
        [ "Incomplete LU Factorization with Threshold (ILUT)", "manual-algorithms.html#manual-algorithms-preconditioners-ilut", null ],
        [ "Incomplete LU Factorization with Static Pattern (ILU0)", "manual-algorithms.html#manual-algorithms-preconditioners-ilu0", null ],
        [ "Incomplete Cholesky Factorization with Static Pattern (IChol0)", "manual-algorithms.html#manual-algorithms-preconditioners-icc0", null ],
        [ "Block-ILU", "manual-algorithms.html#manual-algorithms-preconditioners-block-ilu", null ],
        [ "Jacobi Preconditioner", "manual-algorithms.html#manual-algorithms-preconditioners-jacobi", null ],
        [ "Row-Scaling Preconditioner", "manual-algorithms.html#manual-algorithms-preconditioners-row-scaling", null ],
        [ "Algebraic Multigrid Preconditioners", "manual-algorithms.html#manual-algorithms-preconditioners-amg", null ]
      ] ],
      [ "Eigenvalue Computations", "manual-algorithms.html#manual-algorithms-eigenvalues", [
        [ "Power Iteration", "manual-algorithms.html#manual-algorithms-eigenvalues-power", null ],
        [ "The Lanczos Algorithm", "manual-algorithms.html#manual-algorithms-eigenvalues-lanczos", null ]
      ] ],
      [ "QR Factorization", "manual-algorithms.html#manual-algorithms-qr-factorization", null ]
    ] ],
    [ "Interfacing Other Libraries", "manual-interfacing.html", [
      [ "Boost.uBLAS", "manual-interfacing.html#manual-interfacing-ublas", null ],
      [ "Armadillo", "manual-interfacing.html#manual-interfacing-armadillo", null ],
      [ "Eigen", "manual-interfacing.html#manual-interfacing-eigen", null ],
      [ "MTL 4", "manual-interfacing.html#manual-interfacing-mtl", null ]
    ] ],
    [ "Memory Model", "manual-memory.html", [
      [ "Memory Handle Operations", "manual-memory.html#manual-memory-handle", null ],
      [ "Querying and Switching Active Memory Domains", "manual-memory.html#manual-memory-switching", null ]
    ] ],
    [ "Shared Library", "manual-shared-library.html", null ],
    [ "Additional Algorithms (Unstable)", "manual-additional-algorithms.html", [
      [ "Additional Preconditioners", "manual-additional-algorithms.html#manual-additional-algorithms-preconditioners", [
        [ "Sparse Approximate Inverses", "manual-additional-algorithms.html#manual-additional-algorithms-preconditioners-spai", null ]
      ] ],
      [ "Additional Eigenvalue Routines", "manual-additional-algorithms.html#manual-additional-algorithms-eigenvalues", [
        [ "Symmetric Tridiagonal Matrices: Bisection", "manual-additional-algorithms.html#manual-additional-algorithms-eigenvalues-bisection", null ],
        [ "Symmetric Tridiagonal Matrices: TQL2", "manual-additional-algorithms.html#manual-additional-algorithms-eigenvalues-tql2", null ],
        [ "QR Method for Symmetric Dense Matrices", "manual-additional-algorithms.html#manual-additional-algorithms-eigenvalues-qr-method-symmetric", null ]
      ] ],
      [ "Fast Fourier Transform", "manual-additional-algorithms.html#manual-additional-algorithms-fft", null ],
      [ "Singular Value Decomposition", "manual-additional-algorithms.html#manual-additional-algorithms-svd", null ],
      [ "Bandwidth Reduction", "manual-additional-algorithms.html#manual-additional-algorithms-bandwidth-reduction", null ],
      [ "Nonnegative Matrix Factorization", "manual-additional-algorithms.html#manual-additional-algorithms-nmf", null ]
    ] ],
    [ "User-Provided OpenCL Contexts", "manual-custom-contexts.html", [
      [ "Passing Contexts to ViennaCL", "manual-custom-contexts.html#manual-custom-contexts-passing", null ],
      [ "Wrapping Existing Memory with ViennaCL Types", "manual-custom-contexts.html#manual-custom-contexts-wrapping", null ]
    ] ],
    [ "Configuring OpenCL Contexts and Devices", "manual-multi-device.html", [
      [ "Context Setup", "manual-multi-device.html#manual-multi-device-context", null ],
      [ "Switching Contexts and Devices", "manual-multi-device.html#manual-multi-device-switching", null ],
      [ "Setting OpenCL Compiler Flags", "manual-multi-device.html#manual-multi-device-compiler", null ]
    ] ],
    [ "Custom OpenCL Compute Kernels", "manual-custom-kernels.html", [
      [ "Setting up the OpenCL Source Code", "manual-custom-kernels.html#manual-custom-kernels-opencl-source", null ],
      [ "Compilation of the OpenCL Source Code", "manual-custom-kernels.html#manual-custom-kernels-opencl-build", null ],
      [ "Launching the OpenCL Kernel", "manual-custom-kernels.html#manual-custom-kernels-opencl-launch", null ]
    ] ],
    [ "Structured Matrix Types", "manual-structured-matrix.html", [
      [ "Circulant Matrix", "manual-structured-matrix.html#manual-structured-matrix-circulant", null ],
      [ "Hankel Matrix", "manual-structured-matrix.html#manual-structured-matrix-hankel", null ],
      [ "Toeplitz Matrix", "manual-structured-matrix.html#manual-structured-matrix-toeplitz", null ],
      [ "Vandermonde Matrix", "manual-structured-matrix.html#manual-structured-matrix-vandermonde", null ]
    ] ],
    [ "Design Decisions", "manual-design.html", [
      [ "Transfer CPU-GPU-CPU for Scalars", "manual-design.html#manual-design-transfer-scalars", null ],
      [ "Transfer CPU-GPU-CPU for Vectors", "manual-design.html#manual-design-transfer-vectors", null ],
      [ "Solver Interface", "manual-design.html#manual-design-solver", null ],
      [ "Iterators", "manual-design.html#manual-design-iterators", null ],
      [ "Initialization of Compute Kernels", "manual-design.html#manual-design-init", null ]
    ] ],
    [ "Authors and Contributors", "manual-authors.html", [
      [ "Project Head", "manual-authors.html#manual-authors-head", null ],
      [ "Authors", "manual-authors.html#manual-authors-authors", null ],
      [ "Additional Contributors", "manual-authors.html#manual-authors-contributors", null ]
    ] ],
    [ "Versioning", "manual-versioning.html", null ],
    [ "Change Log", "changelog.html", null ],
    [ "License", "manual-license.html", null ]
];