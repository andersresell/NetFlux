# NetFlux

A 3D unstructured cell centered finite volume solver for compressible flow.

External depencencies: Eigen and Yaml-cpp

Plan forward:
    - Make working explicit solver for Euler eqs
    - Implement moving grid capability using Arbitrary Eulerian-Lagrangian (ALE) formulation
    - Parallellize with MPI