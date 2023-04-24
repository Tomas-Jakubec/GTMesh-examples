This repository contains several applications presenting the usage of the GTMesh library.
GTMesh library is a basic library for handling unstructured meshes and computation on them.

The computation is demonstrated on the heat conduction problem:
```math
\begin{array}{cll}
\frac{\partial T\left(\boldsymbol{x},t\right)}{\partial t} & =-\nabla T(\boldsymbol{x},t) & \text{ for }\boldsymbol{x}\in\varOmega^{\circ}\\
T\left(\boldsymbol{x},t\right) & =T_{\text{wall}} & \text{ for }\boldsymbol{x}\in\partial\varOmega
\end{array}
```

There are three solutions demonstrated. All use finite volume method and Merson version of the Runge-Kutta explicit solver.
- simple single threaded approach,
- multithreaded approach utilizing graph coloring to prevent race conditions,
- multithreaded approach utilizing storing of auxiliary quantities to prevent race conditions. 

# Dependencies
- C++ compiler with OpenMP support
- CMake
- Paraview or another tool to display output data
- python3 for generation of basic meshes

# How to run examples
1. generate the mesh
    ```bash
     cd Meshes
     python3 simple_mesh_generator.py
     cd ..
    ```
    This will produce `Meshes\mesh3D.vtk`. To produce a finer mesh, increase the number of cells on the line 57.
2. build the application
    ```bash
    mkdir build_dir
    cmake -S . -B build_dir
    cmake --build build_dir -DCMAKE_BUILD_TYPE=Release
    ```
3. execute the example 
    ```bash
    mkdir -p out
    ./build/GTMesh-example "Meshes\mesh3D.vtk" "out"
    mkdir -p out-graph
    ./build/GTMesh-example-omp_parallel_graph "Meshes\mesh3D.vtk" "out-graph"
    mkdir -p out-struct
    ./build/GTMesh-example-omp_parallel_struct "Meshes\mesh3D.vtk" "out-struct"
    ```
    the result is stored in the `out`
