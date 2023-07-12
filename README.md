This repository contains several applications presenting the usage of the GTMesh library.
GTMesh library is a basic library for handling unstructured meshes and computation on them.

The computation is demonstrated on the heat conduction problem:
```math
\begin{array}{cll}
\frac{\partial T\left(\boldsymbol{x},t\right)}{\partial t} & =\Delta T(\boldsymbol{x},t) & \text{ for }\boldsymbol{x}\in\varOmega^{\circ}\\
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
1. fetch submodules `git submodule update --init`
1. generate the mesh
    ```bash
     cd Meshes
     python3 simple_mesh_generator.py
     cd ..
    ```
    This will produce `Meshes\mesh3D.vtk`. To produce a finer mesh, increase the number of cells in the main function.
2. build the application
    ```bash
    mkdir build_dir
    cmake -S . -B build_dir
    cmake --build build_dir -DCMAKE_BUILD_TYPE=Release
    ```
3. execute the example 
    ```bash
    mkdir -p out
    ./build_dir/GTMesh-example "Meshes\mesh3D.vtk" "out"
    mkdir -p out-graph
    ./build_dir/GTMesh-example-omp_parallel_graph "Meshes\mesh3D.vtk" "out-graph"
    mkdir -p out-struct
    ./build_dir/GTMesh-example-omp_parallel_struct "Meshes\mesh3D.vtk" "out-struct"
    ```
    the result is stored in the `out`, `out-graph` and `out-struct` directories.

# Unstructured mesh
There is `mesh3D-unstructured.fpma` file in the `Meshes` directory.
It contains an example of an unstructured mesh.
Try running
```bash
mkdir -p out
./build_dir/GTMesh-example "Meshes\mesh3D-unstructured.fpma" "out"
mkdir -p out-graph
./build_dir/GTMesh-example-omp_parallel_graph "Meshes\mesh3D-unstructured.fpma" "out-graph"
mkdir -p out-struct
./build_dir/GTMesh-example-omp_parallel_struct "Meshes\mesh3D-unstructured.fpma" "out-struct"
```
The results are still exported in VTK format, however, the cells are tessellated to tetrahedrons.
Use the surface with edges to see the tessellation. 

> Note that the computation on the unstructured mesh is more time-consuming.
> The causes might be as follows:
> - The number of cells and faces increases.
> - The regular mesh geometry might be more suitable for the current problem.

# Running in 2D
To instruct the algorithm to compute in 2D, uncomment the lines that commented with `2D version` and comment out the ones with `3D version`. Recompile the project and try running.
```bash
mkdir -p out-2D
./build_dir/GTMesh-example "Meshes\mesh2D.vtk" "out-2D"
mkdir -p out-2D-graph
./build_dir/GTMesh-example-omp_parallel_graph "Meshes\mesh2D.vtk" "out-2D-graph"
mkdir -p out-2D-struct
./build_dir/GTMesh-example-omp_parallel_struct "Meshes\mesh2D.vtk" "out-2D-struct"
```