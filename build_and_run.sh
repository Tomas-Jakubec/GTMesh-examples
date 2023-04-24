#! /bin/bash
cd Meshes
python3 simple_mesh_generator.py
cd ..
mkdir build_dir
cmake -S . -B build_dir
cmake --build build_dir -DCMAKE_BUILD_TYPE=Release
mkdir -p out
./build_dir/GTMesh-example "Meshes\mesh3D.vtk" "out"
mkdir -p out-graph
./build_dir/GTMesh-example-omp_parallel_graph "Meshes\mesh3D.vtk" "out-graph"
mkdir -p out-struct
./build_dir/GTMesh-example-omp_parallel_struct "Meshes\mesh3D.vtk" "out-struct"
