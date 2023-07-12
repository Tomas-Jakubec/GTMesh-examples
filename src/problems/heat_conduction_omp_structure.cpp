#include <memory>
#include <limits>
#include "solvers/RKMSolver_omp.hpp"

#include <GTMesh/Traits/Traits.h>
#include <GTMesh/Traits/TraitsAlgorithm/TraitsAlgorithm.h>
#include <GTMesh/UnstructuredMesh/UnstructuredMesh.h>
#include <GTMesh/UnstructuredMesh/MeshDataContainer/MeshDataIO/VTKMeshDataWriter.h>


struct ComputationData
{
    double T; //!< temperature
};

MAKE_NAMED_ATTRIBUTE_TRAIT(ComputationData, "temperature", T);

struct FaceData
{
    double measureOverCellsDistance;
    double measure;
    // Auxiliary structure to store the intermediate result in the data structure of the edge
    ComputationData flux;
};

struct CellData
{
    double invCellVolume;
};


template<unsigned int ProblemDimension>
struct HeatConductionProblem{
    using MeshType = UnstructuredMesh<ProblemDimension, size_t, double>;
    using ResultType = ComputationData;
    using ProblemDataContainerType = MeshDataContainer<ResultType, ProblemDimension>;
    std::shared_ptr<MeshReader<ProblemDimension>> meshReader;
    VTKMeshWriter<ProblemDimension, size_t, double> meshWriter;

    MeshType mesh;
    const double T_wall = 300;
    MeshDataContainer<std::tuple<CellData, FaceData>,ProblemDimension, ProblemDimension-1> meshData;
    std::vector<std::vector<std::size_t>> nonConcurentFaces;


    void calculateRHS(double time, //time is unused in this problem
                      const ProblemDataContainerType &compData,
                      ProblemDataContainerType &outDeltas){
        #pragma omp for
        for (std::size_t i = 0; i < mesh.getFaces().size(); ++i){
            const auto& face = mesh.getFaces()[i];
            const auto cRI = face.getCellRightIndex(), cLI = face.getCellLeftIndex();
            if (!isBoundaryIndex(cRI) and !isBoundaryIndex(cLI)){
                const auto &cR = mesh.getCells().at(cRI), &cL = mesh.getCells()[cLI];
                meshData[face].flux.T = (compData.at(cL).T - compData.at(cR).T) * meshData[face].measureOverCellsDistance;
            } else if (isBoundaryIndex(cLI)) {
                const auto &cR = mesh.getCells().at(cRI);
                meshData[face].flux.T = (T_wall - compData[cR].T) * meshData[face].measureOverCellsDistance;
            } else {
                const auto &cL = mesh.getCells().at(cLI);
                meshData[face].flux.T = (compData[cL].T - T_wall) * meshData[face].measureOverCellsDistance;
            }
        }

        #pragma omp for
        for (std::size_t cellI = 0; cellI < mesh.getCells().size(); ++cellI) {
            const auto& cell = mesh.getCells()[cellI];
            ComputationData auxDelta {0};
            mesh.template apply<ProblemDimension, ProblemDimension - 1>(cellI,
                [&] (std::size_t cellIndex, std::size_t faceIndex){
                    const auto& face = mesh.getFaces()[faceIndex];
                    if (cellIndex == face.getCellLeftIndex()) {
                        auxDelta.T -= meshData[face].flux.T;
                    } else {
                        auxDelta.T += meshData[face].flux.T;
                    }
                });
            outDeltas[cell] = auxDelta;
        }
    }

    ProblemDataContainerType loadMesh(const std::string& meshPath){
        meshReader = mesh.load(meshPath);
        mesh.initializeCenters();
        mesh.setupBoundaryCells();
        mesh.setupBoundaryCellsCenters();
        meshData.allocateData(mesh);
        auto measures = mesh.template computeElementMeasures<METHOD_TESSELLATED>();
        auto cellsDist = computeCellsDistance(mesh);
        for (const auto& cell : mesh.getCells()){
            meshData[cell].invCellVolume = 1.0 / measures[cell];
        }
        for (const auto& face : mesh.getFaces()){
            meshData[face].measure = measures[face];
            meshData[face].measureOverCellsDistance = measures[face] / cellsDist[face];
        }
        return ProblemDataContainerType(mesh, ComputationData{273.15});
    }

    void exportMeshAndData(ProblemDataContainerType& compData,
                           const std::string& outputPath){
        std::ofstream ofst(outputPath);
        meshWriter.writeHeader(ofst, "heat-conduction");
        meshWriter.writeToStream(ofst, mesh,meshReader->getCellTypes());
        VTKMeshDataWriter<ProblemDimension> dataWriter;
        dataWriter.writeToStream(ofst, compData, meshWriter);
    }
};

int main(int argc, char** argv) {
    HeatConductionProblem<3> hcp;       // 3D version
//     HeatConductionProblem<2> hcp;    // 2D version
    std::string defaultMeshPath = "../Meshes/mesh3D.vtk";       // 3D version
//     std::string defaultMeshPath = "../Meshes/mesh2D.vtk";    // 2D version
    std::string meshPath = argc <= 1 ? defaultMeshPath : argv[1];
    auto compData = hcp.loadMesh(meshPath);
    std::string outPath = argc <= 2 ? "out_omp_struct" : argv[2];
    hcp.exportMeshAndData(compData, outPath + "/heat_conduction-omp_struct-t_0s.vtk");
    for (int i = 0; i < 10; ++i) {
        RKMSolverOMP(hcp, compData, 1e-3, i, i + 1.0, 1e-4);
        hcp.exportMeshAndData(compData, outPath + "/heat_conduction-omp_struct-t_" + std::to_string(i+1) + "s.vtk");
    }
}
