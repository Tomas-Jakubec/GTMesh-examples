#include <memory>
#include <limits>
#include "RKMSolver.hpp"

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

    MeshType mesh;
    const double T_wall = 300;
    MeshDataContainer<std::tuple<CellData, FaceData>,ProblemDimension, ProblemDimension-1> meshData;
    std::vector<std::vector<std::size_t>> nonConcurentFaces;


    void calculateRHS(double time, //time is unused in this problem
                      const ProblemDataContainerType &compData,
                      ProblemDataContainerType &outDeltas){
#pragma omp parallel for
        for (const auto& cell : mesh.getCells()){
            outDeltas[cell].T = 0;
        }
        for (const auto& faces : nonConcurentFaces){
#pragma omp parallel for
            for (std::size_t i = 0; i < faces.size(); ++i){
                const auto& face = mesh.getFaces()[faces[i]];
                const auto cRI = face.getCellRightIndex(), cLI = face.getCellLeftIndex();
                if (!isBoundaryIndex(cRI) and !isBoundaryIndex(cLI)){
                    const auto &cR = mesh.getCells().at(cRI), &cL = mesh.getCells()[cLI];
                    const auto dT_dn = (compData.at(cL).T - compData.at(cR).T) * meshData[face].measureOverCellsDistance;
                    outDeltas.at(cL).T -= dT_dn;
                    outDeltas.at(cR).T += dT_dn;
                } else if (isBoundaryIndex(cLI)) {
                    const auto &cR = mesh.getCells().at(cRI);
                    const auto dT_dn = (T_wall - compData[cR].T) * meshData[face].measureOverCellsDistance;
                    outDeltas.at(cR).T += dT_dn;
                } else {
                    const auto &cL = mesh.getCells().at(cLI);
                    const auto dT_dn = (compData[cL].T - T_wall) * meshData[face].measureOverCellsDistance;
                    outDeltas.at(cL).T -= dT_dn;
                }
            }
        }
    }

    void meshFaceColoring(){
        // Randomized coloring creates more even distribution of colors in unstrucutred meshes.
        DBGCHECK;
        auto faceColoring = mesh.template coloring<ProblemDimension-1, ProblemDimension, ColoringMethod::METHOD_RANDOM>();DBGCHECK;
        std::size_t colorCount =0;DBGCHECK;
        for (auto colorIndex : faceColoring.template getDataByPos<0>()) {
            if (colorCount < colorIndex) {
                colorCount = colorIndex+1;
            }
        }DBGCHECK;
        nonConcurentFaces.resize(colorCount);DBGCHECK;
        for (const auto& face : mesh.getFaces()){
            unsigned int colorIndex = faceColoring[face];
            nonConcurentFaces[colorIndex].emplace_back(face.getIndex());
        }
        DBGCHECK;
    }


    ProblemDataContainerType loadMesh(const std::string& meshPath){
        meshReader = mesh.load(meshPath);
        mesh.initializeCenters();
        mesh.setupBoundaryCells();
        mesh.setupBoundaryCellsCenters();
        meshData.allocateData(mesh);
        auto measures = mesh.computeElementMeasures();
        auto cellsDist = computeCellsDistance(mesh);
        for (const auto& cell : mesh.getCells()){
            meshData[cell].invCellVolume = 1.0 / measures[cell];
        }
        for (const auto& face : mesh.getFaces()){
            meshData[face].measure = measures[face];
            meshData[face].measureOverCellsDistance = measures[face] / cellsDist[face];
        }
        meshFaceColoring();
        return ProblemDataContainerType(mesh, ComputationData{273.15});
    }

    void exportMeshAndData(ProblemDataContainerType& compData,
                           const std::string& outputPath){
        std::ofstream ofst(outputPath);
        VTKMeshWriter<ProblemDimension, size_t, double> meshWriter;
        meshWriter.writeHeader(ofst, "heat-conduction");
        meshWriter.writeToStream(ofst, mesh,meshReader->getCellTypes());
        VTKMeshDataWriter<ProblemDimension> dataWriter;
        dataWriter.writeToStream(ofst, compData, meshWriter);
    }
};

int main() {
    constexpr unsigned int Dim = 3;
    HeatConductionProblem<3> hcp;
    auto compData = hcp.loadMesh("../Meshes/mesh3D.vtk");
    hcp.exportMeshAndData(compData, "../out/heat_conduction_t_0s.vtk");
    for (int i = 0; i < 10; ++i) {
        RKMSolver(hcp, compData, 1e-3, i, i + 1.0, 1e-4);
        hcp.exportMeshAndData(compData, "../out/heat_conduction_t_" + std::to_string(i+1) + "s.vtk");
    }
}
