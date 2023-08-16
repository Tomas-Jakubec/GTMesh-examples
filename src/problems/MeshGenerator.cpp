#include <memory>
#include <limits>
#include <iomanip>
#include <assert.h>
#include "solvers/RKMSolver.hpp"

#include <GTMesh/Traits/Traits.h>
#include <GTMesh/Traits/TraitsAlgorithm/TraitsAlgorithm.h>
#include <GTMesh/UnstructuredMesh/UnstructuredMesh.h>
#include <GTMesh/UnstructuredMesh/MeshDataContainer/MeshDataIO/VTKMeshDataWriter.h>


struct ComputationData
{
    double f; //!< temperature
    size_t cellIndex;
};

MAKE_NAMED_ATTRIBUTE_TRAIT(ComputationData, "F(x)", f, "cellIndex", cellIndex);


template<size_t Dimension, typename IndexType = size_t, typename RealType = double>
class VirtualGrid {

    using VertexType = Vertex<Dimension, RealType>;
    const VertexType origin;
    const RealType diameter;
public:

    const IndexType nElements;
    VirtualGrid(const VertexType& origin, IndexType nElements, RealType diameter): origin(origin), nElements(nElements), diameter(diameter) {}

    static VirtualGrid witCenterOrigin(const VertexType& origin, IndexType nElements, RealType diameter) {
        VertexType shift;
        for(auto & coordinate: shift) {
            coordinate = diameter;
        }
        VertexType originShifted = origin - shift;
        return VirtualGrid(originShifted, nElements, diameter);
    }

    double stepSize() const {
        return (2 * diameter / nElements);
    }

    VertexType at(const std::array<IndexType, Dimension>& indices) const {
        VertexType result = origin;
        for(unsigned int i = 0; i < Dimension; ++i) {
            assert(indices[i] <= nElements);
            result[i] += stepSize() * indices[i];
        }
        return result;
    }
};

template<size_t Dimension, typename IndexType = size_t, typename RealType = double>
class VertexScoreCalculator {

public:
    using VertexType = Vertex<Dimension, RealType>;
    const std::vector<VertexType> vertices;

    const RealType scaleFactor;
    const RealType treshold = 0.1;
    const RealType valueTreshold = 0.25;

    VertexScoreCalculator(
            const std::vector<VertexType> vertices,
            RealType scaleFactor = 1.0,
            RealType treshold = 0.1,
            RealType valueTreshold = 0.25
            ): vertices(vertices), scaleFactor(scaleFactor), treshold(treshold), valueTreshold(valueTreshold){}

    double scoreVertex(const VertexType& xCell, const VertexType& x) const {
//        return std::pow(2.7, (xCell - x).normEuclid() - 1);
//        return 1.0 / std::pow(scaleFactor * (xCell - x).normEuclid() + 1, 2);
        return 1.0 / ( scaleFactor * (xCell - x).normEuclid() + 1);
    }

    struct VertexScoreData {
        RealType vertexScore;
        std::vector<IndexType> neighboringCells;
    };

    VertexScoreData calculateScore(const VertexType& x) const {
        VertexScoreData result = {0.0, std::vector<size_t>()};
        result.neighboringCells.reserve(5);
        std::vector<double> f_xs(vertices.size());
        double maxF = 0;
        for (size_t i = 0; i < vertices.size(); ++i) {
            const auto& vert = vertices[i];

            auto f_x = scoreVertex(x, vert);
            f_xs[i] = f_x;

            if(f_x > maxF) {
                maxF = f_x;
            }
        }
        for (size_t i = 0; i < vertices.size(); ++i){
            auto fraction = f_xs[i] / maxF;
            if (/*f_xs[i] > valueTreshold && */fraction > treshold) {
//                result.vertexScore = maxF;
                result.vertexScore -= std::pow(fraction, 4);
//                result.vertexScore -= fraction;
                result.neighboringCells.emplace_back(i);
            }
        }
        return result;
    }

    RealType calculateScoreSum(const VertexType& x) const {
        RealType result = 0.0;
        std::vector<double> f_xs(vertices.size());
        double maxF = 0;
        double minF = 1e8;
        for (size_t i = 0; i < vertices.size(); ++i) {
            const auto& vert = vertices[i];

            auto f_x = scoreVertex(x, vert);
            f_xs[i] = f_x;

            if(f_x > maxF) {
                maxF = f_x;
            }
            if(f_x < minF) {
                minF = f_x;
            }
        }
        return maxF - minF;
        for (size_t i = 0; i < vertices.size(); ++i){
//            result = - maxF;
//            result -= f_xs[i] / maxF;
//              result -= std::pow(1+f_xs[i] / maxF, 4);
//            result -= std::pow(1.0 + f_xs[i] / maxF, 2);
//            result += maxF - f_xs[i];
        }
//        return result;
    }
};

MAKE_ATTRIBUTE_TRAIT(VertexScoreCalculator<2>::VertexScoreData,
                             vertexScore,neighboringCells );


template<unsigned int ProblemDimension>
struct MeshGenerationProblem{
    using MeshType = UnstructuredMesh<ProblemDimension, size_t, double>;
    using ResultType = ComputationData;
    using ProblemDataContainerType = MeshDataContainer<ResultType, ProblemDimension>;
    std::shared_ptr<MeshReader<ProblemDimension>> meshReader;
    VTKMeshWriter<ProblemDimension, size_t, double> meshWriter;

    MeshType mesh;

    void calculateFx(ProblemDataContainerType& data, const std::vector<Vertex<ProblemDimension, double>> vertices) {
        std::vector<double> f_xs(vertices.size());
        VertexScoreCalculator<ProblemDimension> calculator(vertices);
        for (const auto& cell : mesh.getCells()) {
            auto score = calculator.calculateScore(cell.getCenter(), 5.0, 0.5);
            data[cell].f = score.vertexScore;
        }
    }

    ProblemDataContainerType loadMesh(const std::string& meshPath){
        meshReader = mesh.load(meshPath);
        mesh.initializeCenters();
        mesh.setupBoundaryCells();
        mesh.setupBoundaryCellsCenters();
        // meshData.allocateData(mesh);
        auto measures = mesh.template computeElementMeasures<METHOD_TESSELLATED>();
        auto cellsDist = computeCellsDistance(mesh);
//        for (const auto& cell : mesh.getCells()){
//            meshData[cell].invCellVolume = 1.0 / measures[cell];
//        }
//        for (const auto& face : mesh.getFaces()){
//            meshData[face].measure = measures[face];
//            meshData[face].measureOverCellsDistance = measures[face] / cellsDist[face];
//        }
        return ProblemDataContainerType(mesh, ComputationData{0,0});
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

template<typename F>
Vertex<2> steepestGradient(Vertex<2> x, F&& f, double maxDist = 0.01) {
    double delta = 1e-5;
    double invDelta = 1.0 / (2*delta);
    auto tmpX = x;
    double fx = f(x);
    int cnt = 0;
    while(true){
        double df_dx = (f(tmpX-Vertex<2>({-delta, 0})) - f(tmpX-Vertex<2>({delta, 0}))) * invDelta;
        double df_dy = (f(tmpX-Vertex<2>({0, -delta})) - f(tmpX-Vertex<2>({0, delta}))) * invDelta;
        Vector<2> df = {-df_dx, -df_dy};
        if (df.normEuclid() > maxDist) {
            df *= maxDist / df.normEuclid();
        }
        if (cnt > 2000 || df.normEuclid() < 1e-5) {
            DBGVAR(df.normEuclid(), cnt);
            return tmpX;
        }
        double gamma = 1.0;//maxDist / (df).normEuclid();
        double auxFx = 0;
        Vertex<2> auxX;
        // line search
        do {
            auxX = tmpX + gamma * df;
            auxFx = f(auxX);
            gamma *= 0.7;
        } while (auxFx > fx);
        cnt++;
        fx = auxFx;
        tmpX = auxX;
    }
}


Vector<2> multiply(const Vector<2,Vector<2>>& A, const Vector<2>& x){
    Vector<2> result{};
    for(unsigned int i = 0; i < 2; ++i) {
        for(unsigned int j = 0; j < 2; ++j) {
            result[i] += A[i][j] * x[j];
        }
    }
    return result;
}

double multiply(const Vector<2>& x, const Vector<2>& y){
    double result{};
    for(unsigned int i = 0; i < 2; ++i) {
        result += x[i] * y[i];
    }
    return result;
}

Vector<2,Vector<2>> tenzorMultiply(const Vector<2>& x, const Vector<2>& y_t) {
    Vector<2,Vector<2>> result{{},{}};
    for(unsigned int i = 0; i < 2; ++i) {
        for(unsigned int j = 0; j < 2; ++j) {
            result[i][j] += x[i] * y_t[j];
        }
    }
    return result;
}

Vector<2,Vector<2>> operator*(double alpha, const Vector<2,Vector<2>>& A) {
    auto result = A;
    for(unsigned int i = 0; i < 2; ++i) {
        for(unsigned int j = 0; j < 2; ++j) {
            result[i][j] *= alpha;
        }
    }
    return result;
}

template<typename F>
Vertex<2> optimizationDFP(Vertex<2> x, F&& f, double maxDist = 0.01) {
    double delta = 1e-5;
    double invDelta = 1.0 / (2*delta);

    auto x_i = x;

    double fx = f(x);
    double df_dx = (f(x_i-Vertex<2>({-delta, 0})) - f(x_i-Vertex<2>({delta, 0}))) * invDelta;
    double df_dy = (f(x_i-Vertex<2>({0, -delta})) - f(x_i-Vertex<2>({0, delta}))) * invDelta;
    Vector<2> g_i = {df_dx, df_dy};
    DBGVAR(g_i, g_i.normEuclid());
    Vector<2,Vector<2>> S {{1,0},{0,1}};
    int cnt = 0;
    while(true){
        if(++cnt > 100 || g_i.normEuclid() < 10e-5) {
            DBGVAR("success", cnt);
            return x_i;
        }
        auto d = -1.0 * multiply(S,  g_i);
        DBGVAR(S, d, g_i);

        double gamma = 1.0; //maxDist / d.normEuclid();
        double auxFx = 0;
        Vertex<2> x_i1;
        // line search
        do {
            x_i1 = x_i + gamma * d;
            auxFx = f(x_i1);
            DBGVAR_HTML(x_i1, auxFx);
            gamma *= 0.7;
        } while (auxFx > fx);
        fx = auxFx;
        DBGVAR_HTML(x_i1, auxFx);
        double df_dx = (f(x_i1-Vertex<2>({-delta, 0})) - f(x_i1-Vertex<2>({delta, 0}))) * invDelta;
        double df_dy = (f(x_i1-Vertex<2>({0, -delta})) - f(x_i1-Vertex<2>({0, delta}))) * invDelta;
        auto g_i1 = Vector<2>{df_dx, df_dy};
        auto p = x_i1 - x_i;
        auto q = g_i1 - g_i;
        DBGVAR(g_i, g_i1, x_i, x_i1);
        x_i = x_i1;
        g_i = g_i1;
        auto p_q = multiply(p,q);
        DBGVAR(p_q);
        if (abs(p_q) < 1e-30) {
            DBGMSG("crash", cnt);
            return x_i;
        }
        auto S_q = multiply(S,q);
        S = S - tenzorMultiply(S_q, S_q) / multiply(q, S_q) - tenzorMultiply(p,p) / multiply(p, p);
    }
}

template<typename F>
Vertex<2> optimizationBFGS(Vertex<2> x, F&& f, double maxDist = 0.01) {
    double delta = 1e-8;
    double invDelta = 1.0 / (2*delta);

    auto x_i = x;

    double fx = f(x);
    double df_dx = (f(x_i-Vertex<2>({-delta, 0})) - f(x_i-Vertex<2>({delta, 0}))) * invDelta;
    double df_dy = (f(x_i-Vertex<2>({0, -delta})) - f(x_i-Vertex<2>({0, delta}))) * invDelta;
    Vector<2> g_i = {df_dx, df_dy};
    DBGVAR(g_i, g_i.normEuclid());
    Vector<2,Vector<2>> S {{1,0},{0,1}};
    int cnt = 0;
    while(true){
        if(++cnt > 100 || g_i.normEuclid() < 10e-5) {
            DBGVAR("success", cnt);
            return x_i;
        }
        auto d = -1.0 * multiply(S,  g_i);
        DBGVAR(S, d, g_i);

        double gamma = maxDist / d.normEuclid();
        double auxFx = 0;
        Vertex<2> x_i1;
        // line search
        do {
            x_i1 = x_i + gamma * d;
            auxFx = f(x_i1);
            DBGVAR_HTML(x_i1, auxFx);
            gamma *= 0.7;
        } while (auxFx > fx);
        fx = auxFx;
        DBGVAR_HTML(x_i1, auxFx);
        double df_dx = (f(x_i1-Vertex<2>({-delta, 0})) - f(x_i1-Vertex<2>({delta, 0}))) * invDelta;
        double df_dy = (f(x_i1-Vertex<2>({0, -delta})) - f(x_i1-Vertex<2>({0, delta}))) * invDelta;
        auto g_i1 = Vector<2>{df_dx, df_dy};
        auto p = x_i1 - x_i;
        auto q = g_i1 - g_i;
        DBGVAR(g_i, g_i1, x_i, x_i1);
        x_i = x_i1;
        g_i = g_i1;
        auto p_q = multiply(p,q);
        DBGVAR(p, q, p_q);
        if (abs(p_q) < 1e-30) {
            DBGMSG("crash", cnt);
            return x_i;
        }
        auto S_q = multiply(S,q);
        auto deltaS = (tenzorMultiply(p, S_q) + tenzorMultiply(S_q, p) - (1 + (multiply(q, S_q)/p_q)) * tenzorMultiply(p,p)) / p_q;
        DBGVAR(deltaS);
        S -= deltaS;
    }
}



template<typename T>
std::vector<T> slice(const std::vector<T>& vec, const std::vector<size_t>& indices) {
    std::vector<T> result;
    result.reserve(indices.size());
    for (auto index : indices){
        result.emplace_back(vec[index]);
    }
    return result;
}

void preciseScan2(const VirtualGrid<2>& grid, const std::vector<Vertex<2>>& vertices) {
    VertexScoreCalculator<2> calculator(vertices);
    double minFx = 1000;
    Vertex<2> argMin = {0.0, 0.0};
    for (size_t i = 0; i < grid.nElements +1; ++i){
        for (size_t j = 0; j < grid.nElements +1; ++j){
            auto score = calculator.calculateScoreSum(grid.at({i,j}));
            std::cout << "X:\t" << std::setprecision(8) << grid.at({i,j})[0]
                      << "\tY:\t" << std::setprecision(8) << grid.at({i,j})[1]
                      << "\tZ:\t" << std::setprecision(8) << score << std::endl;
            if (minFx > score) {
                minFx = score;
                argMin = grid.at({i,j});
            }
        }
    }
    DBGVAR(minFx, argMin);
//    try {
//        steepestGradient(argMin, calculator, -double(vertices.size()), 2*grid.stepSize());
//    } catch (...) {
//        DBGMSG("continue");
//    }
}

void preciseScan(const VirtualGrid<2>& grid, const std::vector<Vertex<2>>& vertices) {
    std::vector<std::vector<VertexScoreCalculator<2>::VertexScoreData>> gridData;
    VertexScoreCalculator<2> calculator(vertices, 5.0, 0.9);
    for (size_t i = 0; i < grid.nElements +1; ++i){
        std::vector<VertexScoreCalculator<2>::VertexScoreData> row;
        gridData.emplace_back(row);
        for (size_t j = 0; j < grid.nElements +1; ++j){
            auto score = calculator.calculateScore(grid.at({i,j}));
            gridData.back().emplace_back(score);
        }
    }
    for (size_t i = 1; i < grid.nElements; ++i){
        for (size_t j = 1; j < grid.nElements; ++j){
            auto score = gridData[i][j].vertexScore;
            if (score < -2 &&
                score < gridData[i-1][j].vertexScore &&
                score < gridData[i+1][j].vertexScore &&
                score < gridData[i][j-1].vertexScore &&
                score < gridData[i][j+1].vertexScore &&
                score < gridData[i-1][j-1].vertexScore &&
                score < gridData[i+1][j+1].vertexScore &&
                score < gridData[i+1][j-1].vertexScore &&
                score < gridData[i-1][j+1].vertexScore) {
                VertexScoreCalculator<2> calculator(slice(vertices, gridData[i][j].neighboringCells), 1.0 / grid.stepSize());
                auto F = [&calculator](Vertex<2> x) {
                    return calculator.calculateScoreSum(x);
                };
                auto vertex = steepestGradient(grid.at({i,j}), F, grid.stepSize());
                auto vertex2 = optimizationBFGS(grid.at({i,j}), F, grid.stepSize());
                DBGVAR(grid.at({i,j}), score, vertex, calculator.calculateScore(vertex), vertex2);
            }
        }
    }

}


void scanArea(const VirtualGrid<2>& grid, const std::vector<Vertex<2>>& vertices) {
    VertexScoreCalculator<2> calculator(vertices, 1.0 / grid.stepSize(), 0.5);
    for (size_t i = 0; i < grid.nElements +1; ++i){
        for (size_t j = 0; j < grid.nElements +1; ++j){
            auto score = calculator.calculateScore(grid.at({i,j}));
            auto affectedVertices = slice(vertices, score.neighboringCells);
            auto refinedGrid = VirtualGrid<2>::witCenterOrigin(grid.at({i,j}), 11, 1.1*0.5 * grid.stepSize());
            preciseScan(refinedGrid, affectedVertices);
        }
    }
}

int main(int argc, char** argv) {


    std::vector<Vertex<2>> vertices = {
        {0.5, 0.55},
        {0.75, 0.5},
        {0.5, 0.75},
        {0.25, 0.5},
        {0.5, 0.25},

        {0.9, 0.1},
        {0.9, 0.9},
        {0.1, 0.1},
        {0.1, 0.9},

//        {0.2, 0.9},
    };

//    VirtualGrid<2> grid = VirtualGrid<2>::witCenterOrigin({0.5,0.5}, 10, 0.5);
//    scanArea(grid, vertices);

//    auto F = [](Vertex<2> x) {
//        return std::pow(x.normEuclid(),2);
//    };
    VertexScoreCalculator<2> calculator(slice(vertices, {1,2, 6}), 100);
    auto F = [&](Vertex<2> x) {
        return calculator.calculateScoreSum(x);
    };
    Vertex<2> vert = { 0.73, 0.73 };

    auto x1 = optimizationBFGS(vert, F);
    auto x2 = steepestGradient(vert, F);
    DBGVAR(x1, F(x1), x2, F(x2));

//    Vector<2> vec = {2, 3};
//    DBGVAR(
//       multiply(tenzorMultiply(vec, vec), vec)
//           );

//    auto vertSlice = slice(vertices, {1,2, 6});
//    VertexScoreCalculator<2> calc(vertSlice);
//    Vertex<2> vert = { 0.73, 0.73 };
//    Vertex<2> vert2 = { 0.734, 0.734 };
//    auto x1 = steepestGradient(vert, calc, 0.5);
//    auto x2 = optimizationDFP({ 0.70, 0.70 }, calc, 0.4);
//    DBGVAR(x1, x2);

//    auto vertAvg = (vertSlice[0] + vertSlice[1] + vertSlice[2]) / 3.0;
//    DBGVAR(calc.calculateScoreSum(vertAvg), calc.calculateScore(vertAvg).vertexScore);

//    MeshGenerationProblem<2> mgp;
//    std::string defaultMeshPath = "../Meshes/mesh2D.vtk";    // 2D version
//    std::string meshPath = argc <= 1 ? defaultMeshPath : argv[1];
//    auto compData = mgp.loadMesh(meshPath);
//    std::string outPath = argc <= 2 ? "../out" : argv[2];
//    std::vector<Vertex<2>> vertices = {
//        {0.5, 0.5},
//        {0.75, 0.5},
//        {0.5, 0.75},
//        {0.25, 0.5},
//        {0.5, 0.25},

//        {0.9, 0.1},
//        {0.9, 0.9},
//        {0.1, 0.1},
//        {0.1, 0.9},

//        {0.2, 0.9},
//    };
//    mgp.calculateFx(compData, vertices);
//    DBGMSG("exporting mesh");
//    mgp.exportMeshAndData(compData, outPath + "/MeshGenerator.vtk");
}
