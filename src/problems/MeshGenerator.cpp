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
    double minf;
    double F;
    double tildeF;
    size_t cellIndex;
};

MAKE_NAMED_ATTRIBUTE_TRAIT(ComputationData, "f_{min}", minf, "F^{~}=", tildeF, "F=f(x)/f_min(x)", F, "cellIndex", cellIndex);


template<size_t Dimension, typename IndexType = size_t, typename RealType = double>
class VirtualGrid {

    using VertexType = Vertex<Dimension, RealType>;
    const VertexType origin;
    const RealType diameter;
public:

    const IndexType nElements;
    VirtualGrid(const VertexType& origin, IndexType nElements, RealType diameter): origin(origin), nElements(nElements), diameter(diameter) {}

    static VirtualGrid withCenterOrigin(const VertexType& origin, IndexType nElements, RealType diameter) {
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

    const RealType treshold = 0.1;

    VertexScoreCalculator(
            const std::vector<VertexType> vertices,
            RealType treshold = 0.1
            ): vertices(vertices), treshold(treshold){
        if (vertices.empty()) {
            throw std::range_error("Vertex calculator cannot be initialized with empty vertex set.");
        }
    }

    double scoreVertex(const VertexType& xCell, const VertexType& x) const {
        return (xCell - x).normEuclid();
    }

    struct VertexScoreData {
        RealType vertexScore;
        std::vector<IndexType> neighboringCells;
    };

    VertexScoreData calculateScore(const VertexType& x) const {
        VertexScoreData result = {0.0, std::vector<size_t>()};
        result.neighboringCells.reserve(5);
        std::vector<double> f_xs(vertices.size());
        double minF = scoreVertex(x, vertices[0]);
        f_xs[0] = minF;
        size_t minI = 0;
        for (size_t i = 1; i < vertices.size(); ++i) {
            auto f_x = scoreVertex(x, vertices[i]);
            f_xs[i] = f_x;

            if(f_x < minF) {
                minF = f_x;
                minI = i;
            }
        }
        for (size_t i = 0; i < vertices.size(); ++i){
            if (f_xs[i] - minF <= treshold) {
                auto fraction = i == minI ? 1.0 : (minF / f_xs[i]);
                result.vertexScore -= fraction;
                result.neighboringCells.emplace_back(i);
            }
        }
        return result;
    }

    RealType calculateScoreSum(const VertexType& x) const {
        RealType result = 0.0;
        std::vector<double> f_xs(vertices.size());
        auto f_x0 = scoreVertex(x, vertices[0]);
        auto maxF = f_x0;
        auto minF = f_x0;

        f_xs.emplace_back();
        for (size_t i = 1; i < vertices.size(); ++i) {
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
        return std::pow(maxF - minF, 2);
    }
};

MAKE_ATTRIBUTE_TRAIT(VertexScoreCalculator<2>::VertexScoreData,
                             vertexScore,neighboringCells );


template<typename T>
class VectorView {
    std::vector<T>& data;
    std::vector<size_t> indices;
public:

    VectorView(std::vector<T>& data, std::vector<size_t> indices): data(data), indices(indices) {}

    template<typename T_>
    static VectorView<T_> wrap(std::vector<T_>& data) {
        std::vector<size_t> indices(data.size());
        for (size_t i = 0; i < indices.size(); ++i) {
            indices[i] = i;
        }
        return VectorView<T_>(data, indices);
    }

    VectorView slice(const std::vector<size_t>& indices) const {
        std::vector<size_t> newIndices;
        for (auto index : indices) {
            newIndices.emplace_back(this->indices[index]);
        }
        return VectorView(data, newIndices);
    }

    size_t size() const {
        return indices.size();
    }

    T& operator[](size_t index) {
        return data[indices[index]];
    }

    const T& operator[](size_t index) const {
        return data[indices[index]];
    }

    size_t getGlobalIndexAt(size_t index) const {
        return indices.at(index);
    }

    const std::vector<size_t>& getIndices() const {
        return indices;
    }

    std::vector<T> toVector() const {
        std::vector<T> result;
        result.reserve(indices.size());
        for (auto index : indices){
            result.emplace_back(data[index]);
        }
        return result;
    }
};


template<unsigned int ProblemDimension>
struct MeshGenerationProblem{
    using MeshType = UnstructuredMesh<ProblemDimension, size_t, double>;
    using ResultType = ComputationData;
    using ProblemDataContainerType = MeshDataContainer<ResultType, ProblemDimension>;
    std::shared_ptr<MeshReader<ProblemDimension>> meshReader;
    VTKMeshWriter<ProblemDimension, size_t, double> meshWriter;

    MeshType mesh;

    void calculateFx(ProblemDataContainerType& data, std::vector<Vertex<ProblemDimension, double>> vertices) {
        std::vector<double> f_xs(vertices.size());
        double treshold = 0.025;
        VertexScoreCalculator<ProblemDimension> calculator(vertices, treshold);
        for (const auto& cell : mesh.getCells()) {
            auto score = calculator.calculateScore(cell.getCenter());
            data[cell].F = score.vertexScore;
            auto vectorSlice = VectorView<Vertex<ProblemDimension>>(vertices, score.neighboringCells).toVector();
            data[cell].tildeF = VertexScoreCalculator<ProblemDimension>(vectorSlice, treshold).calculateScoreSum(cell.getCenter());
            double minF = calculator.scoreVertex(cell.getCenter(), vertices[0]);
            size_t minI = 0;
            for (size_t i = 1; i < vertices.size(); ++i) {
                auto vertF = calculator.scoreVertex(cell.getCenter(), vertices[i]);
                if (vertF < minF) {
                    minF = vertF;
                    minI = i;
                }
            }
            data[cell].minf = minF;
            data[cell].cellIndex = minI;
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
    double delta = maxDist * 1e-3;
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

Vector<1> multiply(const Vector<1,Vector<1>>& A, const Vector<1>& x){
    return Vector<1>{A[0][0] * x[0]};
}

double multiply(const Vector<1>& x, const Vector<1>& y){
    return x[0] * y[0];
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

Vector<1,Vector<1>> tenzorMultiply(const Vector<1>& x, const Vector<1>& y_t) {
    return {{x[0] * y_t[0]}};
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

Vector<1,Vector<1>> operator*(double alpha, const Vector<1,Vector<1>>& A) {
    return Vector<1, Vector<1>>{{A[0][0] * alpha}};
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
        S -= tenzorMultiply(S_q, S_q) / multiply(q, S_q) - tenzorMultiply(p,p) / multiply(p, p);
    }
}

template<typename F>
Vertex<2> optimizationBFGS(Vertex<2> x, F&& f, double maxDist = 0.01, const double treshold = 1e-5) {
    double delta = maxDist * 1e-3;
    double invDelta = 1.0 / (2*delta);

    auto x_i = x;

    double fx = f(x);
    double df_dx = (f(x_i-Vertex<2>({-delta, 0})) - f(x_i-Vertex<2>({delta, 0}))) * invDelta;
    double df_dy = (f(x_i-Vertex<2>({0, -delta})) - f(x_i-Vertex<2>({0, delta}))) * invDelta;
    Vector<2> g_i = {df_dx, df_dy};
    Vector<2,Vector<2>> S {{1,0},{0,1}};
    int cnt = 0;
    while(true){
        if(++cnt > 100 || g_i.normEuclid() < treshold) {
            DBGVAR(g_i.normEuclid(), cnt);
            return x_i;
        }
        auto d = -1.0 * multiply(S,  g_i);

        double gamma = maxDist / d.normEuclid();
        double auxFx = 0;
        Vertex<2> x_i1;
        // line search
        do {
            x_i1 = x_i + gamma * d;
            auxFx = f(x_i1);
            gamma *= 0.7;
        } while (auxFx > fx);
        fx = auxFx;
        double df_dx = (f(x_i1-Vertex<2>({-delta, 0})) - f(x_i1-Vertex<2>({delta, 0}))) * invDelta;
        double df_dy = (f(x_i1-Vertex<2>({0, -delta})) - f(x_i1-Vertex<2>({0, delta}))) * invDelta;
        auto g_i1 = Vector<2>{df_dx, df_dy};
        auto p = x_i1 - x_i;
        auto q = g_i1 - g_i;
        x_i = x_i1;
        g_i = g_i1;
        auto p_q = multiply(p,q);
        if (abs(p_q) < 1e-30) {
            DBGVAR(p_q, cnt);
            return x_i;
        }
        auto S_q = multiply(S,q);
        auto deltaS = (tenzorMultiply(p, S_q) + tenzorMultiply(S_q, p) - (1 + (multiply(q, S_q)/p_q)) * tenzorMultiply(p,p)) / p_q;
        S -= deltaS;
    }
}


namespace std {
    template<unsigned int Dimension, typename T>
    class hash<Vertex<Dimension, T>> {
        public:
        size_t operator()(const Vertex<Dimension, T>& vert) const noexcept {
            const std::hash<T> hashT{};
            size_t hash = 0;
            for (size_t i = 0; i < Dimension; ++i) {
                hash ^= hashT(vert[i]) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            }
            return hash;
        }
    };

    template<typename T>
    class hash<std::vector<T>> {
        public:
        size_t operator()(const std::vector<T>& vec) const noexcept {
            const std::hash<T> hashT{};
            size_t hash = 0;
            for (size_t i = 0; i < vec.size(); ++i) {
                hash ^= hashT(vec[i]) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            }
            return hash;
        }
    };
}

std::unordered_map<std::vector<size_t>, Vertex<2,double>> cache;

void locateVertex(const Vertex<2>& startVertex, const VectorView<Vertex<2>>& vertices) {

    auto cachedValueIter = cache.find(vertices.getIndices());
    if (cachedValueIter != cache.end()) {
        DBGVAR(*cachedValueIter, vertices);
        return;
    }

    auto meanVertex = vertices[0];
    for (size_t i = 1; i < vertices.size(); ++i) {
        meanVertex += vertices[i];
    }
    meanVertex /= vertices.size();
    double meanDist = 0.0;
    double maxDist = 0.0;
    for (size_t i = 0; i < vertices.size(); ++i) {
        double dist = (meanVertex - vertices[i]).normEuclid();
        meanDist += dist;
        if (maxDist < dist) {
            maxDist = dist;
        }
    }
    meanDist /= vertices.size();

    if ((startVertex - meanVertex).normEuclid() > maxDist){
        return;
    }

    VertexScoreCalculator<2> calculator(vertices.toVector(), 1.0);
    auto F = [&calculator](Vertex<2> x) {
        return calculator.calculateScoreSum(x);
    };
    auto vertex2 = optimizationBFGS(startVertex, F, meanDist, 1e-8);
    cache[vertices.getIndices()] = vertex2;
    DBGVAR(startVertex, calculator.calculateScore(vertex2), vertex2, vertices.getIndices(), vertices);
}


void preciseScan(const VirtualGrid<2>& grid, const VectorView<Vertex<2>>& vertices) {
    std::vector<std::vector<VertexScoreCalculator<2>::VertexScoreData>> gridData;
    double treshold = (1.0 - grid.stepSize()) / (1.0 + grid.stepSize());
    VertexScoreCalculator<2> calculator(vertices.toVector(), treshold);
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
            auto [score, n] = gridData[i][j];
            if (n.size() >= 3) {
                auto localVertices = vertices.slice(gridData[i][j].neighboringCells);
                locateVertex(grid.at({i,j}), localVertices);
            }
        }
    }

}

template<typename T>
std::vector<T> cat(std::vector<T> v1, const std::vector<T>& v2) {
    v1.insert(v1.end(), v2.begin(), v2.end());
    return v1;
}

void scanArea(const VirtualGrid<2>& grid, const VectorView<Vertex<2>>& vertices, const double maxPrecision = 1e-2) {
    double treshold = 2.0 * (1.414213562373095) * grid.stepSize();
    VertexScoreCalculator<2> calculator(vertices.toVector(), treshold);
    for (size_t i = 0; i < grid.nElements +1; ++i){
        for (size_t j = 0; j < grid.nElements +1; ++j){
            auto score = calculator.calculateScore(grid.at({i,j}));
            if (score.neighboringCells.size() > 2){
                auto affectedVertices = vertices.slice(score.neighboringCells);
                auto refinedGrid = VirtualGrid<2>::withCenterOrigin(grid.at({i,j}), 11, 1.1*0.5 * grid.stepSize());
                if (affectedVertices.size() == 3) {
                    locateVertex(grid.at({i,j}), affectedVertices);
                    continue;
                }
                if (refinedGrid.stepSize() < maxPrecision){
                    preciseScan(refinedGrid, affectedVertices);
                } else {
                    scanArea(refinedGrid, vertices, maxPrecision);
                }
            }
        }
    }
}


template<typename F>
Vertex<1> optimizationBFGSBoundary(Vertex<1> x, F&& f, const std::function<Vertex<2>(const Vertex<1>&)> boundaryMapping, double maxDist = 0.01, const double treshold = 1e-5) {
    double delta = maxDist * 1e-3;
    double invDelta = 1.0 / (2*delta);

    auto x_i = x;

    double fx = f(boundaryMapping(x));
    double df_dx = (f(boundaryMapping(x_i-Vertex<1>({-delta}))) - f(boundaryMapping(x_i-Vertex<1>({delta})))) * invDelta;
    Vector<1> g_i = {df_dx};
    Vector<1,Vector<1>> S {{1}};
    int cnt = 0;
    while(true){
        if(++cnt > 100 || g_i.normEuclid() < treshold) {
            return x_i;
        }
        auto d = -1.0 * multiply(S,  g_i);

        double gamma = maxDist / d.normEuclid();
        double auxFx = 0;
        Vertex<1> x_i1;
        // line search
        do {
            x_i1 = x_i + gamma * d;
            auxFx = f(boundaryMapping(x_i1));
            gamma *= 0.7;
        } while (auxFx > fx);
        fx = auxFx;
        double df_dx = (f(boundaryMapping(x_i1-Vertex<1>({-delta}))) - f(boundaryMapping(x_i1-Vertex<1>({delta})))) * invDelta;
        auto g_i1 = Vector<1>{df_dx};
        auto p = x_i1 - x_i;
        auto q = g_i1 - g_i;
        x_i = x_i1;
        g_i = g_i1;
        auto p_q = multiply(p,q);
        if (abs(p_q) < 1e-30) {
            return x_i;
        }
        auto S_q = multiply(S,q);
        auto deltaS = (tenzorMultiply(p, S_q) + tenzorMultiply(S_q, p) - (1 + (multiply(q, S_q)/p_q)) * tenzorMultiply(p,p)) / p_q;
        S -= deltaS;
    }
}

void locateVertexBoundary(const Vertex<1>& startVertex, const std::function<Vertex<2>(const Vertex<1>&)> boundaryMapping, const VectorView<Vertex<2>>& vertices, const std::vector<size_t>& boundaryIndices) {
    auto cachedValueIter = cache.find(vertices.getIndices());
    if (cachedValueIter != cache.end()) {
        // DBGVAR(*cachedValueIter, vertices);
        return;
    }

    auto meanVertex = vertices[0];
    for (size_t i = 1; i < vertices.size(); ++i) {
        meanVertex += vertices[i];
    }
    meanVertex /= vertices.size();
    double meanDist = 0.0;
    double maxDist = 0.0;
    for (size_t i = 0; i < vertices.size(); ++i) {
        double dist = (meanVertex - vertices[i]).normEuclid();
        meanDist += dist;
        if (maxDist < dist) {
            maxDist = dist;
        }
    }
    meanDist /= vertices.size();

    if ((boundaryMapping(startVertex) - meanVertex).normEuclid() > maxDist){
        return;
    }

    VertexScoreCalculator<2> calculator(vertices.toVector(), 1.0);
    auto F = [&calculator](Vertex<2> x) {
        return calculator.calculateScoreSum(x);
    };
    auto vertex2 = optimizationBFGSBoundary(startVertex, F, boundaryMapping, meanDist, meanDist * 1e-3);
    cache[std::move(cat(vertices.getIndices(), boundaryIndices))] = boundaryMapping(vertex2);
}

void scanAreaBoundary(const VirtualGrid<1>& grid, const std::function<Vertex<2>(const Vertex<1>&)> boundaryMapping, const VectorView<Vertex<2>>& vertices, const std::vector<size_t>& boundaryIndices, const double maxPrecision = 1e-2) {
    double treshold = 2.0 * grid.stepSize();
    VertexScoreCalculator<2> calculator(vertices.toVector(), treshold);
    for (size_t i = 0; i < grid.nElements +1; ++i){
        auto score = calculator.calculateScore(boundaryMapping(grid.at({i})));
        if (score.neighboringCells.size() > 1){
            auto affectedVertices = vertices.slice(score.neighboringCells);
            auto refinedGrid = VirtualGrid<1>::withCenterOrigin(grid.at({i}), 11, 1.1*0.5 * grid.stepSize());
            if (affectedVertices.size() == 2) {
                locateVertexBoundary(grid.at({i}), boundaryMapping, affectedVertices, boundaryIndices);
                continue;
            }
            if (refinedGrid.stepSize() < maxPrecision){
                locateVertexBoundary(grid.at({i}), boundaryMapping, affectedVertices, boundaryIndices);
            } else {
                scanAreaBoundary(refinedGrid, boundaryMapping, vertices, boundaryIndices, maxPrecision);
            }
        }
    }
}


auto buildMeshFromCache(std::unordered_map<std::vector<size_t>, Vertex<2,double>>& cache) {
    std::unordered_map<size_t, std::unordered_map<size_t,std::vector<std::size_t>>> buildCache;
    std::vector<std::reference_wrapper<Vertex<2>>> vertices;
    vertices.reserve(buildCache.size());

    for(auto item : cache) {
        vertices.emplace_back(cache[item.first]);
        for(const auto index : item.first) {
            auto cachedData = buildCache.find(index);
            if (cachedData == buildCache.end()) {
                cachedData = buildCache.emplace(index, std::unordered_map<size_t,std::vector<std::size_t>>()).first;
            }
            auto& cachedNeighbors = cachedData->second;
            for (const auto index2 : item.first){
                if (index2 == index) {
                    continue;
                }
                auto cachedNeighborIter = cachedNeighbors.find(index2);
                if (cachedNeighborIter == cachedNeighbors.end()) {
                    cachedNeighborIter = cachedNeighbors.emplace(index2, std::vector<std::size_t>()).first;
                }
                cachedNeighborIter->second.emplace_back(vertices.size() -1);
            }
        }
    }
    DBGCHECK;
    DBGVAR(buildCache);
    using MeshType = UnstructuredMesh<2, size_t, double>;
    MeshType mesh;
    DBGVAR(vertices);
    for (auto ref : vertices){
        DBGVAR(ref.get());
        mesh.getVertices().emplace_back(MeshType::Vertex(mesh.getVertices().size(), ref.get()));
    }
    std::unordered_map<std::vector<size_t>, std::reference_wrapper<MeshType::Edge>> edgeCache;
    DBGCHECK;
    for (const auto& cellCache : buildCache) {
        if (isBoundaryIndex(cellCache.first)) {
            continue;
        }
        auto cellIndex = mesh.getCells().size();
        mesh.getCells().push_back(cellIndex);
        auto & cell = mesh.getCells().back();
        size_t lastEdgeIndex = INVALID_INDEX(size_t);
        size_t firstEdgeIndex;
        for(const auto& cellBoundary : cellCache.second) {
            if(cellBoundary.second.size() >= 2 ){
                auto cachedEdge = edgeCache.find(cellBoundary.second);
                size_t newIndex;
                if (cachedEdge != edgeCache.end()) {
                    newIndex = cachedEdge->second.get().getIndex();
                } else {
                    newIndex = mesh.getEdges().size();
                    mesh.getEdges().emplace_back(MeshType::Edge(newIndex, cellBoundary.second[0], cellBoundary.second[1]));
                }
                auto & edge = mesh.getEdges()[newIndex];
                if (lastEdgeIndex == INVALID_INDEX(size_t)){
                    cell.setBoundaryElementIndex(newIndex);
                    firstEdgeIndex = newIndex;
                } else {
                    edge.setNextBElem(lastEdgeIndex, cellIndex);
                }
                lastEdgeIndex = newIndex;
            }
        }
        mesh.getEdges()[firstEdgeIndex].setNextBElem(lastEdgeIndex, cellIndex);
    }
    DBGCHECK;
    DBGVAR(mesh.getVertices());
    return mesh;
}

int main(int argc, char** argv) {

    std::vector<Vertex<2>> vertices = {
        {0.5, 0.5},
        {0.75, 0.5},
        {0.5, 0.75},
        {0.25, 0.5},
        {0.5, 0.25},

        {0.9, 0.1},
        {0.9, 0.9},
        {0.1, 0.1},
        {0.1, 0.9},

    //    {0.2, 0.9},
    };
    VirtualGrid<2> grid = VirtualGrid<2>::withCenterOrigin({0.5,0.5}, 10, 0.5);
    auto verticesView = VectorView<Vertex<2>>::wrap(vertices);
    scanArea(grid, verticesView);
    auto gridBoundary = VirtualGrid<1>::withCenterOrigin({0.5}, 10, 0.5);
    auto b1 = [](const Vertex<1>& v){
        return Vertex<2>{v[0], 0.0};
    };
    auto b2 = [](const Vertex<1>& v){
        return Vertex<2>{v[0], 1.0};
    };

    auto b3 = [](const Vertex<1>& v){
        return Vertex<2>{0.0, v[0]};
    };

    auto b4 = [](const Vertex<1>& v){
        return Vertex<2>{1.0, v[0]};
    };
    // makeBoundaryIndex(0)
    scanAreaBoundary(gridBoundary, b1, verticesView, {makeBoundaryIndex(0ul)});
    scanAreaBoundary(gridBoundary, b2, verticesView, {makeBoundaryIndex(1ul)});
    scanAreaBoundary(gridBoundary, b3, verticesView, {makeBoundaryIndex(2ul)});
    scanAreaBoundary(gridBoundary, b4, verticesView, {makeBoundaryIndex(3ul)});
    VertexScoreCalculator<2> calculator(vertices, 1e-4);
    cache[cat(calculator.calculateScore(Vertex<2>{0.0,0.0}).neighboringCells, {makeBoundaryIndex(0ul), makeBoundaryIndex(2ul)})] = Vertex<2>{0.0,0.0};//0 2
    cache[cat(calculator.calculateScore(Vertex<2>{1.0,0.0}).neighboringCells, {makeBoundaryIndex(0ul), makeBoundaryIndex(3ul)})] = Vertex<2>{1.0,0.0};//0 3
    cache[cat(calculator.calculateScore(Vertex<2>{0.0,1.0}).neighboringCells, {makeBoundaryIndex(1ul), makeBoundaryIndex(2ul)})] = Vertex<2>{0.0,1.0};//1 2
    cache[cat(calculator.calculateScore(Vertex<2>{1.0,1.0}).neighboringCells, {makeBoundaryIndex(1ul), makeBoundaryIndex(3ul)})] = Vertex<2>{1.0,1.0};//1 3
    DBGVAR(cache);

    auto mesh = buildMeshFromCache(cache);
    mesh.write("test-mesh.vtk");
//    auto F = [](Vertex<2> x) {
//        return std::pow(x.normEuclid(),2);
//    };
//    VertexScoreCalculator<2> calculator(slice(vertices, {1,2, 6}), 100);
//    auto F = [&](Vertex<2> x) {
//        return calculator.calculateScoreSum(x);
//    };
//    Vertex<2> vert = { 0.73, 0.73 };

//    auto x1 = optimizationBFGS(vert, F);
//    auto x2 = steepestGradient(vert, F);
//    DBGVAR(x1, F(x1), x2, F(x2));

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

//    mgp.calculateFx(compData, vertices);
//    DBGMSG("exporting mesh");
//    mgp.exportMeshAndData(compData, outPath + "/MeshGenerator.vtk");
}
