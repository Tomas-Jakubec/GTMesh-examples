#include <memory>
#include <limits>
#include <iomanip>
#include <assert.h>
#include <list>
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

template<typename T>
class VectorView {
    std::vector<T>& data;
    std::vector<size_t> indices;
public:

    VectorView(std::vector<T>& data, std::vector<size_t> indices): data(data), indices(indices) {}

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

template<typename T>
static VectorView<T> wrapToVectorView(std::vector<T>& data) {
    std::vector<size_t> indices(data.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        indices[i] = i;
    }
    return VectorView<T>(data, indices);
}


template<unsigned int Dimension, typename IndexType = size_t, typename RealType = double>
class VirtualGrid {

    using VertexType = Vertex<Dimension, RealType>;
    const VertexType origin;
public:

    const IndexType nElements;
    const RealType stepSize;
    VirtualGrid(const VertexType& origin, IndexType nElements, RealType stepSize): origin(origin), nElements(nElements), stepSize(stepSize) {}

    static VirtualGrid withCenterOrigin(const VertexType& origin, IndexType nElements, RealType diameter) {
        VertexType shift;
        for(auto & coordinate: shift) {
            coordinate = diameter;
        }
        VertexType originShifted = origin - shift;
        auto stepSize = 2 * diameter / nElements;
        return VirtualGrid(originShifted, nElements + 1, stepSize);
    }

    VertexType at(const std::array<IndexType, Dimension>& indices) const {
        VertexType result = origin;
        for(unsigned int i = 0; i < Dimension; ++i) {
            assert(indices[i] < nElements);
            result[i] += stepSize * indices[i];
        }
        return result;
    }
};

template<size_t Dimension, typename IndexType = size_t, typename RealType = double>
class VertexScoreCalculator {

public:
    using VertexType = Vertex<Dimension, RealType>;
    const VectorView<VertexType>& vertices;

    const RealType treshold = 0.1;

    VertexScoreCalculator(
            const VectorView<VertexType>& vertices,
            RealType treshold = 0.1
            ): vertices(vertices), treshold(treshold){
    }

    VertexScoreCalculator withVertices(const VectorView<VertexType>& vertices) {
        return VertexScoreCalculator(vertices, treshold);
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

    auto calculateNeighbors(const VertexType& x, const double treshold) const {
        auto result = std::vector<size_t>();
        result.reserve(5);
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
                result.emplace_back(i);
            }
        }
        return result;
    }

    RealType calculateScoreSum(const VertexType& x) const {
        RealType result = 0.0;
        auto f_x0 = scoreVertex(x, vertices[0]);
        auto maxF = f_x0;
        auto minF = f_x0;

        for (size_t i = 1; i < vertices.size(); ++i) {
            const auto& vert = vertices[i];

            auto f_x = scoreVertex(x, vert);
            if(f_x > maxF) {
                maxF = f_x;
            }
            if(f_x < minF) {
                minF = f_x;
            }
        }
        return pow(maxF - minF, 2) * maxF;
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

    void calculateFx(ProblemDataContainerType& data, std::vector<Vertex<ProblemDimension, double>> vertices) {
        std::vector<double> f_xs(vertices.size());
        double treshold = 0.025;
        auto verticesView = wrapToVectorView(vertices);
        VertexScoreCalculator<ProblemDimension> calculator(verticesView, treshold);
        for (const auto& cell : mesh.getCells()) {
            auto score = calculator.calculateScore(cell.getCenter());
            data[cell].F = score.vertexScore;
            auto vectorSlice = verticesView.slice(score.neighboringCells);
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

template <unsigned int Dim, unsigned int Index>
constexpr double basisComponent() {
    return (Index == Dim) ? 1.0 : 0.0;
}

template <unsigned int Dim>
constexpr void fillBasis(Vertex<Dim>& basis, std::integral_constant<unsigned int, 0>) {
    basis[0] = basisComponent<Dim, 0>();
}

template <unsigned int Dim, unsigned int Index>
constexpr void fillBasis(Vertex<Dim>& basis, std::integral_constant<unsigned int, Index>) {
    fillBasis<Dim>(basis, std::integral_constant<unsigned int, Index - 1>());
    basis[Index] = basisComponent<Dim, Index>();
}

template <unsigned int Dim, unsigned int index>
constexpr Vertex<Dim> e() {
    Vertex<Dim> basis {};
    fillBasis<Dim>(basis, std::integral_constant<unsigned int, Dim - 1>());
    return basis;
}


template <unsigned int Index, unsigned int Dim, typename F>
double calculateGradientComponent(F&& f, const Vertex<Dim>& x, double delta) {
    Vertex<Dim> deltaVec = e<Dim, Index>() * delta;
    return (f(x + deltaVec) - f(x - deltaVec)) / (2 * delta);
}


template <typename F, unsigned int Dim>
constexpr void calculateGradient(Vector<Dim>& result, F&& f, const Vertex<Dim>& x, double delta, std::integral_constant<unsigned int, 0>) {
    result[0] = calculateGradientComponent<0>(f, x, delta);
}
template <typename F, unsigned int Dim, unsigned int Index = 0>
constexpr void calculateGradient(Vector<Dim>& result, F&& f, const Vertex<Dim>& x, double delta, std::integral_constant<unsigned int, Index>) {
    result[Index] = calculateGradientComponent<Index>(f, x, delta);
    calculateGradient(result, f, x, delta, std::integral_constant<unsigned int, Index - 1>());
}


template <unsigned int Dim, typename F>
Vector<Dim> calculateGradient(F&& f, const Vertex<Dim>& x, double delta = 1e-5) {
    Vector<Dim> gradient;

    calculateGradient(gradient, f, x, delta, std::integral_constant<unsigned int, Dim - 1>());

    return gradient;
}


template <unsigned int Dim>
void unitMatrix(Vector<Dim,Vector<Dim>>& matrix, std::integral_constant<unsigned int, 0>) {
    matrix[0] = e<Dim, 0>();
}

template <unsigned int Dim, unsigned int Index>
void unitMatrix(Vector<Dim,Vector<Dim>>& matrix, std::integral_constant<unsigned int, Index>) {
    unitMatrix(matrix, std::integral_constant<unsigned int, Index -1>());
    matrix[Index] = e<Dim, Index>();
}


template <unsigned int Dim>
Vector<Dim,Vector<Dim>> unitMatrix() {
    Vector<Dim,Vector<Dim>> matrix;
    unitMatrix(matrix, std::integral_constant<unsigned int, Dim - 1>());
    return matrix;
}

template<unsigned int Dim, typename F>
Vertex<Dim> optimizationBFGS(Vertex<Dim> x, F&& f, double maxDist = 0.01, const double treshold = 1e-5) {
    double delta = maxDist * 1e-3;
    double invDelta = 1.0 / (2*delta);

    auto x_i = x;

    double fx = f(x);
    Vector<Dim> g_i = calculateGradient(f, x_i, delta);
    auto S = unitMatrix<Dim>();
    int cnt = 0;
    while(true){
        if(++cnt > 100 || g_i.normEuclid() < treshold) {
            DBGVAR(g_i.normEuclid(), cnt);
            return x_i;
        }
        auto d = -1.0 * multiply(S,  g_i);

        double gamma = maxDist / d.normEuclid();
        double auxFx = 0;
        
        Vertex<Dim> x_i1;
        // line search
        do {
            x_i1 = x_i + gamma * d;
            auxFx = f(x_i1);
            gamma *= 0.7;
        } while (auxFx > fx);
        fx = auxFx;
        auto g_i1 = calculateGradient(f, x_i1, delta);
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


class Cache {
    using cacheElementT = std::pair<std::vector<std::size_t>, Vertex<2>>;
    std::vector<cacheElementT> cacheData;
    std::unordered_map<std::size_t, std::set<std::size_t>> cacheMap;
public:
    void insert(const std::vector<std::size_t>& indices, Vertex<2> vertex) {
        auto newIndex = cacheData.size();
        cacheData.emplace_back(std::make_pair(indices, vertex));
        for (const auto index : indices) {
            auto indexIter = cacheMap.find(index);
            if (indexIter == cacheMap.end()) {
                indexIter = cacheMap.emplace(index, std::set<std::size_t>()).first;
            }
            indexIter->second.emplace(newIndex);
        }
    }

    void insert(std::vector<std::size_t>&& indices, Vertex<2> vertex) {
        auto newIndex = cacheData.size();
        cacheData.emplace_back(std::make_pair(indices, vertex));
        for (const auto index : indices) {
            auto indexIter = cacheMap.find(index);
            if (indexIter == cacheMap.end()) {
                indexIter = cacheMap.emplace(index, std::set<std::size_t>()).first;
            }
            indexIter->second.emplace(newIndex);
        }
    }
    
    auto find(const std::vector<std::size_t>& indices) {
        auto cachedValueIter = cacheMap.find(indices[0]);
        std::set<std::size_t> result;
        if (cachedValueIter != cacheMap.end()) {
            result = cachedValueIter->second;
        } else {
            return cacheData.end();
        }
        for (std::size_t i = 1; i < indices.size(); ++i) {
            auto indexValIter = cacheMap.find(indices[i]);
            if (indexValIter != cacheMap.end()) {
                auto& indexSet = indexValIter->second;
                std::set<std::size_t> resultIntersection;
                std::set_intersection(result.begin(), result.end(), indexSet.begin(), indexSet.end(), std::inserter(resultIntersection, resultIntersection.begin())); 
                result = resultIntersection;
            } else {
                return cacheData.end();
            }   
            
        }
        if (result.size() == 0) {
            return cacheData.end();
        }
        if (result.size() != 1) {
            DBGVAR(indices, cacheData, result);
            throw std::runtime_error("Cache invalid state -> multiple results for indices");
        }
        return cacheData.begin() + *result.begin();
    }


    auto begin() {
        return cacheData.begin();
    }

    auto end() {
        return cacheData.end();
    }

    auto begin() const {
        return cacheData.cbegin();
    }

    auto end() const {
        return cacheData.cend();
    }
    
    auto at(const std::vector<std::size_t>& indices) {
        auto iter = find(indices);
        if (iter == end()) {
            throw std::out_of_range("Indices not found in the cache!");
        }
        return *iter;
    }
};

Cache cache;

template <typename T>
bool isSubset(const std::vector<T>& set1, const std::vector<T>& set2) {
    std::vector<T> intersection;
    std::set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(), std::back_inserter(intersection));
    return intersection.size() == set1.size();
}

template <unsigned int Dim>
void locateVertex(const Vertex<Dim>& startVertex, const VectorView<Vertex<Dim>>& vertices, const VectorView<Vertex<Dim>>& enclosingVertices) {

    auto cachedValueIter = cache.find(vertices.getIndices());
    if (cachedValueIter != cache.end()) {
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

    VertexScoreCalculator<2> calculator(vertices, 1.0);
    auto F = [&calculator](Vertex<2> x) {
        return calculator.calculateScoreSum(x);
    };
    auto vertex2 = optimizationBFGS(meanVertex, F, meanDist, 1e-5);
    // Calculate the real neghborhood
    auto neighbors = calculator.withVertices(enclosingVertices).calculateNeighbors(vertex2, 1e-5);
    auto actualVertices = enclosingVertices.slice(neighbors);
    if (isSubset(vertices.getIndices(), actualVertices.getIndices())){
        cache.insert(neighbors, vertex2);
    }
    DBGVAR(startVertex, meanVertex, calculator.calculateScore(vertex2), vertex2, vertices.getIndices(), vertices, actualVertices.getIndices());
}

template<typename T>
std::vector<T> cat(std::vector<T> v1, const std::vector<T>& v2) {
    v1.insert(v1.end(), v2.begin(), v2.end());
    return v1;
}

template <unsigned int Dim>
void scanArea(const VirtualGrid<Dim>& grid, const VectorView<Vertex<Dim>>& vertices, const double maxPrecision = 1e-2) {
    double treshold = 2.0 * (1.414213562373095) * grid.stepSize;
    VertexScoreCalculator<Dim> calculator(vertices, treshold);
    std::array<size_t, Dim> indices {};
    for (size_t n = 0; n < pow(grid.nElements, Dim); ++n) {
        auto neighboringCells = calculator.calculateNeighbors(grid.at(indices), treshold);
        if (neighboringCells.size() > Dim){
            DBGVAR(grid.at(indices), neighboringCells, indices);
            auto affectedVertices = vertices.slice(neighboringCells);
            if (affectedVertices.size() == Dim + 1) {
                locateVertex(grid.at(indices), affectedVertices, vertices);
            } else {
                auto refinedGrid = VirtualGrid<2>::withCenterOrigin(grid.at(indices), 10, 0.5 * grid.stepSize);
                if (refinedGrid.stepSize < maxPrecision){
                    locateVertex(grid.at(indices), affectedVertices, vertices);
                } else {
                    scanArea(refinedGrid, vertices, maxPrecision);
                }
            }
        }

        indices[0]++;
        for(size_t i = 0; i < Dim - 1; ++i){
            if(indices[i] == grid.nElements) {
                indices[i] = 0;
                indices[i+1]++;
            }
        }
    }
}

void locateVertexBoundary(const Vertex<1>& startVertex, const std::function<Vertex<2>(const Vertex<1>&)> boundaryMapping, const VectorView<Vertex<2>>& vertices, const std::vector<size_t>& boundaryIndices) {
    auto cachedValueIter = cache.find(cat(vertices.getIndices(), boundaryIndices));
    if (cachedValueIter != cache.end()) {
        DBGVAR(*cachedValueIter);
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

    VertexScoreCalculator<2> calculator(vertices, 1.0);
    
    auto F = [calculator, boundaryMapping](const Vertex<1>& x) {
        return calculator.calculateScoreSum(boundaryMapping(x));
    };
    
    auto vertex2 = optimizationBFGS(startVertex, F, meanDist, meanDist * 1e-3);
    
    cache.insert(std::move(cat(vertices.getIndices(), boundaryIndices)), boundaryMapping(vertex2));
}

void scanAreaBoundary(const VirtualGrid<1>& grid, const std::function<Vertex<2>(const Vertex<1>&)> boundaryMapping, const VectorView<Vertex<2>>& vertices, const std::vector<size_t>& boundaryIndices, const double maxPrecision = 1e-2) {
    double treshold = 2.0 * grid.stepSize;
    VertexScoreCalculator<2> calculator(vertices, treshold);
    for (size_t i = 0; i < grid.nElements; ++i){
        auto score = calculator.calculateScore(boundaryMapping(grid.at({i})));
        if (score.neighboringCells.size() > 1){
            auto affectedVertices = vertices.slice(score.neighboringCells);
            if (affectedVertices.size() == 2) {
                locateVertexBoundary(grid.at({i}), boundaryMapping, affectedVertices, boundaryIndices);
                continue;
            }
            auto refinedGrid = VirtualGrid<1>::withCenterOrigin(grid.at({i}), 10, 0.5 * grid.stepSize);
            if (refinedGrid.stepSize < maxPrecision){
                locateVertexBoundary(grid.at({i}), boundaryMapping, affectedVertices, boundaryIndices);
            } else {
                scanAreaBoundary(refinedGrid, boundaryMapping, vertices, boundaryIndices, maxPrecision);
            }
        }
    }
}


auto buildMeshFromCache(Cache& cache) {
    std::unordered_map<size_t, std::unordered_map<size_t,std::vector<std::size_t>>> buildCache;
    std::vector<Vertex<2>> vertices;
    vertices.reserve(buildCache.size());

    for(auto item : cache) {
        DBGVAR(item);
        vertices.emplace_back(item.second);
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
    DBGVAR(buildCache);
    using MeshType = UnstructuredMesh<2, size_t, double>;
    MeshType mesh;
    for (auto vert : vertices){
        mesh.getVertices().emplace_back(MeshType::Vertex(mesh.getVertices().size(), vert));
    }
    std::unordered_map<std::vector<size_t>, std::reference_wrapper<MeshType::Edge>> edgeCache;
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
    DBGVAR(mesh.getVertices());
    return mesh;
}

void testGenerate() {
    // std::vector<Vertex<2>> vertices = {
    //     {0.5, 0.5},
    //     {0.75, 0.5},
    //     {0.5, 0.75},
    //     {0.25, 0.5},
    //     {0.5, 0.25},

    //     {0.9, 0.1},
    //     {0.9, 0.9},
    //     {0.1, 0.1},
    //     {0.1, 0.9},

    //    {0.2, 0.9},
    // };
    std::vector<Vertex<2>> vertices;
    auto wg1 = VirtualGrid<2>::withCenterOrigin({0.5009,0.5015}, 5, 0.45);

    for (size_t i = 0; i < wg1.nElements; ++i){ 
        for (size_t j = 0; j < wg1.nElements; ++j){
            vertices.emplace_back(wg1.at({i, j}));
        }
    }

    VirtualGrid<2> grid = VirtualGrid<2>::withCenterOrigin({0.5,0.5}, 100, 0.5);
    auto verticesView = wrapToVectorView(vertices);
    scanArea(grid, verticesView);
    auto gridBoundary = VirtualGrid<1>::withCenterOrigin({0.5}, 100, 0.5);
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
    scanAreaBoundary(gridBoundary, b1, verticesView, {makeBoundaryIndex(0ul)});
    scanAreaBoundary(gridBoundary, b2, verticesView, {makeBoundaryIndex(1ul)});
    scanAreaBoundary(gridBoundary, b3, verticesView, {makeBoundaryIndex(2ul)});
    scanAreaBoundary(gridBoundary, b4, verticesView, {makeBoundaryIndex(3ul)});
    VertexScoreCalculator<2> calculator(verticesView, 1e-4);
    cache.insert(cat(calculator.calculateScore(Vertex<2>{0.0,0.0}).neighboringCells, {makeBoundaryIndex(0ul), makeBoundaryIndex(2ul)}), Vertex<2>{0.0,0.0});//0 2
    cache.insert(cat(calculator.calculateScore(Vertex<2>{1.0,0.0}).neighboringCells, {makeBoundaryIndex(0ul), makeBoundaryIndex(3ul)}), Vertex<2>{1.0,0.0});//0 3
    cache.insert(cat(calculator.calculateScore(Vertex<2>{0.0,1.0}).neighboringCells, {makeBoundaryIndex(1ul), makeBoundaryIndex(2ul)}), Vertex<2>{0.0,1.0});//1 2
    cache.insert(cat(calculator.calculateScore(Vertex<2>{1.0,1.0}).neighboringCells, {makeBoundaryIndex(1ul), makeBoundaryIndex(3ul)}), Vertex<2>{1.0,1.0});//1 3
    DBGVAR(cache);

    // std::vector<decltype(cache)::iterator> rmVert;
    // for (auto x = cache.begin(); x != cache.end(); ++x) {
    //     auto vert = x->second;
    //     double minDist = (vertices[0] - vert).normEuclid();
    //     size_t argMin = 0;
    //     for (size_t i = 1; i < vertices.size(); ++i) {
    //         auto dist = (vert - vertices[i]).normEuclid();
    //         if (dist < minDist) {
    //             minDist = dist;
    //             argMin = i;
    //         }
    //     }
    //     if (std::find(x->first.begin(), x->first.end(), argMin) == x->first.end()) {
    //         rmVert.emplace_back(x);
    //     }
    // }
    // for (auto x: rmVert) {
    //     DBGVAR(*x);
    //     cache.erase(x);
    // }

    auto mesh = buildMeshFromCache(cache);
    mesh.write("test-mesh.vtk");
}

void testView() {
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

       {0.2, 0.9},
    };
   MeshGenerationProblem<2> mgp;
   std::string meshPath = "../Meshes/mesh2D.vtk";    // 2D version
   auto compData = mgp.loadMesh(meshPath);
   std::string outPath = "../out";

   mgp.calculateFx(compData, vertices);
   DBGMSG("exporting mesh");
   mgp.exportMeshAndData(compData, outPath + "/MeshGenerator.vtk");
}

void testCache() {
    auto vCache = Cache();
    // { [ 3, 4, 9, 10 ]: [ 0.14, 0.68 ]}, { [ 4, 5, 10, 11 ]: [ 0.14, 0.86 ]} ]
    vCache.insert({0, 1, 6, 7}, {0.14, 0.14});
    vCache.insert({1, 2, 7, 8}, {0.14, 0.32});
    vCache.insert({2, 3, 8, 9}, {0.14, 0.5});
    auto val = vCache.find({0, 1, 6, 7});
    DBGVAR(*val, vCache);
}

int main(int argc, char** argv) {
    testGenerate();
    // testCache();
}
