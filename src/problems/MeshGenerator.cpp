#include <GTMesh/Traits/Traits.h>
#include <GTMesh/Traits/TraitsAlgorithm/TraitsAlgorithm.h>
#include <GTMesh/UnstructuredMesh/MeshDataContainer/MeshDataIO/VTKMeshDataWriter.h>
#include <GTMesh/UnstructuredMesh/UnstructuredMesh.h>
#include <assert.h>

#include <iomanip>
#include <lib/Gradient.hpp>
#include <lib/VectorView.hpp>
#include <lib/VirtualGrid.hpp>
#include <lib/optimisation.hpp>
#include <limits>
#include <list>
#include <memory>
#include <random>

#include "solvers/RKMSolver.hpp"

// struct ComputationData
// {
//     double minf;
//     // double F;
//     double tildeF;
//     size_t cellIndex;
// };

// MAKE_NAMED_ATTRIBUTE_TRAIT(ComputationData, "f_{min}", minf, "F^{~}=",
// tildeF, "cellIndex", cellIndex);

/**
 * Calculates value represending the distance from x and the equidistant vertex
 * of vertices.
 */
template<unsigned int Dimension, typename RealType>
RealType calculateScoreSum(const VectorView<Vertex<Dimension, RealType> >& vertices, const Vertex<Dimension, RealType>& x) {
    RealType result = 0.0;
    auto f_x0 = (x - vertices[0]).normEuclid();
    auto maxF = f_x0;
    auto minF = f_x0;

    for (size_t i = 1; i < vertices.size(); ++i) {
        auto f_x = (x - vertices[i]).normEuclid();
        if (f_x > maxF) {
            maxF = f_x;
        }
        if (f_x < minF) {
            minF = f_x;
        }
    }
    return pow(maxF - minF, 2) * maxF;
}

template<unsigned int Dimension, typename RealType>
auto calculateNeighbors(const VectorView<Vertex<Dimension, RealType> >& vertices, const Vertex<Dimension, RealType>& x, const RealType treshold) {
    auto result = std::vector<size_t>();
    result.reserve(5);
    std::vector<double> f_xs(vertices.size());
    double minF = (x - vertices[0]).normEuclid();
    f_xs[0] = minF;
    for (size_t i = 1; i < vertices.size(); ++i) {
        auto f_x = (x - vertices[i]).normEuclid();
        f_xs[i] = f_x;

        if (f_x < minF) {
            minF = f_x;
        }
    }
    for (size_t i = 0; i < vertices.size(); ++i) {
        if (f_xs[i] - minF <= treshold) {
            result.emplace_back(i);
        }
    }
    return result;
}

// template<unsigned int ProblemDimension>
// struct MeshGenerationProblem{
//     using MeshType = UnstructuredMesh<ProblemDimension, size_t, double>;
//     using ResultType = ComputationData;
//     using ProblemDataContainerType = MeshDataContainer<ResultType,
//     ProblemDimension>; std::shared_ptr<MeshReader<ProblemDimension>>
//     meshReader; VTKMeshWriter<ProblemDimension, size_t, double> meshWriter;

//     MeshType mesh;

//     void calculateFx(ProblemDataContainerType& data,
//     std::vector<Vertex<ProblemDimension, double>> vertices) {
//         std::vector<double> f_xs(vertices.size());
//         double treshold = 0.025;
//         auto verticesView = VectorView<>::of(vertices);
//         for (const auto& cell : mesh.getCells()) {
//             auto neighborhood = calculateNeighbors(verticesView,
//             cell.getCenter(), treshold); data[cell].tildeF =
//             calculateScoreSum(verticesView.slice(neighborhood),
//             cell.getCenter());
//         }
//     }

//     ProblemDataContainerType loadMesh(const std::string& meshPath){
//         meshReader = mesh.load(meshPath);
//         mesh.initializeCenters();
//         mesh.setupBoundaryCells();
//         mesh.setupBoundaryCellsCenters();
//         // meshData.allocateData(mesh);
//         auto measures = mesh.template
//         computeElementMeasures<METHOD_TESSELLATED>(); auto cellsDist =
//         computeCellsDistance(mesh);
// //        for (const auto& cell : mesh.getCells()){
// //            meshData[cell].invCellVolume = 1.0 / measures[cell];
// //        }
// //        for (const auto& face : mesh.getFaces()){
// //            meshData[face].measure = measures[face];
// //            meshData[face].measureOverCellsDistance = measures[face] /
// cellsDist[face];
// //        }
//         return ProblemDataContainerType(mesh, ComputationData{0,0});
//     }

//     void exportMeshAndData(ProblemDataContainerType& compData,
//                            const std::string& outputPath){
//         std::ofstream ofst(outputPath);
//         meshWriter.writeHeader(ofst, "heat-conduction");
//         meshWriter.writeToStream(ofst, mesh,meshReader->getCellTypes());
//         VTKMeshDataWriter<ProblemDimension> dataWriter;
//         dataWriter.writeToStream(ofst, compData, meshWriter);
//     }
// };

namespace std {
template<unsigned int Dimension, typename T> class hash<Vertex<Dimension, T> > {
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

template<typename T> class hash<std::vector<T> > {
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

template<typename T, size_t Size> class hash<std::array<T, Size> > {
public:
    size_t operator()(const std::array<T, Size>& vec) const noexcept {
        const std::hash<T> hashT{};
        size_t hash = 0;
        for (size_t i = 0; i < vec.size(); ++i) {
            hash ^= hashT(vec[i]) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};
}  // namespace std

template<typename T> bool isSubset(const std::vector<T>& set1, const std::vector<T>& set2) {
    std::vector<T> intersection;
    std::set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(), std::back_inserter(intersection));
    return intersection.size() == set1.size();
}

template<typename KeyType, typename ValueType>
ValueType& getOrPut(std::unordered_map<KeyType, ValueType>& map, const KeyType& key, const ValueType& value) {
    auto it = map.find(key);
    if (it == map.end()) {
        // Key not found, insert the key-value pair
        it = map.emplace(key, value).first;
    }
    return it->second;
}

/**
 *
 */
template<unsigned int Dim> class Cache {
    using cacheElementT = std::pair<std::vector<std::size_t>, Vertex<Dim> >;
    std::vector<cacheElementT> cacheData;
    std::unordered_map<std::size_t, std::set<std::size_t> > cacheMap;
    static constexpr auto invalidOffset = std::numeric_limits<size_t>::max();

public:
    void insert(const std::vector<std::size_t>& indices, Vertex<Dim> vertex) {
        auto newIndex = cacheData.size();
        cacheData.emplace_back(std::make_pair(indices, vertex));
        for (const auto index : indices) {
            auto& cachedItems = getOrPut(cacheMap, index, std::set<std::size_t>());
            cachedItems.emplace(newIndex);
        }
    }

    void insert(std::vector<std::size_t>&& indices, Vertex<Dim> vertex) {
        auto newIndex = cacheData.size();
        cacheData.emplace_back(std::make_pair(indices, vertex));
        for (const auto index : indices) {
            auto& cachedItems = getOrPut(cacheMap, index, std::set<std::size_t>());
            cachedItems.emplace(newIndex);
        }
    }

    void insertInvalid(std::vector<std::size_t>&& indices) {
        for (const auto index : indices) {
            auto& cachedItems = getOrPut(cacheMap, index, std::set<std::size_t>());
            cachedItems.emplace(invalidOffset);
        }
    }

    struct FindResult {
        const typename std::vector<cacheElementT>::const_iterator iter;
        bool valid;
    };

    const auto find(const std::vector<std::size_t>& indices) const {
        std::set<std::size_t> result;
        auto cachedSets = std::vector<std::pair<std::set<size_t>::iterator, std::set<size_t>::iterator> >();
        cachedSets.reserve(indices.size());

        for (std::size_t index : indices) {
            auto indexValIter = cacheMap.find(index);
            if (indexValIter != cacheMap.end()) {
                auto& cachedOffests = indexValIter->second;
                cachedSets.emplace_back(std::make_pair(cachedOffests.begin(), cachedOffests.end()));
            } else {
                return FindResult{cacheData.end(), true};
            }
        }

        bool process = true;
        // calculate intersection of multiple sets at once
        while (process) {
            auto minValue = std::numeric_limits<size_t>::max();
            for (auto& range : cachedSets) {
                if (minValue > *range.first) {
                    minValue = *range.first;
                }
            }
            bool addValue = true;
            for (auto& range : cachedSets) {
                if (minValue != *range.first) {
                    addValue = false;
                } else {
                    ++range.first;
                    if (range.first == range.second) {
                        process = false;
                    }
                }
            }
            if (addValue) {
                result.insert(minValue);
            }
        }
        if (result.size() == 0) {
            return FindResult{cacheData.end(), true};
        } else if (*result.begin() == invalidOffset) {
            return FindResult{cacheData.end(), false};
        } else {
            return FindResult{cacheData.begin() + *result.begin(), true};
        }
    }

    void consolide() {
        Cache<Dim> result;
        std::sort(cacheData.begin(), cacheData.end(),
                  [](const cacheElementT& x1, const cacheElementT& x2) { return x1.first.size() > x2.first.size(); });
        for (const cacheElementT& value : cacheData) {
            if (!result.contains(value.first)) {
                result.insert(value.first, value.second);
            }
        }
        cacheData = std::move(result.cacheData);
        cacheMap = std::move(result.cacheMap);
    }

    bool contains(const std::vector<std::size_t>& indices) const {
        auto cachedSets = std::vector<std::pair<std::set<size_t>::iterator, std::set<size_t>::iterator> >();
        cachedSets.reserve(indices.size());

        for (std::size_t index : indices) {
            auto indexValIter = cacheMap.find(index);
            if (indexValIter != cacheMap.end()) {
                auto& cachedOffests = indexValIter->second;
                cachedSets.emplace_back(std::make_pair(cachedOffests.begin(), cachedOffests.end()));
            } else {
                return false;
            }
        }

        bool process = true;
        // calculate intersection of multiple sets at once
        while (process) {
            auto minValue = std::numeric_limits<size_t>::max();
            for (auto& range : cachedSets) {
                if (minValue > *range.first) {
                    minValue = *range.first;
                }
            }
            bool addValue = true;
            for (auto& range : cachedSets) {
                if (minValue != *range.first) {
                    addValue = false;
                } else {
                    ++range.first;
                    if (range.first == range.second) {
                        process = false;
                    }
                }
            }
            if (addValue) {
                return true;
            }
        }
        return false;
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
        auto findResult = find(indices);
        if (findResult.valid || findResult.iter == end()) {
            throw std::out_of_range("Indices not found in the cache!");
        }
        return *findResult.iter;
    }
};

Cache<2> cache;

template<typename T> std::vector<T> cat(std::vector<T> v1, const std::vector<T>& v2) {
    v1.insert(v1.end(), v2.begin(), v2.end());
    return v1;
}

template<unsigned int Dim, typename TArea>
void scanArea(const VirtualGrid<Dim>& grid, const VectorView<Vertex<Dim> >& vertices, const TArea& area, const double maxPrecision = 1e-2) {
    scanArea(
        grid, vertices, [](const Vertex<Dim>& x) { return x; }, area, {}, maxPrecision);
}

/**
 * @param virtualBoundaryCells vector of virtual cell indices that are added to
 * all identified vertices. Each manifold of Dim -1 has virtual cell assigned.
 * Virtual cells identifying boundaries of lower dimensions are derived from the
 * neighboring boundaries.
 */
template<unsigned int GridDim, unsigned int ProblemDimension, typename F, typename TArea>
void scanArea(const VirtualGrid<GridDim>& grid, const VectorView<Vertex<ProblemDimension> >& vertices, F&& boundaryMapping, const TArea& area,
              std::vector<size_t> virtualBoundaryCells, const double maxPrecision = 1e-2) {
    double treshold = 2.0 * std::pow(double(GridDim), 0.5) * grid.stepSize;
    for (auto vertex : grid) {
        auto neighboringCells = calculateNeighbors(vertices, boundaryMapping(vertex), treshold);
        if (neighboringCells.size() > GridDim &&
            !cache.contains(virtualBoundaryCells.empty() ? neighboringCells : cat(neighboringCells, virtualBoundaryCells))) {
            if (neighboringCells.size() == GridDim + 1 || grid.stepSize < maxPrecision) {
                auto affectedVertices = vertices.slice(neighboringCells);
                auto lossFunction = [&boundaryMapping, &affectedVertices](const Vertex<GridDim>& x) {
                    return calculateScoreSum(affectedVertices, boundaryMapping(x));
                };
                auto optimisedVertex = boundaryMapping(optimisation::BFGS(vertex, lossFunction, maxPrecision, maxPrecision * 1e-3));
                if (!area.contains(optimisedVertex)) {
                    continue;
                }
                // Calculate the real neighborhood
                auto neighbors = calculateNeighbors(vertices, optimisedVertex, maxPrecision);
                if ((optimisedVertex - boundaryMapping(vertex)).normEuclid() <= treshold && isSubset(neighboringCells, neighbors)) {
                    cache.insert(std::move(virtualBoundaryCells.empty() ? neighbors : cat(neighbors, virtualBoundaryCells)), optimisedVertex);
                }
            } else {
                auto refinedGrid = VirtualGrid<GridDim>::withCenterOrigin(vertex, 10, 0.5 * grid.stepSize);
                scanArea(refinedGrid, vertices, boundaryMapping, area, virtualBoundaryCells, maxPrecision);
            }
        }
    }
}

template<typename IndexType> class Histogram {
public:
    void update(const std::vector<IndexType>& numbers) {
        // Count occurrences of each number
        for (IndexType num : numbers) {
            auto& auxCount = getOrPut(counts, num, IndexType());
            auxCount++;
        }
    }

    const std::unordered_map<IndexType, IndexType>& getHistogram() const {
        return counts;
    }

private:
    std::unordered_map<IndexType, IndexType> counts;
};

auto buildMeshFromCache(Cache<2>& cache) {
    // cell index -> { neigboring cell index -> vertex index (they have in common) }
    std::unordered_map<size_t, std::unordered_map<size_t, std::vector<std::size_t> > > buildCache;
    std::vector<Vertex<2> > vertices;
    // vertices.reserve(cache.cacheData.size());

    for (auto item : cache) {
        vertices.emplace_back(item.second);
        for (const auto index : item.first) {
            // get or put empty value
            auto cachedData = buildCache.find(index);
            if (cachedData == buildCache.end()) {
                cachedData = buildCache.emplace(index, std::unordered_map<size_t, std::vector<std::size_t> >()).first;
            }

            auto& cachedNeighbors = cachedData->second;

            for (const auto index2 : item.first) {
                if (index2 == index) {
                    continue;
                }
                // get or put
                auto cachedNeighborIter = cachedNeighbors.find(index2);
                if (cachedNeighborIter == cachedNeighbors.end()) {
                    cachedNeighborIter = cachedNeighbors.emplace(index2, std::vector<std::size_t>()).first;
                }
                cachedNeighborIter->second.emplace_back(vertices.size() - 1);
            }
        }
    }
    DBGVAR(buildCache);
    using MeshType = UnstructuredMesh<2, size_t, double>;
    MeshType mesh;
    for (auto vert : vertices) {
        mesh.getVertices().emplace_back(MeshType::Vertex(mesh.getVertices().size(), vert));
    }
    std::unordered_map<std::vector<size_t>, size_t> edgeCache;
    for (const auto& cellCache : buildCache) {
        if (isBoundaryIndex(cellCache.first)) {
            continue;
        }
        auto cellIndex = mesh.getCells().size();
        mesh.getCells().push_back(cellIndex);
        auto& cell = mesh.getCells().back();
        size_t lastEdgeIndex = INVALID_INDEX(size_t);
        size_t firstEdgeIndex;
        for (const auto& cellBoundary : cellCache.second) {
            if (cellBoundary.second.size() >= 2) {
                auto cachedEdge = edgeCache.find(cellBoundary.second);
                size_t newIndex;
                if (cachedEdge != edgeCache.end()) {
                    newIndex = cachedEdge->second;
                } else {
                    newIndex = mesh.getEdges().size();
                    mesh.getEdges().emplace_back(MeshType::Edge(newIndex, cellBoundary.second[0], cellBoundary.second[1]));
                    edgeCache.emplace(cellBoundary.second, newIndex);
                }
                auto& edge = mesh.getEdges()[newIndex];
                if (lastEdgeIndex == INVALID_INDEX(size_t)) {
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

std::vector<Vertex<2> > generateHexagonalGridVertices(int nElements) {
    std::vector<Vertex<2> > vertices;

    auto diameter = 1.0;
    auto stepSize = diameter / (nElements - 1);
    for (auto vertex : VirtualGrid<2>({0.0, 0.0}, nElements, stepSize)) {
        vertices.emplace_back(std::move(vertex));
    }

    for (auto vertex : VirtualGrid<2>({(stepSize * 0.5), (stepSize * 0.5)}, nElements - 1, stepSize)) {
        vertices.emplace_back(std::move(vertex));
    }

    return vertices;
}

std::vector<Vertex<2> > generate_hexagonal_grid(int rows, int columns) {
    std::vector<Vertex<2> > grid;

    auto horizontal_distance_between_hexagons = 2.0 / (2.0 * columns - 1);
    auto vertical_distance_between_hexagons = 1.0 / (rows - 1);

    for (int row = 0; row < rows; ++row) {
        for (int column = 0; column < columns; ++column) {
            double x = column * horizontal_distance_between_hexagons;
            double y = row * vertical_distance_between_hexagons;

            if (row % 2 == 1) {
                x += horizontal_distance_between_hexagons * 0.5;
            }
            Vertex<2> hexagon = Vertex<2>{x, y};
            grid.push_back(hexagon);
        }
    }

    return grid;
}

std::vector<Vertex<2> > generateRandomVertices(int nElements) {
    std::vector<Vertex<2> > vertices;
    std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist(0, 1);
    for (int i = 0; i < nElements; ++i) {
        vertices.emplace_back(Vector<2>{dist(e2), dist(e2)});
    }

    return vertices;
}

class Square {
private:
    Vertex<2, double> center;
    double sideLength;
    Vertex<2, double> min;
    Vertex<2, double> max;

public:
    // Constructor
    Square(Vertex<2, double> center, double sideLength) : center(center), sideLength(sideLength) {
        double halfSide = sideLength / 2.0;
        min = center - Vertex<2, double>{halfSide, halfSide};
        max = center + Vertex<2, double>{halfSide, halfSide};
    }

    // Method to check if a point is inside the square
    bool contains(const Vertex<2, double>& point) const {
        return (point[0] >= min[0] && point[0] <= max[0] && point[1] >= min[1] && point[1] <= max[1]);
    }
};

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

    auto vertices = generate_hexagonal_grid(50, 50);

    VirtualGrid<2> grid = VirtualGrid<2>::withCenterOrigin({0.5, 0.5}, 100, 0.5);
    auto area = Square({0.5, 0.5}, 1.0);
    auto verticesView = VectorView<>::of(vertices);
    DBGVAR(verticesView);
    scanArea(grid, verticesView, area, 1e-3);
    DBGVAR(cache);
    auto gridBoundary = VirtualGrid<1>::withCenterOrigin({0.5}, 100, 0.5);
    auto b1 = [](const Vertex<1>& v) { return Vertex<2>{v[0], 0.0}; };
    auto b2 = [](const Vertex<1>& v) { return Vertex<2>{v[0], 1.0}; };

    auto b3 = [](const Vertex<1>& v) { return Vertex<2>{0.0, v[0]}; };

    auto b4 = [](const Vertex<1>& v) { return Vertex<2>{1.0, v[0]}; };
    scanArea(gridBoundary, verticesView, b1, area, {makeBoundaryIndex(0ul)});
    scanArea(gridBoundary, verticesView, b2, area, {makeBoundaryIndex(1ul)});
    scanArea(gridBoundary, verticesView, b3, area, {makeBoundaryIndex(2ul)});
    scanArea(gridBoundary, verticesView, b4, area, {makeBoundaryIndex(3ul)});
    cache.insert(cat(calculateNeighbors(verticesView, Vertex<2>{0.0, 0.0}, 1e-4), {makeBoundaryIndex(0ul), makeBoundaryIndex(2ul)}),
                 Vertex<2>{0.0, 0.0});  // 0 2
    cache.insert(cat(calculateNeighbors(verticesView, Vertex<2>{1.0, 0.0}, 1e-4), {makeBoundaryIndex(0ul), makeBoundaryIndex(3ul)}),
                 Vertex<2>{1.0, 0.0});  // 0 3
    cache.insert(cat(calculateNeighbors(verticesView, Vertex<2>{0.0, 1.0}, 1e-4), {makeBoundaryIndex(1ul), makeBoundaryIndex(2ul)}),
                 Vertex<2>{0.0, 1.0});  // 1 2
    cache.insert(cat(calculateNeighbors(verticesView, Vertex<2>{1.0, 1.0}, 1e-4), {makeBoundaryIndex(1ul), makeBoundaryIndex(3ul)}),
                 Vertex<2>{1.0, 1.0});  // 1 3
    cache.consolide();
    DBGVAR(cache);

    auto mesh = buildMeshFromCache(cache);
    mesh.write("test-mesh.vtk");
}

// void testView() {
//     std::vector<Vertex<2>> vertices = {
//         {0.5, 0.5},
//         {0.75, 0.5},
//         {0.5, 0.75},
//         {0.25, 0.5},
//         {0.5, 0.25},

//         {0.9, 0.1},
//         {0.9, 0.9},
//         {0.1, 0.1},
//         {0.1, 0.9},

//        {0.2, 0.9},
//     };
//    MeshGenerationProblem<2> mgp;
//    std::string meshPath = "../Meshes/mesh2D.vtk";    // 2D version
//    auto compData = mgp.loadMesh(meshPath);
//    std::string outPath = "../out";

//    mgp.calculateFx(compData, vertices);
//    DBGMSG("exporting mesh");
//    mgp.exportMeshAndData(compData, outPath + "/MeshGenerator.vtk");
// }

void testCache() {
    auto vCache = Cache<2>();
    // { [ 3, 4, 9, 10 ]: [ 0.14, 0.68 ]}, { [ 4, 5, 10, 11 ]: [ 0.14, 0.86 ]} ]
    vCache.insert({0, 1, 6, 7}, {0.14, 0.14});
    vCache.insert({1, 2, 7, 8}, {0.14, 0.32});
    vCache.insert({2, 3, 8, 9}, {0.14, 0.5});
    auto val = vCache.find({0, 1, 6, 7});
    DBGVAR(*val.iter, vCache);
}

void testGradientCalculation() {
    Vertex<2> x{1.0, 2.0};
    auto f = [](const Vector<2>& x) { return multiply(x, x); };
    auto gradF = gradientOf(f, 1e-5);
    auto grad = gradF(x);
    DBGVAR(grad, (e<2, 0>()), (e<2, 1>()));
}

int main(int argc, char** argv) {
    testGenerate();
    // testGradientCalculation();
    // testCache();
}
