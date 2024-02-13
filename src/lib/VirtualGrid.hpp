#pragma once

#include <GTMesh/NumericStaticArray/Vector.h>

template<unsigned int Dimension, typename IndexType = size_t, typename RealType = double>
class VirtualGrid {
public:
    using VertexType = Vertex<Dimension, RealType>;
    const VertexType origin;

    class Iterator {
        const VirtualGrid& virtualGrid;
        std::array<IndexType, Dimension> indices;
        public:

        Iterator(const VirtualGrid& virtualGrid, std::array<IndexType, Dimension> indices = {}): virtualGrid(virtualGrid), indices(indices) {}

        Iterator operator++() {
            indices[0]++;
            for(size_t i = 0; i < Dimension - 1; ++i){
                if(indices[i] == virtualGrid.nElements) {
                    indices[i] = 0;
                    indices[i+1]++;
                }
            }
            return *this;
        }

        VertexType operator*(){
            return virtualGrid.at(indices);
        }

        bool operator==(const Iterator& other) {
            return indices == other.indices;
        }

        bool operator!=(const Iterator& other) {
            return indices != other.indices;
        }
    };

    const IndexType nElements;
    const RealType stepSize;
    VirtualGrid(const VertexType& origin, IndexType nElements, RealType stepSize): origin(origin), nElements(nElements), stepSize(stepSize) {}

    static VirtualGrid withCenterOrigin(const VertexType& origin, IndexType nElements, RealType radius) {
        VertexType shift;
        for(auto & coordinate: shift) {
            coordinate = radius;
        }
        VertexType originShifted = origin - shift;
        auto stepSize = 2 * radius / nElements;
        return VirtualGrid(originShifted, nElements + 1, stepSize);
    }

    static VirtualGrid withCenterOriginStep(const VertexType& origin, IndexType nElements, RealType stepSize) {
        VertexType shift;
        auto radius = 0.5 * stepSize * (nElements - 1);
        for(auto & coordinate: shift) {
            coordinate = radius;
        }
        VertexType originShifted = origin - shift;
        return VirtualGrid(originShifted, nElements, stepSize);
    }

    Iterator begin() const {
        return Iterator(*this);
    }

    Iterator end() const {
        std::array<IndexType, Dimension> indices{};
        indices.back() = nElements;
        return Iterator(*this, indices);
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
