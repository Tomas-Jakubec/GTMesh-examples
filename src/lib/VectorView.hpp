#pragma once

#include <vector>

template<typename T = int>
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

    template<typename TData>
    static VectorView<TData> of(std::vector<TData>& data) {
        std::vector<size_t> indices(data.size());
        for (size_t i = 0; i < indices.size(); ++i) {
            indices[i] = i;
        }
        return VectorView<TData>(data, indices);
    }
};
