//
// Created by Simon on 27.07.2017.
//

#ifndef SDF_GRID_HPP
#define SDF_GRID_HPP

#include <memory>
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include "common.hpp"

// Define 2d Grid class
template<int W, int H>
class Grid {
private:
    // Grid points
    std::unique_ptr<float[]> data;
    // Grid bounds
    float x_min, x_max, y_min, y_max;
    // Delta between points
    float delta_x, delta_y;

    // Convert 2D index to linear index, using periodic boundaries
    inline int LinearIndex(int i, int j) const {
        if (i == W) { i = 0; }
        else if (i == -1) { i = W - 1; }

        if (j == H) { j = 0; }
        else if (j == -1) { j = H - 1; }

        return (j * W + i);
    }

public:
    // Default constructor
    Grid(float x_min, float x_max, float y_min, float y_max);

    // Element access operator
    inline float operator()(int i, int j) const {
        return data[LinearIndex(i, j)];
    }

    inline float &operator()(int i, int j) {
        return data[LinearIndex(i, j)];
    }

    // Get grid dimensions
    inline constexpr int Width() const { return W; }

    inline constexpr int Height() const { return H; }

    // Get boundaries
    inline float XMin() const { return x_min; }

    inline float XMax() const { return x_max; }

    inline float YMin() const { return y_min; }

    inline float YMax() const { return y_max; }

    inline float DeltaX() const { return delta_x; }

    inline float DeltaY() const { return delta_y; }

    // Get maximum and minimum
    inline float Max() const {
        float max = -INFINITY;
        for (int i = 0; i < W * H; i++) {
            if (data[i] > max) { max = data[i]; }
        }
        return max;
    }

    inline float Min() const {
        float min = INFINITY;
        for (int i = 0; i < W * H; i++) {
            if (data[i] < min) { min = data[i]; }
        }
        return min;
    }

    // Compute value at given point (x,y) using bilinear interpolation
    inline float At(float x, float y) const {
        // Find indices
        const int x_i = static_cast<int>((x - x_min) / delta_x);
        const int y_i = static_cast<int>((y - y_min) / delta_y);

        // Interpolate along x
        const float t_x = x - std::floor(x);
        const float t_y = y - std::floor(y);
        const float c0 = Lerp(t_x, this->operator()(x_i, y_i), this->operator()(x_i + 1, y_i));
        const float c1 = Lerp(t_x, this->operator()(x_i, y_i + 1), this->operator()(x_i + 1, y_i + 1));

        return Lerp(t_y, c0, c1);
    }

    // Compute coordinates of point given indices
    inline void CoordsAt(int i, int j, float *x, float *y) const {
        *x = x_min + delta_x * i;
        *y = y_min + delta_y * j;
    }

    // Print grid to file
    void Print(const std::string &file_name) const;
};

template<int W, int H>
Grid<W, H>::Grid(float x_min, float x_max, float y_min, float y_max)
        : data(new float[W * H]),
          x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max),
          delta_x((x_max - x_min) / static_cast<float>(W - 1)),
          delta_y((y_max - y_min) / static_cast<float>(H - 1)) {
    for (int j = 0; j < H; j++) {
        for (int i = 0; i < W; i++) {
            this->operator()(i, j) = 0.f;
        }
    }
}

template<int W, int H>
void Grid<W, H>::Print(const std::string &file_name) const {
    // Create and open file
    std::ofstream file;
    file.open(file_name);

    file << "[";
    for (int j = 0; j < H; j++) {
        for (int i = 0; i < W; i++) {
            file << this->operator()(i, j);
            if ((i * j) != ((H - 1) * (W - 1))) { file << ","; }
        }
    }
    file << "]" << std::endl;

    // Close file
    file.close();
}

#endif //SDF_GRID_HPP
