//
// Created by Simone Raimondi on 22.07.17.
//

#ifndef SDF_SDF_HPP
#define SDF_SDF_HPP

#include <cmath>
#include "grid.hpp"

// Compute squared magnitude of the gradient at given point
template<int W, int H>
float GradNorm2(const Grid<W, H> &sdf_grid, int i, int j) {
    // Compute derivatives using finite difference
    const float grad_x = (sdf_grid(i + 1, j) - sdf_grid(i, j)) / sdf_grid.DeltaX();
    const float grad_y = (sdf_grid(i, j + 1) - sdf_grid(i, j)) / sdf_grid.DeltaY();

    return grad_x * grad_x + grad_y * grad_y;
}

// Initialise the sign grid
template<int W, int H>
void InitialiseSignGrid(const Grid<W, H> &sdf_grid, Grid<W, H> &sign_grid) {
    // This procedure uses formula 7.5 of the Level Set book to initialise the sign grid used after in the PDE solving
    for (int j = 0; j < sign_grid.Height(); j++) {
        for (int i = 0; i < sign_grid.Width(); i++) {
            const float sdf_v = sdf_grid(i, j);
            sign_grid(i, j) = sdf_v / std::sqrt(sdf_v * sdf_v + sdf_grid.DeltaX() * sdf_grid.DeltaX());
        }
    }
}

// Compute the sign grid
template<int W, int H>
void ComputeSignGrid(const Grid<W, H> &sdf_grid, Grid<W, H> &sign_grid) {
    // This procedure uses formula 7.6 of the Level Set book to initialise the sign grid used after in the PDE solving
    for (int j = 0; j < sign_grid.Height(); j++) {
        for (int i = 0; i < sign_grid.Width(); i++) {
            // Get SDF value
            const float sdf_v = sdf_grid(i, j);
            // Gradient squared norm
            const float grad_norm_sq = GradNorm2(sdf_grid, i, j);
            // Compute sign
            sign_grid(i, j) = sdf_v / std::sqrt(sdf_v * sdf_v + grad_norm_sq * sdf_grid.DeltaX() * sdf_grid.DeltaX());
        }
    }
}

// Compute the right term of the reinitialisation PDE
template<int W, int H>
float RightTermPDE(const Grid<W, H> &sdf_grid, const Grid<W, H> &sign_grid, int i, int j) {
    // Compute gradient norm term using first order Upwind
    const float dphi_dx_plus = (sdf_grid(i + 1, j) - sdf_grid(i, j)) / sdf_grid.DeltaX();
    const float dphi_dx_minus = (sdf_grid(i, j) - sdf_grid(i - 1, j)) / sdf_grid.DeltaX();

    const float dphi_dy_plus = (sdf_grid(i, j + 1) - sdf_grid(i, j)) / sdf_grid.DeltaY();
    const float dphi_dy_minus = (sdf_grid(i, j) - sdf_grid(i, j - 1)) / sdf_grid.DeltaY();

    // Check sign and compute gradient norm
    float grad_x_sq, grad_y_sq;
    if (sign_grid(i, j) > 0.f) {
        grad_x_sq = Max(PositiveSQ(dphi_dx_minus), NegativeSQ(dphi_dx_plus));
        grad_y_sq = Max(PositiveSQ(dphi_dy_minus), NegativeSQ(dphi_dy_plus));
    } else {
        grad_x_sq = Max(PositiveSQ(dphi_dx_plus), NegativeSQ(dphi_dx_minus));
        grad_y_sq = Max(PositiveSQ(dphi_dy_plus), NegativeSQ(dphi_dy_minus));
    }

    const float grad_norm = std::sqrt(grad_x_sq + grad_y_sq);

    return sign_grid(i, j) * (1.f - grad_norm);
}

// Construct SDF on given initial configuration
template<int W, int H>
void ConstructSDF(Grid<W, H> &sdf_grid, float band, float dt, int max_iters, float tolerance) {
    // Sign grid
    Grid<W, H> sign_grid(sdf_grid.XMin(), sdf_grid.XMax(), sdf_grid.YMin(), sdf_grid.YMax());
    InitialiseSignGrid(sdf_grid, sign_grid);

    // Grid storing dphi_dt
    Grid<W, H> sdf_dt(sdf_grid.XMin(), sdf_grid.XMax(), sdf_grid.YMin(), sdf_grid.YMax());
    float max_dt = tolerance + 1.f;
    int iters = 0;

    while (iters < max_iters && max_dt > tolerance) {
        // Reset max_dt
        max_dt = 0.f;

        // Loop over all grid
        for (int j = 0; j < H; j++) {
            for (int i = 0; i < W; i++) {
                if (std::abs(sdf_grid(i, j)) >= band) {
                    sdf_dt(i, j) = 0.f;
                } else {
                    // Compute right part of the PDE
                    sdf_dt(i, j) = RightTermPDE(sdf_grid, sign_grid, i, j);
                    // Update max dt
                    if (std::abs(sdf_dt(i, j)) > max_dt) { max_dt = std::abs(sdf_dt(i, j)); }
                }
            }
        }

        // Update SDF
        for (int j = 0; j < H; j++) {
            for (int i = 0; i < W; i++) {
                sdf_grid(i, j) += dt * sdf_dt(i, j);
            }
        }

        iters++;
    }

    std::cout << "Done constructing SDF in " << iters << " iterations." << std::endl;
}

// Construct SDF on given initial configuration
template<int W, int H>
void ConstructSDF2(Grid<W, H> &sdf_grid, float dt, int max_iters, float tolerance, float band) {
    // Initialise sign grid
    Grid<W, H> sign_grid(sdf_grid.XMin(), sdf_grid.XMax(), sdf_grid.YMin(), sdf_grid.YMax());
    ComputeSignGrid(sdf_grid, sign_grid);

    // Grid storing dphi_dt
    Grid<W, H> sdf_dt(sdf_grid.XMin(), sdf_grid.XMax(), sdf_grid.YMin(), sdf_grid.YMax());
    float max_dt = tolerance + 1.f;
    int iters = 0;

    while (iters < max_iters && max_dt > tolerance) {
        // Reset max_dt
        max_dt = 0.f;

        // Loop over all grid
        for (int j = 0; j < H; j++) {
            for (int i = 0; i < W; i++) {
                if (std::abs(sdf_grid(i, j)) >= band) {
                    sdf_dt(i, j) = 0.f;
                } else {
                    // Compute right part of the PDE
                    sdf_dt(i, j) = RightTermPDE(sdf_grid, sign_grid, i, j);
                    // Update max dt
                    if (std::abs(sdf_dt(i, j)) > max_dt) { max_dt = std::abs(sdf_dt(i, j)); }
                }
            }
        }

        // Update SDF
        for (int j = 0; j < H; j++) {
            for (int i = 0; i < W; i++) {
                sdf_grid(i, j) += dt * sdf_dt(i, j);
            }
        }

        // Update sign grid for the next iteration
        ComputeSignGrid(sdf_grid, sign_grid);

        iters++;
    }

    std::cout << "Done constructing SDF in " << iters << " iterations." << std::endl;
}

#endif //SDF_SDF_HPP
