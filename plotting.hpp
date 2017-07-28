//
// Created by Simone Raimondi on 22.07.17.
//

#ifndef SDF_PLOTTING_HPP
#define SDF_PLOTTING_HPP

#include <vector>
#include "grid.hpp"
#include "lodepng.hpp"

template<int W, int H>
void CreatePNG(const Grid<W, H> &grid,
               const std::string &file_name, int x_res, int y_res,
               const std::vector<float> &isolines, float iso_tol = 0.01f) {
    // Final vector containing the image data
    std::vector<unsigned char> image(x_res * y_res * 4, 0);

    // Initialise alpha channel
    for (int i = 3; i < x_res * y_res * 4; i += 4) {
        image[i] = 255;
    }

    const float x_pixel_delta = (grid.XMax() - grid.XMin()) / (float) x_res;
    const float y_pixel_delta = (grid.YMax() - grid.YMin()) / (float) y_res;

    // Go through the image pixels and interpolate
    for (int j = y_res - 1; j >= 0; j--) {
        for (int i = 0; i < x_res; i++) {
            // Compute position in grid
            const float x_pos = grid.XMin() + i * x_pixel_delta;
            const float y_pos = grid.YMin() + j * y_pixel_delta;

            // Compute grid value
            const float c = grid.At(x_pos, y_pos);

            const int image_index = j * x_res + i;
            // Plot isoline at 0 in red
            if (std::abs(c) < iso_tol) {
                image[4 * image_index] = 255;
            } else {
                for (auto const iso : isolines) {
                    if (std::abs(c - iso) < iso_tol) {
                        image[4 * image_index] = 255;
                        image[4 * image_index + 1] = 255;
                    }
                }
            }
        }
    }

    // Encode the image
    const unsigned error = lodepng::encode(file_name, image, x_res, y_res);
    // Check if there is an error
    if (error) {
        std::cerr << "Encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
    }
}


#endif //SDF_PLOTTING_HPP
