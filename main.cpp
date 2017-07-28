#include <iostream>
#include "plotting.hpp"
#include "sdf.hpp"

const int WIDTH = 512;
const int HEIGHT = 512;

int main() {
    // Create grid
    Grid<WIDTH, HEIGHT> grid_ellipse(-5, 5, -5, 5);
    Grid<WIDTH, HEIGHT> grid_circle(-5, 5, -5, 5);

    // Initialize grid using Ellipse function and circle
    const float a = 2.f;
    const float b = 1.f;
    const float radius = 1.f;
    // Initialise grid
    for (int j = 0; j < HEIGHT; j++) {
        for (int i = 0; i < WIDTH; i++) {
            float x, y;
            // Compute coordinates of the point
            grid_ellipse.CoordsAt(i, j, &x, &y);

            // Test with a function that is not a SDF
            grid_ellipse(i, j) = (x * x) / (a * a) + (y * y) / (b * b) - 1.f;

            // Test with a SDF
            grid_circle(i, j) = std::sqrt(x * x + y * y) - radius;
        }
    }

    // Isolines
    std::vector<float> isolines({-1.f, -0.5f, 0.5f, 1.f, 1.5f, 2.f, 2.5f, 3.f, 3.5f, 4.f, 4.5f, 5.f});

    // Create starting images
    CreatePNG(grid_ellipse, "start_ellipse.png", 512, 512, isolines, 0.02f);
    CreatePNG(grid_circle, "start_circle.png", 512, 512, isolines);

    // Try to construct SDF using formula 7.6
    ConstructSDF2(grid_ellipse, 0.001f, 20000, 0.01f, 100.f);
    CreatePNG(grid_ellipse, "final_ellipse.png", 512, 512, isolines);

    ConstructSDF2(grid_circle, 0.001f, 20000, 0.01f, 100.f);
    CreatePNG(grid_circle, "final_circle.png", 512, 512, isolines);
}