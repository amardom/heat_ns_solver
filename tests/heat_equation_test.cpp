#include "solver.hpp"

#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>

double sine_decay_error(int nx, int ny) {
    Config config;
    config.nx = nx;
    config.ny = ny;
    config.length = 2.0;
    config.height = 1.0;
    config.diffusivity = 0.1;
    config.dt = 0.0001;
    config.steps = 1000;
    config.fluid = false;
    config.initial_temperature = 0.0;
    config.left = {0.0, ThermalBoundaryType::temperature, 0.0};
    config.right = config.left;
    config.bottom = config.left;
    config.top = config.left;
    Solver solver(config);
    double pi = std::acos(-1.0);
    for (int y = 1; y <= ny; ++y) {
        for (int x = 1; x <= nx; ++x) {
            solver.temperature(x, y) = std::sin(pi * (x - 0.5) * solver.dx / config.length) *
                                       std::sin(pi * (y - 0.5) * solver.dy / config.height);
        }
    }
    for (int step = 0; step < config.steps; ++step) {
        solver.step();
    }
    double time = config.steps * config.dt;
    double decay = std::exp(
        -config.diffusivity * pi * pi *
        (1 / (config.length * config.length) + 1 / (config.height * config.height)) * time);
    double error = 0;
    for (int y = 1; y <= ny; ++y) {
        for (int x = 1; x <= nx; ++x) {
            double exact = decay * std::sin(pi * (x - 0.5) * solver.dx / config.length) *
                           std::sin(pi * (y - 0.5) * solver.dy / config.height);
            error = std::max(error, std::abs(solver.temperature(x, y) - exact));
        }
    }
    return error;
}

TEST(HeatEquation, MatchesAnalyticalSineDecayOnRectangle) {
    EXPECT_LT(sine_decay_error(48, 24), 2e-4);
}

TEST(HeatEquation, HalvingCellSizeReducesSpatialErrorByFour) {
    double coarse_error = sine_decay_error(24, 12);
    double fine_error = sine_decay_error(48, 24);
    EXPECT_GT(coarse_error, 0.0);
    EXPECT_GT(fine_error, 0.0);
    EXPECT_GT(coarse_error / fine_error, 3.5);
    EXPECT_LT(coarse_error / fine_error, 4.5);
}

TEST(HeatEquation, PreservesSteadyLinearConduction) {
    Config config;
    config.nx = 21;
    config.ny = 13;
    config.fluid = false;
    config.left = {0.0, ThermalBoundaryType::temperature, 303.15};
    config.right = {0.0, ThermalBoundaryType::temperature, 293.15};
    Solver solver(config);
    for (int y = 1; y <= config.ny; ++y) {
        for (int x = 1; x <= config.nx; ++x) {
            solver.temperature(x, y) = 303.15 - 10.0 * (x - 0.5) / config.nx;
        }
    }
    for (int step = 0; step < 100; ++step) {
        solver.step();
    }
    for (int y = 1; y <= config.ny; ++y) {
        for (int x = 1; x <= config.nx; ++x) {
            EXPECT_NEAR(solver.temperature(x, y), 303.15 - 10.0 * (x - 0.5) / config.nx, 1e-11);
        }
    }
}
