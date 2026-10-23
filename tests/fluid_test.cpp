#include "solver.hpp"

#include <gtest/gtest.h>
#include <cmath>

TEST(Fluid, UniformBuoyancyMatchesHydrostaticPressure) {
    Config config;
    config.nx = 12;
    config.ny = 8;
    config.buoyancy = true;
    config.initial_temperature = config.reference_temperature + 10;
    config.divergence_tolerance = 1e-10;
    Solver solver(config);
    solver.step();
    double force = config.gravity * config.expansion * 10;
    for (int y = 1; y <= config.ny; ++y) {
        for (int x = 1; x <= config.nx; ++x) {
            EXPECT_NEAR(solver.pressure(x, y), config.density * force * (y - 1) * solver.dy, 1e-7);
        }
    }
    EXPECT_LE(solver.maximum_divergence(), config.divergence_tolerance * 1.01);
}

TEST(Fluid, MovingLidPreservesIncompressibilityOnRectangle) {
    Config config;
    config.nx = 18;
    config.ny = 12;
    config.thermal = false;
    config.top.tangential_velocity = 1.0;
    Solver solver(config);
    for (int step = 0; step < 20; ++step) {
        solver.step();
    }
    EXPECT_GT(solver.u(config.nx / 2, config.ny), 0.0);
    EXPECT_LE(solver.maximum_divergence(), config.divergence_tolerance * 1.01);
}

TEST(Fluid, PeriodicProjectionRemovesNonuniformDivergence) {
    Config config;
    config.length = 0.12;
    config.height = 0.02;
    config.dt = 0.001;
    config.viscosity = 1e-5;
    config.thermal = false;
    config.periodic_x = true;
    config.nx = 16;
    config.ny = 12;
    Solver solver(config);
    double pi = std::acos(-1.0);
    for (int y = 1; y <= config.ny; ++y) {
        for (int x = 0; x < config.nx; ++x) {
            solver.u(x, y) =
                0.001 * std::sin(2 * pi * x / config.nx) * std::sin(pi * (y - 0.5) / config.ny);
        }
    }
    solver.apply_velocity_boundaries();
    EXPECT_GT(solver.maximum_divergence(), 1e-3);
    solver.step();
    EXPECT_LE(solver.maximum_divergence(), config.divergence_tolerance * 1.01);
    for (int y = 1; y <= config.ny; ++y) {
        EXPECT_DOUBLE_EQ(solver.u(0, y), solver.u(config.nx, y));
    }
}
