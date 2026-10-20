#include "solver.hpp"

#include <gtest/gtest.h>

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
