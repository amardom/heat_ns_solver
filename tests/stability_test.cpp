#include "solver.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>

TEST(Stability, CalculatesCourantAndFourierOnRectangle) {
    Config config;
    config.nx = 4;
    config.ny = 3;
    config.length = 2.0;
    config.height = 3.0;
    config.dt = 0.1;
    config.viscosity = 0.2;
    config.diffusivity = 0.1;
    Solver solver(config);

    for (double &velocity : solver.u.values) {
        velocity = 0.3;
    }
    for (double &velocity : solver.v.values) {
        velocity = -0.2;
    }

    TimeStepLimits limits = solver.time_step_limits();
    EXPECT_NEAR(limits.courant_x, 0.06, 1e-14);
    EXPECT_NEAR(limits.courant_y, 0.02, 1e-14);
    EXPECT_NEAR(limits.courant, 0.08, 1e-14);
    EXPECT_NEAR(limits.viscous_fourier, 0.1, 1e-14);
    EXPECT_NEAR(limits.thermal_fourier, 0.05, 1e-14);
    EXPECT_NEAR(limits.momentum, 0.28, 1e-14);
    EXPECT_NEAR(limits.thermal, 0.18, 1e-14);
    EXPECT_NEAR(limits.maximum_dt, config.dt / 0.28, 1e-14);
}

TEST(Stability, RejectsCombinedLimitEvenWhenSeparateLimitsPass) {
    Config config;
    config.nx = 4;
    config.ny = 4;
    config.length = 1.0;
    config.height = 1.0;
    config.dt = 0.1;
    config.viscosity = 0.1;
    config.thermal = false;
    Solver solver(config);
    for (double &velocity : solver.u.values) {
        velocity = 1.0;
    }

    TimeStepLimits limits = solver.time_step_limits();
    EXPECT_LT(limits.courant, 1.0);
    EXPECT_LT(limits.viscous_fourier, 0.5);
    EXPECT_GT(limits.combined, 1.0);
    EXPECT_THROW(solver.step(), std::runtime_error);
    EXPECT_DOUBLE_EQ(solver.time, 0.0);
    EXPECT_DOUBLE_EQ(solver.u(1, 1), 1.0);
}

TEST(Stability, RejectsExcessiveThermalDiffusionBeforeAdvancing) {
    Config config;
    config.nx = 4;
    config.ny = 4;
    config.length = 1.0;
    config.height = 1.0;
    config.dt = 0.1;
    config.diffusivity = 0.2;
    config.fluid = false;
    Solver solver(config);

    EXPECT_GT(solver.time_step_limits().thermal_fourier, 0.5);
    EXPECT_THROW(solver.step(), std::runtime_error);
    EXPECT_DOUBLE_EQ(solver.time, 0.0);
    EXPECT_DOUBLE_EQ(solver.temperature(2, 2), config.initial_temperature);
}

TEST(Stability, DisabledMomentumDoesNotRestrictThermalTimeStep) {
    Config config;
    config.nx = 4;
    config.ny = 4;
    config.length = 1.0;
    config.height = 1.0;
    config.dt = 0.1;
    config.viscosity = 1000.0;
    config.diffusivity = 0.01;
    config.fluid = false;
    Solver solver(config);

    TimeStepLimits limits = solver.time_step_limits();
    EXPECT_DOUBLE_EQ(limits.viscous_fourier, 0.0);
    EXPECT_DOUBLE_EQ(limits.momentum, 0.0);
    EXPECT_LT(limits.thermal, 1.0);
    EXPECT_NO_THROW(solver.step());
}

TEST(Stability, NoTransportHasNoTimeStepRestriction) {
    Config config;
    config.viscosity = 0.0;
    config.diffusivity = 0.0;
    Solver solver(config);

    TimeStepLimits limits = solver.time_step_limits();
    EXPECT_DOUBLE_EQ(limits.combined, 0.0);
    EXPECT_TRUE(std::isinf(limits.maximum_dt));
    EXPECT_NO_THROW(solver.step());
}

TEST(Stability, RejectsNonfiniteVelocity) {
    Config config;
    Solver solver(config);
    solver.u(1, 1) = std::numeric_limits<double>::quiet_NaN();

    EXPECT_THROW(solver.time_step_limits(), std::runtime_error);
    EXPECT_THROW(solver.step(), std::runtime_error);
    EXPECT_DOUBLE_EQ(solver.time, 0.0);
}

TEST(Stability, RejectsNonfiniteTemperature) {
    Config config;
    config.fluid = false;
    Solver solver(config);
    solver.temperature(2, 2) = std::numeric_limits<double>::quiet_NaN();

    EXPECT_THROW(solver.step(), std::runtime_error);
    EXPECT_DOUBLE_EQ(solver.time, 0.0);
}

TEST(Stability, RejectsNonfinitePressure) {
    Config config;
    config.thermal = false;
    Solver solver(config);
    solver.pressure(2, 2) = std::numeric_limits<double>::quiet_NaN();

    EXPECT_THROW(solver.step(), std::runtime_error);
    EXPECT_DOUBLE_EQ(solver.time, 0.0);
}

TEST(Stability, RejectsNonfiniteConfiguration) {
    Config config;
    config.dt = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(Solver solver(config), std::invalid_argument);

    config.dt = std::numeric_limits<double>::infinity();
    EXPECT_THROW(Solver solver(config), std::invalid_argument);
}
