#include "solver.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <gtest/gtest.h>

TEST(Inlet, ParabolicProfilePreservesVolumeFlux) {
    Config config;
    config.nx = 24;
    config.ny = 16;
    config.length = 0.06;
    config.height = 0.04;
    config.thermal = false;
    config.left.flow_type = FlowBoundaryType::inlet;
    config.left.inlet_velocity = 0.01;
    config.left.inlet_center = 0.025;
    config.left.inlet_half_width = 0.01;
    config.right.flow_type = FlowBoundaryType::pressure;
    Solver solver(config);
    EXPECT_LE(solver.maximum_divergence(), config.divergence_tolerance * 1.01);
    for (int step = 0; step < 20; ++step) {
        solver.step();
    }
    double inlet_flux = 0.0;
    double outlet_flux = 0.0;
    double inlet_moment = 0.0;
    for (int y = 1; y <= config.ny; ++y) {
        double height = (y - 0.5) * solver.dy;
        double position = (height - 0.025) / 0.01;
        double expected = 0.0;
        if (std::abs(position) < 1.0) {
            expected = 0.01 * (1 - position * position);
        }
        EXPECT_DOUBLE_EQ(solver.u(0, y), expected);
        inlet_flux += solver.u(0, y) * solver.dy;
        inlet_moment += height * solver.u(0, y) * solver.dy;
        outlet_flux += solver.u(config.nx, y) * solver.dy;
    }
    EXPECT_NEAR(inlet_moment / inlet_flux, config.left.inlet_center, 1e-12);
    EXPECT_NEAR(inlet_flux, outlet_flux, config.length * config.height * config.divergence_tolerance);
    EXPECT_LE(solver.maximum_divergence(), config.divergence_tolerance * 1.01);
    for (int x = 1; x <= config.nx; ++x) {
        EXPECT_DOUBLE_EQ(solver.v(x, 0), 0.0);
        EXPECT_DOUBLE_EQ(solver.v(x, config.ny), 0.0);
    }
}

TEST(Inlet, UniformVelocityIsPrescribedWhenWidthIsZero) {
    Config config;
    config.thermal = false;
    config.left.flow_type = FlowBoundaryType::inlet;
    config.left.inlet_velocity = 0.01;
    config.right.flow_type = FlowBoundaryType::pressure;
    Solver solver(config);
    solver.step();
    for (int y = 1; y <= config.ny; ++y) {
        EXPECT_DOUBLE_EQ(solver.u(0, y), 0.01);
    }
    EXPECT_LE(solver.maximum_divergence(), config.divergence_tolerance * 1.01);
}

TEST(Inlet, RejectsInvalidOrUnsupportedBoundaries) {
    Config config;
    config.thermal = false;
    config.left.flow_type = FlowBoundaryType::inlet;
    EXPECT_THROW(Solver solver(config), std::invalid_argument);
    config.right.flow_type = FlowBoundaryType::pressure;
    config.left.inlet_half_width = -0.01;
    EXPECT_THROW(Solver solver(config), std::invalid_argument);
    config.left.inlet_half_width = 0.02;
    config.left.inlet_center = 0.0;
    EXPECT_THROW(Solver solver(config), std::invalid_argument);
    config.left.inlet_half_width = 0.0;
    config.left.inlet_velocity = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(Solver solver(config), std::invalid_argument);
    config.left.inlet_velocity = 0.01;
    config.thermal = true;
    EXPECT_THROW(Solver solver(config), std::invalid_argument);
}

TEST(Inlet, OscillatingDirectionPreservesFluxAndIncompressibilityAtEveryStep) {
    Config config;
    config.nx = 24;
    config.ny = 16;
    config.length = 0.06;
    config.height = 0.04;
    config.dt = 0.002;
    config.thermal = false;
    config.left.flow_type = FlowBoundaryType::inlet;
    config.left.inlet_velocity = 0.015;
    config.left.inlet_center = 0.025;
    config.left.inlet_half_width = 0.01;
    config.left.inlet_angle_amplitude_degrees = 10.0;
    config.left.inlet_oscillation_period = 0.08;
    config.right.flow_type = FlowBoundaryType::pressure;
    Solver solver(config);
    double pi = std::acos(-1.0);
    for (int step = 1; step <= 40; ++step) {
        solver.step();
        double angle = (10 * pi / 180) * std::sin(2 * pi * solver.time / 0.08);
        double inlet_flux = 0.0;
        double outlet_flux = 0.0;
        for (int y = 1; y <= config.ny; ++y) {
            double speed = config.left.velocity_inlet((y - 0.5) * solver.dy);
            EXPECT_NEAR(solver.u(0, y), speed * std::cos(angle), 1e-14);
            inlet_flux += solver.u(0, y) * solver.dy;
            outlet_flux += solver.u(config.nx, y) * solver.dy;
        }
        for (int y = 1; y < config.ny; ++y) {
            double speed = config.left.velocity_inlet(y * solver.dy);
            double boundary_v = 0.5 * (solver.v(0, y) + solver.v(1, y));
            EXPECT_NEAR(boundary_v, speed * std::sin(angle), 1e-14);
        }
        EXPECT_LE(solver.maximum_divergence(), config.divergence_tolerance * 1.01);
        EXPECT_NEAR(inlet_flux, outlet_flux,
                    config.length * config.height * config.divergence_tolerance);
    }
}
