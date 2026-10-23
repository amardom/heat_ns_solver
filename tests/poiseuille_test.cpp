#include "cases.hpp"

#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>

double analytical_poiseuille(double y, double time, const Config &config) {
    if (time == 0.0) {
        return 0.0;
    }
    double acceleration = (config.left.pressure - config.right.pressure) /
                          (config.density * config.length);
    double pi = std::acos(-1.0);
    double velocity = acceleration * y * (config.height - y) / (2 * config.viscosity);
    for (int n = 1; n <= 201; n += 2) {
        double wave_number = n * pi / config.height;
        velocity -= 4 * acceleration * config.height * config.height /
                    (config.viscosity * pi * pi * pi * n * n * n) * std::sin(wave_number * y) *
                    std::exp(-config.viscosity * wave_number * wave_number * time);
    }
    return velocity;
}

TEST(Poiseuille, MatchesAnalyticalStartupAtPhysicalTimes) {
    Config config = poiseuille();
    double acceleration = (config.left.pressure - config.right.pressure) /
                          (config.density * config.length);
    config.nx = 8;
    config.ny = 64;
    config.dt = 0.0005;
    Solver solver(config);
    for (double velocity : solver.u.values) {
        EXPECT_DOUBLE_EQ(velocity, 0.0);
    }

    int current_step = 0;
    for (double time : {0.1, 0.5, 2.0, 5.0, 8.0, 15.0, 30.0}) {
        SCOPED_TRACE(::testing::Message() << "time = " << time << " s");
        int target_step = static_cast<int>(std::lround(time / config.dt));
        while (current_step < target_step) {
            solver.step();
            ++current_step;
        }
        double maximum_error = 0;
        double maximum_exact = 0;
        for (int y = 1; y <= config.ny; ++y) {
            double position = (y - 0.5) * solver.dy;
            double exact = analytical_poiseuille(position, time, config);
            maximum_exact = std::max(maximum_exact, exact);
            for (int x = 0; x <= config.nx; ++x) {
                maximum_error = std::max(maximum_error, std::abs(solver.u(x, y) - exact));
            }
        }
        EXPECT_LT(maximum_error / maximum_exact, 0.015);
        EXPECT_LE(solver.maximum_divergence(), config.divergence_tolerance);
        for (double velocity : solver.v.values) {
            EXPECT_NEAR(velocity, 0.0, 1e-12);
        }

        for (int y = 1; y <= config.ny; ++y) {
            for (int x = 1; x <= config.nx; ++x) {
                double fraction = (x - 0.5) * solver.dx / config.length;
                double expected_pressure = config.left.pressure +
                    fraction * (config.right.pressure - config.left.pressure);
                EXPECT_NEAR(solver.pressure(x, y), expected_pressure, 1e-10);
            }
        }
        if (time == 2.0) {
            double position = (config.ny / 2 - 0.5) * solver.dy;
            double steady = acceleration * position * (config.height - position) /
                            (2 * config.viscosity);
            EXPECT_LT(solver.u(0, config.ny / 2), 0.6 * steady);
        }
    }
}

TEST(Poiseuille, PressureBoundariesRemoveNonuniformDivergence) {
    Config config = poiseuille();
    config.nx = 16;
    config.ny = 12;
    config.left.pressure = config.right.pressure;
    Solver solver(config);
    double pi = std::acos(-1.0);
    for (int y = 1; y <= config.ny; ++y) {
        for (int x = 0; x < config.nx; ++x) {
            solver.u(x, y) =
                0.001 * std::sin(2 * pi * x / config.nx) * std::sin(pi * (y - 0.5) / config.ny);
        }
    }
    solver.apply_velocity_boundaries(solver.time);
    EXPECT_GT(solver.maximum_divergence(), 1e-3);
    solver.step();
    EXPECT_LE(solver.maximum_divergence(), config.divergence_tolerance * 1.01);
    double inlet_flow = 0.0;
    double outlet_flow = 0.0;
    for (int y = 1; y <= config.ny; ++y) {
        inlet_flow += solver.u(0, y) * solver.dy;
        outlet_flow += solver.u(config.nx, y) * solver.dy;
        EXPECT_DOUBLE_EQ(solver.v(0, y), -solver.v(1, y));
        EXPECT_DOUBLE_EQ(solver.v(config.nx + 1, y), -solver.v(config.nx, y));
    }
    EXPECT_NEAR(inlet_flow, outlet_flow,
                config.length * config.height * config.divergence_tolerance);
}

TEST(Poiseuille, ApproachesSteadyParabolaAfterThirtySeconds) {
    Config config = poiseuille();
    double acceleration = (config.left.pressure - config.right.pressure) /
                          (config.density * config.length);
    Solver solver(config);
    for (int step = 0; step < config.steps; ++step) {
        solver.step();
    }

    double maximum_error = 0;
    double maximum_steady =
        acceleration * config.height * config.height / (8 * config.viscosity);
    for (int y = 1; y <= config.ny; ++y) {
        double position = (y - 0.5) * solver.dy;
        double steady =
            acceleration * position * (config.height - position) / (2 * config.viscosity);
        for (int x = 0; x <= config.nx; ++x) {
            maximum_error = std::max(maximum_error, std::abs(solver.u(x, y) - steady));
        }
    }
    EXPECT_LT(maximum_error / maximum_steady, 0.002);
    EXPECT_LE(solver.maximum_divergence(), config.divergence_tolerance);
}

TEST(Poiseuille, PressureReferenceChangesPressureButNotVelocity) {
    Config config = poiseuille();
    config.nx = 8;
    config.ny = 12;
    Solver original(config);
    double offset = 7.0;
    config.left.pressure += offset;
    config.right.pressure += offset;
    Solver shifted(config);
    for (int step = 0; step < 20; ++step) {
        original.step();
        shifted.step();
    }
    for (int y = 1; y <= config.ny; ++y) {
        for (int x = 0; x <= config.nx; ++x) {
            EXPECT_NEAR(original.u(x, y), shifted.u(x, y), 1e-12);
        }
        for (int x = 1; x <= config.nx; ++x) {
            EXPECT_NEAR(shifted.pressure(x, y) - original.pressure(x, y), offset, 1e-10);
        }
    }
    EXPECT_LE(shifted.maximum_divergence(), config.divergence_tolerance);
}

TEST(Poiseuille, EqualPressuresKeepFluidAtRest) {
    Config config = poiseuille();
    config.nx = 8;
    config.ny = 12;
    config.left.pressure = 5.0;
    config.right.pressure = 5.0;
    Solver solver(config);
    for (int step = 0; step < 20; ++step) {
        solver.step();
    }
    for (double velocity : solver.u.values) {
        EXPECT_NEAR(velocity, 0.0, 1e-12);
    }
    for (double velocity : solver.v.values) {
        EXPECT_NEAR(velocity, 0.0, 1e-12);
    }
    for (int y = 1; y <= config.ny; ++y) {
        for (int x = 1; x <= config.nx; ++x) {
            EXPECT_NEAR(solver.pressure(x, y), 5.0, 1e-12);
        }
    }
}
