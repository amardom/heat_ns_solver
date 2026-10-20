#include "cases.hpp"

#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>

double analytical_poiseuille(double y, double time, const Config &config) {
    if (time == 0.0) {
        return 0.0;
    }
    double pi = std::acos(-1.0);
    double velocity = config.acceleration_x * y * (config.height - y) / (2 * config.viscosity);
    for (int n = 1; n <= 201; n += 2) {
        double wave_number = n * pi / config.height;
        velocity -= 4 * config.acceleration_x * config.height * config.height /
                    (config.viscosity * pi * pi * pi * n * n * n) * std::sin(wave_number * y) *
                    std::exp(-config.viscosity * wave_number * wave_number * time);
    }
    return velocity;
}

TEST(Poiseuille, MatchesAnalyticalStartupAtPhysicalTimes) {
    Config config = poiseuille();
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

        if (time == 2.0) {
            double position = (config.ny / 2 - 0.5) * solver.dy;
            double steady = config.acceleration_x * position * (config.height - position) /
                            (2 * config.viscosity);
            EXPECT_LT(solver.u(0, config.ny / 2), 0.6 * steady);
        }
    }
}

TEST(Poiseuille, PeriodicProjectionRemovesNonuniformDivergence) {
    Config config = poiseuille();
    config.nx = 16;
    config.ny = 12;
    config.acceleration_x = 0;
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

TEST(Poiseuille, ApproachesSteadyParabolaAfterThirtySeconds) {
    Config config = poiseuille();
    Solver solver(config);
    for (int step = 0; step < config.steps; ++step) {
        solver.step();
    }

    double maximum_error = 0;
    double maximum_steady =
        config.acceleration_x * config.height * config.height / (8 * config.viscosity);
    for (int y = 1; y <= config.ny; ++y) {
        double position = (y - 0.5) * solver.dy;
        double steady =
            config.acceleration_x * position * (config.height - position) / (2 * config.viscosity);
        for (int x = 0; x <= config.nx; ++x) {
            maximum_error = std::max(maximum_error, std::abs(solver.u(x, y) - steady));
        }
    }
    EXPECT_LT(maximum_error / maximum_steady, 0.002);
    EXPECT_LE(solver.maximum_divergence(), config.divergence_tolerance);
}
