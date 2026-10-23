#include "cases.hpp"

#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>

double analytical_couette(double y, double time, const Config &config) {
    if (time == 0.0) {
        return 0.0;
    }
    double pi = std::acos(-1.0);
    double wall_velocity = config.top.tangential_velocity;
    double velocity = wall_velocity * y / config.height;
    for (int n = 1; n <= 201; ++n) {
        double wave_number = n * pi / config.height;
        double sign = -1.0;
        if (n % 2 == 0) {
            sign = 1.0;
        }
        velocity += 2 * wall_velocity * sign / (n * pi) * std::sin(wave_number * y) *
                    std::exp(-config.viscosity * wave_number * wave_number * time);
    }
    return velocity;
}

TEST(Couette, MatchesAnalyticalStartupAtPhysicalTimes) {
    Config config = couette();
    config.nx = 8;
    config.ny = 64;
    config.dt = 0.0005;
    Solver solver(config);
    for (int y = 1; y <= config.ny; ++y) {
        for (int x = 0; x <= config.nx; ++x) {
            EXPECT_DOUBLE_EQ(solver.u(x, y), 0.0);
        }
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
        for (int y = 1; y <= config.ny; ++y) {
            double position = (y - 0.5) * solver.dy;
            double exact = analytical_couette(position, time, config);
            for (int x = 0; x <= config.nx; ++x) {
                maximum_error = std::max(maximum_error, std::abs(solver.u(x, y) - exact));
            }
        }
        EXPECT_LT(maximum_error / config.top.tangential_velocity, 0.005);
        EXPECT_LE(solver.maximum_divergence(), config.divergence_tolerance);
        for (double velocity : solver.v.values) {
            EXPECT_NEAR(velocity, 0.0, 1e-12);
        }
        if (time == 2.0) {
            EXPECT_LT(solver.u(0, config.ny / 2), 0.3 * config.top.tangential_velocity);
        }
    }
}

TEST(Couette, ApproachesSteadyLinearProfileAfterThirtySeconds) {
    Config config = couette();
    Solver solver(config);
    for (int step = 0; step < config.steps; ++step) {
        solver.step();
    }

    double maximum_error = 0;
    for (int y = 1; y <= config.ny; ++y) {
        double position = (y - 0.5) * solver.dy;
        double steady = config.top.tangential_velocity * position / config.height;
        for (int x = 0; x <= config.nx; ++x) {
            maximum_error = std::max(maximum_error, std::abs(solver.u(x, y) - steady));
        }
    }
    EXPECT_NEAR(0.5 * (solver.u(0, 0) + solver.u(0, 1)), 0.0, 1e-12);
    EXPECT_NEAR(0.5 * (solver.u(0, config.ny) + solver.u(0, config.ny + 1)),
                config.top.tangential_velocity, 1e-12);
    EXPECT_LT(maximum_error / config.top.tangential_velocity, 0.002);
    EXPECT_LE(solver.maximum_divergence(), config.divergence_tolerance);
}

TEST(Couette, PeriodicProjectionRemovesNonuniformDivergence) {
    Config config = couette();
    config.nx = 16;
    config.ny = 12;
    config.top.tangential_velocity = 0.0;
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
