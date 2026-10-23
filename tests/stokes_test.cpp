#include "cases.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <gtest/gtest.h>

double analytical_stokes(double y, double time, const Config &config) {
    if (time == 0.0) {
        return 0.0;
    }
    double pi = std::acos(-1.0);
    double amplitude = config.bottom.oscillation_amplitude;
    double omega = 2 * pi / config.bottom.oscillation_period;
    double velocity = amplitude * (1 - y / config.height) * std::sin(omega * time);
    for (int n = 1; n <= 301; ++n) {
        double wave_number = n * pi / config.height;
        double decay = config.viscosity * wave_number * wave_number;
        double coefficient = -2 * amplitude * omega / (n * pi) *
                             (decay * std::cos(omega * time) + omega * std::sin(omega * time) -
                              decay * std::exp(-decay * time)) /
                             (decay * decay + omega * omega);
        velocity += coefficient * std::sin(wave_number * y);
    }
    return velocity;
}

TEST(Stokes, WaterMatchesAnalyticalStartupAndReversal) {
    Config config = stokes();
    config.nx = 8;
    Solver solver(config);
    for (double velocity : solver.u.values) {
        EXPECT_DOUBLE_EQ(velocity, 0.0);
    }

    int current_step = 0;
    for (double time : {0.1, 0.5, 1.0, 1.5, 2.0, 3.5, 5.0, 7.5, 8.0}) {
        SCOPED_TRACE(::testing::Message() << "time = " << time << " s");
        int target_step = static_cast<int>(std::lround(time / config.dt));
        while (current_step < target_step) {
            solver.step();
            ++current_step;
        }
        double maximum_error = 0;
        for (int y = 1; y <= config.ny; ++y) {
            double exact = analytical_stokes((y - 0.5) * solver.dy, time, config);
            for (int x = 0; x <= config.nx; ++x) {
                maximum_error = std::max(maximum_error, std::abs(solver.u(x, y) - exact));
            }
        }
        EXPECT_LT(maximum_error / config.bottom.oscillation_amplitude, 0.005);
        EXPECT_NEAR(0.5 * (solver.u(0, 0) + solver.u(0, 1)), config.bottom.velocity(time), 1e-11);
        EXPECT_LE(solver.maximum_divergence(), config.divergence_tolerance);
        for (double velocity : solver.v.values) {
            EXPECT_NEAR(velocity, 0.0, 1e-12);
        }
        double left_flow = 0.0;
        double right_flow = 0.0;
        for (int y = 1; y <= config.ny; ++y) {
            EXPECT_NEAR(solver.u(0, y), solver.u(1, y), 1e-12);
            EXPECT_NEAR(solver.u(config.nx, y), solver.u(config.nx - 1, y), 1e-12);
            left_flow += solver.u(0, y) * solver.dy;
            right_flow += solver.u(config.nx, y) * solver.dy;
            for (int x = 1; x <= config.nx; ++x) {
                EXPECT_NEAR(solver.pressure(x, y), 0.0, 1e-10);
            }
        }
        EXPECT_NEAR(left_flow, right_flow, 1e-12);
        if (time == 0.5) {
            EXPECT_GT(solver.u(0, 1), 0.0);
        }
        if (time == 1.5) {
            EXPECT_LT(solver.u(0, 1), 0.0);
        }
    }
}

TEST(Stokes, ApproachesTimePeriodicFiniteChannelSolution) {
    Config config = stokes();
    config.nx = 8;
    Solver solver(config);
    for (int step = 0; step < config.steps; ++step) {
        solver.step();
    }

    double omega = 2 * std::acos(-1.0) / config.bottom.oscillation_period;
    std::complex<double> wave_number = std::sqrt(std::complex<double>(0, omega / config.viscosity));
    std::complex<double> phase = std::exp(std::complex<double>(0, omega * solver.time));
    double maximum_error = 0;
    for (int y = 1; y <= config.ny; ++y) {
        double position = (y - 0.5) * solver.dy;
        double exact = std::imag(config.bottom.oscillation_amplitude * phase *
                                 std::sinh(wave_number * (config.height - position)) /
                                 std::sinh(wave_number * config.height));
        maximum_error = std::max(maximum_error, std::abs(solver.u(0, y) - exact));
    }
    EXPECT_LT(maximum_error / config.bottom.oscillation_amplitude, 0.003);
}
