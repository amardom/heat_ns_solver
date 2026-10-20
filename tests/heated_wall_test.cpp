#include "cases.hpp"

#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>

double analytical_heated_wall(double x, double time, const Config &config) {
    if (time == 0.0) {
        return config.initial_temperature;
    }
    double pi = std::acos(-1.0);
    double difference = config.left.temperature - config.right.temperature;
    double temperature = config.right.temperature + difference * (1 - x / config.length);
    for (int n = 1; n <= 201; ++n) {
        double wave_number = n * pi / config.length;
        temperature -= 2 * difference / (n * pi) * std::sin(wave_number * x) *
                       std::exp(-config.diffusivity * wave_number * wave_number * time);
    }
    return temperature;
}

TEST(HeatedWall, WaterMatchesAnalyticalTransientAtPhysicalTimes) {
    Config config = heated_wall();
    Solver solver(config);
    for (int y = 1; y <= config.ny; ++y) {
        for (int x = 1; x <= config.nx; ++x) {
            EXPECT_DOUBLE_EQ(solver.temperature(x, y), config.initial_temperature);
        }
    }

    int current_step = 0;
    for (double time : {0.1, 0.5, 1.0, 2.0, 5.0, 12.0}) {
        SCOPED_TRACE(::testing::Message() << "time = " << time << " s");
        int target_step = static_cast<int>(std::lround(time / config.dt));
        while (current_step < target_step) {
            solver.step();
            ++current_step;
        }
        double maximum_error = 0;
        for (int y = 1; y <= config.ny; ++y) {
            for (int x = 1; x <= config.nx; ++x) {
                double position = (x - 0.5) * solver.dx;
                double exact = analytical_heated_wall(position, time, config);
                maximum_error = std::max(maximum_error, std::abs(solver.temperature(x, y) - exact));
                EXPECT_GE(solver.temperature(x, y), config.right.temperature);
                EXPECT_LE(solver.temperature(x, y), config.left.temperature);
            }
        }
        EXPECT_LT(maximum_error / (config.left.temperature - config.right.temperature), 0.005);
        for (double velocity : solver.u.values) {
            EXPECT_DOUBLE_EQ(velocity, 0.0);
        }
        for (double velocity : solver.v.values) {
            EXPECT_DOUBLE_EQ(velocity, 0.0);
        }
    }
}
