#pragma once

#include <vector>

struct Field {
    int nx;
    int ny;
    std::vector<double> values;

    Field(int nx, int ny, double value = 0.0);
    double &operator()(int x, int y);
    double operator()(int x, int y) const;
};

enum class FlowBoundaryType { wall, pressure };

enum class ThermalBoundaryType { temperature, insulated };

struct Boundary {
    double tangential_velocity = 0.0;
    ThermalBoundaryType thermal_type = ThermalBoundaryType::insulated;
    double temperature = 293.15;
    double oscillation_amplitude = 0.0;
    double oscillation_period = 0.0;
    FlowBoundaryType flow_type = FlowBoundaryType::wall;
    double pressure = 0.0;

    double velocity(double time) const;
};

struct Config {
    int nx = 32;
    int ny = 24;
    double length = 0.1;
    double height = 0.075;
    double dt = 0.0005;
    int steps = 1000;
    int save_every = 50;
    bool fluid = true;
    bool thermal = true;
    bool buoyancy = false;
    double density = 1.204;
    double viscosity = 1.516e-5;
    double diffusivity = 2.14e-5;
    double expansion = 3.43e-3;
    double gravity = 9.81;
    double reference_temperature = 293.15;
    double initial_temperature = 293.15;
    double temperature_perturbation = 0.0;
    double bubble_x = 0.0;
    double bubble_y = 0.0;
    double bubble_radius = 0.0;
    double bubble_temperature_excess = 0.0;
    int pressure_iterations = 20000;
    double pressure_relaxation = 1.0;
    double divergence_tolerance = 1e-7;
    Boundary left;
    Boundary right;
    Boundary bottom;
    Boundary top;
};

struct TimeStepLimits {
    double courant_x = 0.0;
    double courant_y = 0.0;
    double courant = 0.0;
    double viscous_fourier = 0.0;
    double thermal_fourier = 0.0;
    double momentum = 0.0;
    double thermal = 0.0;
    double combined = 0.0;
    double maximum_dt = 0.0;
};

struct Solver {
    Config config;
    double dx;
    double dy;
    double time = 0.0;
    Field u;
    Field v;
    Field pressure;
    Field temperature;

    explicit Solver(const Config &config);
    void apply_velocity_boundaries();
    void apply_temperature_boundaries();
    void advance_temperature();
    void advance_fluid();
    void predict_velocity(Field &predicted_u, Field &predicted_v);
    void solve_pressure(const Field &predicted_u, const Field &predicted_v);
    void correct_velocity(const Field &predicted_u, const Field &predicted_v);
    TimeStepLimits time_step_limits() const;
    void step();
    double maximum_divergence() const;
};
