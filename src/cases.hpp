#pragma once

#include "solver.hpp"

inline Config air_cavity() {
    Config config;
    config.top.tangential_velocity = 1.0;
    config.left.thermal_type = ThermalBoundaryType::temperature;
    config.left.temperature = 303.15;
    config.right.thermal_type = ThermalBoundaryType::temperature;
    config.right.temperature = 293.15;
    return config;
}

inline Config water_cavity() {
    Config config = air_cavity();
    config.density = 998.2;
    config.viscosity = 1.004e-6;
    config.diffusivity = 1.43e-7;
    config.expansion = 2.07e-4;
    return config;
}

inline Config poiseuille() {
    Config config;
    config.nx = 24;
    config.ny = 48;
    config.length = 0.12;
    config.height = 0.02;
    config.dt = 0.001;
    config.steps = 30000;
    config.save_every = 100;
    config.thermal = false;
    config.density = 1000.0;
    config.viscosity = 1e-5;
    config.left.flow_type = FlowBoundaryType::pressure;
    config.left.pressure = 1.2;
    config.right.flow_type = FlowBoundaryType::pressure;
    config.right.pressure = 0.0;
    return config;
}

inline Config couette() {
    Config config;
    config.nx = 24;
    config.ny = 48;
    config.length = 0.12;
    config.height = 0.02;
    config.dt = 0.001;
    config.steps = 30000;
    config.save_every = 100;
    config.thermal = false;
    config.density = 1000.0;
    config.viscosity = 1e-5;
    config.left.flow_type = FlowBoundaryType::pressure;
    config.left.pressure = 0.0;
    config.right.flow_type = FlowBoundaryType::pressure;
    config.right.pressure = 0.0;
    config.top.tangential_velocity = 0.05;
    return config;
}

inline Config lid_cavity() {
    Config config;
    config.nx = 32;
    config.ny = 32;
    config.length = 0.01;
    config.height = 0.01;
    config.dt = 0.005;
    config.steps = 8000;
    config.save_every = 40;
    config.fluid = true;
    config.thermal = false;
    config.buoyancy = false;
    config.density = 998.2;
    config.viscosity = 1.004e-6;
    config.top.tangential_velocity = 0.01004;
    return config;
}

inline Config buoyancy() {
    Config config;
    config.nx = 32;
    config.ny = 16;
    config.length = 0.04;
    config.height = 0.02;
    config.dt = 0.002;
    config.steps = 4000;
    config.save_every = 50;
    config.fluid = true;
    config.thermal = true;
    config.buoyancy = true;
    config.density = 1.204;
    config.viscosity = 1.516e-5;
    config.diffusivity = 2.14e-5;
    config.expansion = 3.43e-3;
    config.initial_temperature = 293.15;
    config.temperature_perturbation = 0.02;
    config.bottom.thermal_type = ThermalBoundaryType::temperature;
    config.bottom.temperature = 313.15;
    config.top.thermal_type = ThermalBoundaryType::temperature;
    config.top.temperature = 293.15;
    return config;
}

inline Config heated_wall() {
    Config config;
    config.nx = 48;
    config.ny = 12;
    config.length = 0.002;
    config.height = 0.001;
    config.dt = 0.001;
    config.steps = 12000;
    config.save_every = 50;
    config.fluid = false;
    config.thermal = true;
    config.buoyancy = false;
    config.density = 998.2;
    config.viscosity = 1.004e-6;
    config.diffusivity = 1.43e-7;
    config.initial_temperature = 293.15;
    config.left.thermal_type = ThermalBoundaryType::temperature;
    config.left.temperature = 313.15;
    config.right.thermal_type = ThermalBoundaryType::temperature;
    config.right.temperature = 293.15;
    return config;
}

inline Config stokes() {
    Config config;
    config.nx = 24;
    config.ny = 64;
    config.length = 0.016;
    config.height = 0.004;
    config.dt = 0.0005;
    config.steps = 16000;
    config.save_every = 40;
    config.fluid = true;
    config.thermal = false;
    config.buoyancy = false;
    config.density = 998.2;
    config.viscosity = 1.004e-6;
    config.left.flow_type = FlowBoundaryType::pressure;
    config.left.pressure = 0.0;
    config.right.flow_type = FlowBoundaryType::pressure;
    config.right.pressure = 0.0;
    config.bottom.oscillation_amplitude = 0.01;
    config.bottom.oscillation_period = 2.0;
    return config;
}

inline Config thermal_bubble() {
    Config config;
    config.nx = 48;
    config.ny = 72;
    config.length = 0.02;
    config.height = 0.03;
    config.dt = 0.002;
    config.steps = 2000;
    config.save_every = 20;
    config.pressure_relaxation = 1.7;
    config.fluid = true;
    config.thermal = true;
    config.buoyancy = true;
    config.density = 998.2;
    config.viscosity = 1.004e-6;
    config.diffusivity = 1.43e-7;
    config.expansion = 2.07e-4;
    config.initial_temperature = 293.15;
    config.reference_temperature = 293.15;
    config.bubble_x = 0.01;
    config.bubble_y = 0.006;
    config.bubble_radius = 0.003;
    config.bubble_temperature_excess = 15.0;
    return config;
}
