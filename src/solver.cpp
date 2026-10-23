#include "solver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

Field::Field(int nx, int ny, double value) : nx(nx), ny(ny), values(nx * ny, value) {
}

double &Field::operator()(int x, int y) {
    return values[y * nx + x];
}
double Field::operator()(int x, int y) const {
    return values[y * nx + x];
}

double Boundary::velocity_wall(double time) const {
    if (oscillation_period_stokes <= 0.0) {
        return tangential_velocity;
    }
    double pi = std::acos(-1.0);
    return tangential_velocity +
           oscillation_amplitude_stokes * std::sin(2 * pi * time / oscillation_period_stokes);
}

double Boundary::velocity_inlet(double height) const {
    if (inlet_half_width == 0.0) {
        return inlet_velocity;
    }
    double position = (height - inlet_center) / inlet_half_width;
    if (std::abs(position) >= 1.0) {
        return 0.0;
    }
    return inlet_velocity * (1.0 - position * position);
}

double Boundary::angle_inlet(double time) const {
    if (inlet_oscillation_period == 0.0) {
        return 0.0;
    }
    double pi = std::acos(-1.0);
    double amplitude = inlet_angle_amplitude_degrees * pi / 180.0;
    return amplitude * std::sin(2 * pi * time / inlet_oscillation_period);
}

Config validate_config(const Config &parameters) {
    for (double value :
         {parameters.length, parameters.height, parameters.dt, parameters.density,
          parameters.viscosity, parameters.diffusivity, parameters.expansion, parameters.gravity,
          parameters.reference_temperature, parameters.initial_temperature,
          parameters.temperature_perturbation, parameters.bubble_x, parameters.bubble_y,
          parameters.bubble_radius, parameters.bubble_temperature_excess,
          parameters.pressure_relaxation, parameters.divergence_tolerance}) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Simulation parameters must be finite");
        }
    }
    for (const Boundary &boundary :
         {parameters.left, parameters.right, parameters.bottom, parameters.top}) {
        for (double value : {boundary.pressure, boundary.tangential_velocity, boundary.temperature,
                             boundary.inlet_velocity, boundary.inlet_center, boundary.inlet_half_width,
                             boundary.inlet_angle_amplitude_degrees, boundary.inlet_oscillation_period,
                             boundary.oscillation_amplitude_stokes, boundary.oscillation_period_stokes}) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument("Boundary parameters must be finite");
            }
        }
    }
    if (parameters.nx < 2 || parameters.ny < 2 || parameters.length <= 0 ||
        parameters.height <= 0 || parameters.dt <= 0 || parameters.density <= 0 ||
        parameters.viscosity < 0 || parameters.diffusivity < 0 || parameters.steps < 0 ||
        parameters.save_every < 1 || parameters.pressure_iterations < 1 ||
        parameters.divergence_tolerance <= 0 || parameters.pressure_relaxation <= 0 ||
        parameters.pressure_relaxation >= 2) {
        throw std::invalid_argument("Invalid simulation parameters");
    }
    if (parameters.top.flow_type != FlowBoundaryType::wall ||
        parameters.bottom.flow_type != FlowBoundaryType::wall) {
        throw std::invalid_argument("Pressure boundaries are supported on left and right only");
    }
    if (parameters.right.flow_type == FlowBoundaryType::inlet ||
        (parameters.left.flow_type == FlowBoundaryType::inlet &&
         parameters.right.flow_type != FlowBoundaryType::pressure)) {
        throw std::invalid_argument("Velocity inlet requires left inlet and right pressure outlet");
    }
    if (parameters.left.flow_type == FlowBoundaryType::inlet) {
        if (!parameters.fluid || parameters.thermal) {
            throw std::invalid_argument("Velocity inlet currently supports fluid-only simulations");
        }
        const Boundary &inlet = parameters.left;
        if (inlet.inlet_angle_amplitude_degrees < 0.0 || inlet.inlet_angle_amplitude_degrees >= 90.0 ||
            inlet.inlet_oscillation_period < 0.0 ||
            (inlet.inlet_angle_amplitude_degrees > 0.0 && inlet.inlet_oscillation_period == 0.0)) {
            throw std::invalid_argument("Inlet oscillation requires amplitude below 90 degrees and positive period");
        }
        if (inlet.inlet_half_width < 0.0 ||
            (inlet.inlet_half_width > 0.0 &&
             (inlet.inlet_center - inlet.inlet_half_width < 0.0 ||
              inlet.inlet_center + inlet.inlet_half_width > parameters.height))) {
            throw std::invalid_argument("Inlet profile must fit inside the channel height");
        }
    }
    return parameters;
}

Solver::Solver(const Config &parameters)
    : config(validate_config(parameters)), dx(parameters.length / parameters.nx),
      dy(parameters.height / parameters.ny), u(parameters.nx + 1, parameters.ny + 2),
      v(parameters.nx + 2, parameters.ny + 1), pressure(parameters.nx + 2, parameters.ny + 2),
      temperature(parameters.nx + 2, parameters.ny + 2, parameters.initial_temperature) {
    if (parameters.left.flow_type == FlowBoundaryType::pressure &&
        parameters.right.flow_type == FlowBoundaryType::pressure) {
        for (int y = 1; y <= parameters.ny; ++y) {
            for (int x = 1; x <= parameters.nx; ++x) {
                double fraction = (x - 0.5) * dx / parameters.length;
                pressure(x, y) = parameters.left.pressure +
                                fraction * (parameters.right.pressure - parameters.left.pressure);
            }
        }
    }
    if (parameters.temperature_perturbation != 0.0) {
        double pi = std::acos(-1.0);
        for (int y = 1; y <= parameters.ny; ++y) {
            for (int x = 1; x <= parameters.nx; ++x) {
                temperature(x, y) += parameters.temperature_perturbation *
                                     std::sin(2 * pi * (x - 0.5) * dx / parameters.length) *
                                     std::sin(pi * (y - 0.5) * dy / parameters.height);
            }
        }
    }
    if (parameters.bubble_radius > 0.0) {
        for (int y = 1; y <= parameters.ny; ++y) {
            for (int x = 1; x <= parameters.nx; ++x) {
                double distance = std::hypot((x - 0.5) * dx - parameters.bubble_x,
                                             (y - 0.5) * dy - parameters.bubble_y);
                double shape = 0.5 * (1 - std::tanh((distance - parameters.bubble_radius) /
                                                    (0.15 * parameters.bubble_radius)));
                temperature(x, y) += parameters.bubble_temperature_excess * shape;
            }
        }
    }
    apply_velocity_boundaries(u, v, time);
    apply_temperature_boundaries();
    if (config.left.flow_type == FlowBoundaryType::inlet) {
        Field initial_u = u;
        Field initial_v = v;
        solve_pressure(initial_u, initial_v);
        correct_velocity(initial_u, initial_v);
        apply_velocity_boundaries(u, v, time);
        std::fill(pressure.values.begin(), pressure.values.end(), 0.0);
    }
}

void Solver::apply_velocity_boundaries(Field &boundary_u, Field &boundary_v, double boundary_time) {
    const Config &parameters = config;
    double inlet_angle = parameters.left.angle_inlet(boundary_time);
    for (int y = 1; y <= parameters.ny; ++y) {
        if (parameters.left.flow_type == FlowBoundaryType::wall) {
            boundary_u(0, y) = 0;
        }
        if (parameters.left.flow_type == FlowBoundaryType::inlet) {
            double height = (y - 0.5) * dy;
            boundary_u(0, y) = parameters.left.velocity_inlet(height) * std::cos(inlet_angle);
        }
        if (parameters.right.flow_type == FlowBoundaryType::wall) {
            boundary_u(parameters.nx, y) = 0;
        }
    }
    for (int x = 0; x <= parameters.nx; ++x) {
        boundary_u(x, 0) = 2 * parameters.bottom.velocity_wall(boundary_time) - boundary_u(x, 1);
        boundary_u(x, parameters.ny + 1) = 2 * parameters.top.velocity_wall(boundary_time) - boundary_u(x, parameters.ny);
    }
    for (int x = 1; x <= parameters.nx; ++x) {
        boundary_v(x, 0) = 0;
        boundary_v(x, parameters.ny) = 0;
    }
    for (int y = 0; y <= parameters.ny; ++y) {
        boundary_v(0, y) = 2 * parameters.left.velocity_wall(boundary_time) - boundary_v(1, y);
        if (parameters.left.flow_type == FlowBoundaryType::inlet) {
            double height = y * dy;
            double inlet_v = parameters.left.velocity_inlet(height) * std::sin(inlet_angle);
            boundary_v(0, y) = 2 * inlet_v - boundary_v(1, y);
        }
        boundary_v(parameters.nx + 1, y) = 2 * parameters.right.velocity_wall(boundary_time) - boundary_v(parameters.nx, y);
    }
}

double ghost_temperature(double inside, const Boundary &wall) {
    if (wall.thermal_type == ThermalBoundaryType::insulated) {
        return inside;
    }
    return 2 * wall.temperature - inside;
}

void Solver::apply_temperature_boundaries() {
    const Config &parameters = config;
    for (int y = 1; y <= parameters.ny; ++y) {
        temperature(0, y) = ghost_temperature(temperature(1, y), parameters.left);
        temperature(parameters.nx + 1, y) =
            ghost_temperature(temperature(parameters.nx, y), parameters.right);
    }
    for (int x = 1; x <= parameters.nx; ++x) {
        temperature(x, 0) = ghost_temperature(temperature(x, 1), parameters.bottom);
        temperature(x, parameters.ny + 1) =
            ghost_temperature(temperature(x, parameters.ny), parameters.top);
    }
}

void Solver::update_temperature() {
    Field next_temperature = temperature;
    for (int y = 1; y <= config.ny; ++y) {
        for (int x = 1; x <= config.nx; ++x) {
            double horizontal_velocity = 0.5 * (u(x - 1, y) + u(x, y));
            double vertical_velocity = 0.5 * (v(x, y - 1) + v(x, y));
            double dT_dx;
            if (horizontal_velocity >= 0) {
                dT_dx = (temperature(x, y) - temperature(x - 1, y)) / dx;
            } else {
                dT_dx = (temperature(x + 1, y) - temperature(x, y)) / dx;
            }

            double dT_dy;
            if (vertical_velocity >= 0) {
                dT_dy = (temperature(x, y) - temperature(x, y - 1)) / dy;
            } else {
                dT_dy = (temperature(x, y + 1) - temperature(x, y)) / dy;
            }
            double d2T_dx2 =
                (temperature(x + 1, y) - 2 * temperature(x, y) + temperature(x - 1, y)) / (dx * dx);
            double d2T_dy2 =
                (temperature(x, y + 1) - 2 * temperature(x, y) + temperature(x, y - 1)) / (dy * dy);
            double laplacian = d2T_dx2 + d2T_dy2;
            next_temperature(x, y) +=
                config.dt * (config.diffusivity * laplacian - horizontal_velocity * dT_dx -
                             vertical_velocity * dT_dy);
            if (!std::isfinite(next_temperature(x, y))) {
                throw std::runtime_error("Temperature is not finite");
            }
        }
    }
    temperature = next_temperature;
}

void Solver::predict_velocity(Field &predicted_u, Field &predicted_v) {
    const Config &parameters = config;
    for (int y = 1; y <= parameters.ny; ++y) {
        for (int x = 1; x < parameters.nx; ++x) {
            double vertical_velocity =
                0.25 * (v(x, y - 1) + v(x + 1, y - 1) + v(x, y) + v(x + 1, y));
            double horizontal_velocity = u(x, y);
            double du_dx;
            if (horizontal_velocity >= 0) {
                du_dx = (horizontal_velocity - u(x - 1, y)) / dx;
            } else {
                du_dx = (u(x + 1, y) - horizontal_velocity) / dx;
            }
            double du_dy;
            if (vertical_velocity >= 0) {
                du_dy = (horizontal_velocity - u(x, y - 1)) / dy;
            } else {
                du_dy = (u(x, y + 1) - horizontal_velocity) / dy;
            }
            double d2u_dx2 = (u(x + 1, y) - 2 * horizontal_velocity + u(x - 1, y)) / (dx * dx);
            double d2u_dy2 = (u(x, y + 1) - 2 * horizontal_velocity + u(x, y - 1)) / (dy * dy);
            double laplacian = d2u_dx2 + d2u_dy2;
            predicted_u(x, y) +=
                parameters.dt * (-horizontal_velocity * du_dx - vertical_velocity * du_dy +
                                 parameters.viscosity * laplacian);
        }
    }
    for (int y = 1; y < parameters.ny; ++y) {
        for (int x = 1; x <= parameters.nx; ++x) {
            double horizontal_velocity =
                0.25 * (u(x - 1, y) + u(x, y) + u(x - 1, y + 1) + u(x, y + 1));
            double vertical_velocity = v(x, y);
            double dv_dx;
            if (horizontal_velocity >= 0) {
                dv_dx = (vertical_velocity - v(x - 1, y)) / dx;
            } else {
                dv_dx = (v(x + 1, y) - vertical_velocity) / dx;
            }
            double dv_dy;
            if (vertical_velocity >= 0) {
                dv_dy = (vertical_velocity - v(x, y - 1)) / dy;
            } else {
                dv_dy = (v(x, y + 1) - vertical_velocity) / dy;
            }
            double d2v_dx2 = (v(x + 1, y) - 2 * vertical_velocity + v(x - 1, y)) / (dx * dx);
            double d2v_dy2 = (v(x, y + 1) - 2 * vertical_velocity + v(x, y - 1)) / (dy * dy);
            double laplacian = d2v_dx2 + d2v_dy2;
            double buoyancy_force = 0.0;
            if (parameters.buoyancy && parameters.thermal) {
                double face_temperature = 0.5 * (temperature(x, y) + temperature(x, y + 1));
                buoyancy_force = parameters.gravity * parameters.expansion *
                                 (face_temperature - parameters.reference_temperature);
            }
            predicted_v(x, y) +=
                parameters.dt * (-horizontal_velocity * dv_dx - vertical_velocity * dv_dy +
                                 parameters.viscosity * laplacian + buoyancy_force);
        }
    }
    for (int y = 1; y <= parameters.ny; ++y) {
        if (parameters.left.flow_type == FlowBoundaryType::pressure) {
            double d2u_dx2 = 2 * (u(1, y) - u(0, y)) / (dx * dx);
            double d2u_dy2 = (u(0, y + 1) - 2 * u(0, y) + u(0, y - 1)) / (dy * dy);
            predicted_u(0, y) = u(0, y) + parameters.dt * parameters.viscosity * (d2u_dx2 + d2u_dy2);
        }
        if (parameters.right.flow_type == FlowBoundaryType::pressure) {
            int x = parameters.nx;
            double d2u_dx2 = 2 * (u(x - 1, y) - u(x, y)) / (dx * dx);
            double d2u_dy2 = (u(x, y + 1) - 2 * u(x, y) + u(x, y - 1)) / (dy * dy);
            predicted_u(x, y) = u(x, y) + parameters.dt * parameters.viscosity * (d2u_dx2 + d2u_dy2);
        }
    }
}

void Solver::solve_pressure(const Field &predicted_u, const Field &predicted_v) {
    const Config &parameters = config;
    Field pressure_rhs(parameters.nx + 2, parameters.ny + 2);
    for (int y = 1; y <= parameters.ny; ++y) {
        for (int x = 1; x <= parameters.nx; ++x) {
            pressure_rhs(x, y) = parameters.density / parameters.dt *
                                 ((predicted_u(x, y) - predicted_u(x - 1, y)) / dx +
                                  (predicted_v(x, y) - predicted_v(x, y - 1)) / dy);
        }
    }

    double inverse_dx_squared = 1 / (dx * dx);
    double inverse_dy_squared = 1 / (dy * dy);
    bool converged = false;
    for (int iteration = 0; iteration < parameters.pressure_iterations; ++iteration) {
        for (int y = 1; y <= parameters.ny; ++y) {
            for (int x = 1; x <= parameters.nx; ++x) {
                double sum = 0;
                double diagonal = 0;
                if (x > 1) {
                    sum += inverse_dx_squared * pressure(x - 1, y);
                    diagonal += inverse_dx_squared;
                }
                if (x < parameters.nx) {
                    sum += inverse_dx_squared * pressure(x + 1, y);
                    diagonal += inverse_dx_squared;
                }
                if (x == 1 && parameters.left.flow_type == FlowBoundaryType::pressure) {
                    sum += 2 * inverse_dx_squared * parameters.left.pressure;
                    diagonal += 2 * inverse_dx_squared;
                }
                if (x == parameters.nx && parameters.right.flow_type == FlowBoundaryType::pressure) {
                    sum += 2 * inverse_dx_squared * parameters.right.pressure;
                    diagonal += 2 * inverse_dx_squared;
                }
                if (y > 1) {
                    sum += inverse_dy_squared * pressure(x, y - 1);
                    diagonal += inverse_dy_squared;
                }
                if (y < parameters.ny) {
                    sum += inverse_dy_squared * pressure(x, y + 1);
                    diagonal += inverse_dy_squared;
                }
                double next_pressure = (sum - pressure_rhs(x, y)) / diagonal;
                pressure(x, y) += parameters.pressure_relaxation * (next_pressure - pressure(x, y));
            }
        }
        if (parameters.left.flow_type != FlowBoundaryType::pressure &&
            parameters.right.flow_type != FlowBoundaryType::pressure) {
            double offset = pressure(1, 1);
            for (int y = 1; y <= parameters.ny; ++y) {
                for (int x = 1; x <= parameters.nx; ++x) {
                    pressure(x, y) -= offset;
                }
            }
        }
        bool check_residual = iteration == 0 || (iteration + 1) % 10 == 0;
        bool last_iteration = iteration + 1 == parameters.pressure_iterations;
        if (!check_residual && !last_iteration) {
            continue;
        }
        double residual = 0;
        for (int y = 1; y <= parameters.ny; ++y) {
            for (int x = 1; x <= parameters.nx; ++x) {
                double laplacian = 0;
                if (x > 1) {
                    laplacian += inverse_dx_squared * (pressure(x - 1, y) - pressure(x, y));
                }
                if (x < parameters.nx) {
                    laplacian += inverse_dx_squared * (pressure(x + 1, y) - pressure(x, y));
                }
                if (x == 1 && parameters.left.flow_type == FlowBoundaryType::pressure) {
                    laplacian += 2 * inverse_dx_squared * (parameters.left.pressure - pressure(x, y));
                }
                if (x == parameters.nx && parameters.right.flow_type == FlowBoundaryType::pressure) {
                    laplacian += 2 * inverse_dx_squared * (parameters.right.pressure - pressure(x, y));
                }
                if (y > 1) {
                    laplacian += inverse_dy_squared * (pressure(x, y - 1) - pressure(x, y));
                }
                if (y < parameters.ny) {
                    laplacian += inverse_dy_squared * (pressure(x, y + 1) - pressure(x, y));
                }
                double local_residual =
                    std::abs(laplacian - pressure_rhs(x, y)) * parameters.dt / parameters.density;
                if (!std::isfinite(local_residual)) {
                    throw std::runtime_error("Pressure residual is not finite");
                }
                residual = std::max(residual, local_residual);
            }
        }
        if (residual <= parameters.divergence_tolerance) {
            converged = true;
            break;
        }
    }
    if (!converged) {
        throw std::runtime_error("Pressure did not converge");
    }
}

void Solver::correct_velocity(const Field &predicted_u, const Field &predicted_v) {
    const Config &parameters = config;
    for (int y = 1; y <= parameters.ny; ++y) {
        for (int x = 1; x < parameters.nx; ++x) {
            u(x, y) = predicted_u(x, y) - parameters.dt / parameters.density *
                                              (pressure(x + 1, y) - pressure(x, y)) / dx;
        }
    }
    for (int y = 1; y < parameters.ny; ++y) {
        for (int x = 1; x <= parameters.nx; ++x) {
            v(x, y) = predicted_v(x, y) - parameters.dt / parameters.density *
                                              (pressure(x, y + 1) - pressure(x, y)) / dy;
        }
    }
    for (int y = 1; y <= parameters.ny; ++y) {
        if (parameters.left.flow_type == FlowBoundaryType::inlet) {
            u(0, y) = predicted_u(0, y);
        }
        if (parameters.left.flow_type == FlowBoundaryType::pressure) {
            u(0, y) = predicted_u(0, y) - parameters.dt / parameters.density *
                      (pressure(1, y) - parameters.left.pressure) / (0.5 * dx);
        }
        if (parameters.right.flow_type == FlowBoundaryType::pressure) {
            u(parameters.nx, y) = predicted_u(parameters.nx, y) -
                parameters.dt / parameters.density *
                (parameters.right.pressure - pressure(parameters.nx, y)) / (0.5 * dx);
        }
    }
}

TimeStepLimits Solver::time_step_limits() const {
    double maximum_u = 0.0;
    for (double velocity : u.values) {
        if (!std::isfinite(velocity)) {
            throw std::runtime_error("Horizontal velocity is not finite");
        }
        maximum_u = std::max(maximum_u, std::abs(velocity));
    }

    double maximum_v = 0.0;
    for (double velocity : v.values) {
        if (!std::isfinite(velocity)) {
            throw std::runtime_error("Vertical velocity is not finite");
        }
        maximum_v = std::max(maximum_v, std::abs(velocity));
    }

    TimeStepLimits limits;
    limits.courant_x = maximum_u * config.dt / dx;
    limits.courant_y = maximum_v * config.dt / dy;
    limits.courant = limits.courant_x + limits.courant_y;

    double inverse_dx_squared = 1.0 / (dx * dx);
    double inverse_dy_squared = 1.0 / (dy * dy);
    if (config.fluid) {
        limits.viscous_fourier =
            config.viscosity * config.dt * (inverse_dx_squared + inverse_dy_squared);
        limits.momentum = limits.courant + 2.0 * limits.viscous_fourier;
    }
    if (config.thermal) {
        limits.thermal_fourier =
            config.diffusivity * config.dt * (inverse_dx_squared + inverse_dy_squared);
        limits.thermal = limits.courant + 2.0 * limits.thermal_fourier;
    }
    limits.combined = std::max(limits.momentum, limits.thermal);

    limits.maximum_dt = std::numeric_limits<double>::infinity();
    if (limits.combined > 0.0) {
        limits.maximum_dt = config.dt / limits.combined;
    }

    if (!std::isfinite(limits.courant) || !std::isfinite(limits.momentum) ||
        !std::isfinite(limits.thermal) || !std::isfinite(limits.combined)) {
        throw std::runtime_error("Explicit stability estimate is not finite");
    }
    return limits;
}

void Solver::step() {
    TimeStepLimits limits = time_step_limits();
    if (limits.combined > 1.0) {
        throw std::runtime_error("Explicit stability limit exceeded: CFL + 2*Fo > 1. Reduce dt.");
    }
    // Advance temperature.
    if (config.thermal) {
        update_temperature();
        apply_temperature_boundaries();
    }
    // Advance fluid.
    if (config.fluid) {
        Field predicted_u = u;
        Field predicted_v = v;
        predict_velocity(predicted_u, predicted_v);
        apply_velocity_boundaries(predicted_u, predicted_v, time + config.dt);
        solve_pressure(predicted_u, predicted_v);
        correct_velocity(predicted_u, predicted_v);
    }
    time += config.dt;
    apply_velocity_boundaries(u, v, time);
}

double Solver::maximum_divergence() const {
    double maximum = 0;
    for (int y = 1; y <= config.ny; ++y) {
        for (int x = 1; x <= config.nx; ++x) {
            maximum = std::max(
                maximum, std::abs((u(x, y) - u(x - 1, y)) / dx + (v(x, y) - v(x, y - 1)) / dy));
        }
    }
    return maximum;
}
