#include "cases.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

Config select_case(const std::string &name) {
    if (name == "air") {
        return air_cavity();
    }
    if (name == "water") {
        return water_cavity();
    }
    if (name == "poiseuille") {
        return poiseuille();
    }
    if (name == "couette") {
        return couette();
    }
    if (name == "stokes") {
        return stokes();
    }
    if (name == "heated_wall") {
        return heated_wall();
    }
    if (name == "lid_cavity") {
        return lid_cavity();
    }
    if (name == "buoyancy") {
        return buoyancy();
    }
    if (name == "thermal_bubble") {
        return thermal_bubble();
    }
    throw std::invalid_argument("Unknown simulation case: " + name);
}

void save_snapshot(const Solver &solver, const std::filesystem::path &directory, int step) {
    std::ofstream file(directory / ("step_" + std::to_string(step) + ".csv"));
    if (!file) {
        throw std::runtime_error("Cannot write results");
    }
    file << std::setprecision(16);
    file << "time,x,y,u,v,p,T\n";
    for (int y = 1; y <= solver.config.ny; ++y) {
        for (int x = 1; x <= solver.config.nx; ++x) {
            file << step * solver.config.dt << ',' << (x - 0.5) * solver.dx << ','
                 << (y - 0.5) * solver.dy << ',' << 0.5 * (solver.u(x - 1, y) + solver.u(x, y))
                 << ',' << 0.5 * (solver.v(x, y - 1) + solver.v(x, y)) << ','
                 << solver.pressure(x, y) << ',' << solver.temperature(x, y) << '\n';
        }
    }
    if (!file) {
        throw std::runtime_error("Cannot write results");
    }
}

void save_metadata(const Config &config, const std::string &case_name,
                   const std::filesystem::path &output_directory) {
    std::ofstream metadata(output_directory / "run.json");
    if (!metadata) {
        throw std::runtime_error("Cannot write run metadata");
    }
    metadata << std::setprecision(16) << "{\n"
             << "  \"case\": \"" << case_name << "\",\n"
             << "  \"length\": " << config.length << ",\n"
             << "  \"height\": " << config.height << ",\n"
             << "  \"viscosity\": " << config.viscosity << ",\n"
             << "  \"density\": " << config.density << ",\n"
             << "  \"left_pressure\": " << config.left.pressure << ",\n"
             << "  \"right_pressure\": " << config.right.pressure << ",\n"
             << "  \"top_velocity\": " << config.top.tangential_velocity << ",\n"
             << "  \"diffusivity\": " << config.diffusivity << ",\n"
             << "  \"initial_temperature\": " << config.initial_temperature << ",\n"
             << "  \"left_temperature\": " << config.left.temperature << ",\n"
             << "  \"right_temperature\": " << config.right.temperature << ",\n"
             << "  \"oscillation_amplitude\": " << config.bottom.oscillation_amplitude << ",\n"
             << "  \"oscillation_period\": " << config.bottom.oscillation_period << "\n"
             << "}\n";
    metadata.close();
    if (!metadata) {
        throw std::runtime_error("Cannot write run metadata");
    }
}

void print_simulation_checks(const Solver &solver, const std::string &case_name) {
    std::cout << std::scientific << std::setprecision(3);
    std::cout << "Case: " << case_name << '\n';
    std::cout << "Grid: " << solver.config.nx << " x " << solver.config.ny << "  dx=" << solver.dx
              << " m  dy=" << solver.dy << " m  dt=" << solver.config.dt << " s\n";
    std::cout << "nu=" << solver.config.viscosity << " m^2/s"
              << "  alpha=" << solver.config.diffusivity << " m^2/s" << '\n';
    std::cout << std::boolalpha << "fluid=" << solver.config.fluid
              << "  thermal=" << solver.config.thermal << '\n';
    std::cout << "CFL_x = max|u| * dt/dx; CFL_y = max|v| * dt/dy; CFL = CFL_x + CFL_y\n";
    std::cout << "Fo_nu = nu*dt*(1/dx^2 + 1/dy^2); Fo_alpha = alpha*dt*(1/dx^2 + 1/dy^2)\n";
    std::cout << "S_fluid = CFL + 2*Fo_nu; S_heat = CFL + 2*Fo_alpha.\n";
    std::cout << "Explicit check each step: S <= 1 for each active equation.\n";
    std::cout << "Stored velocity bounds include ghost cells for a conservative estimate.\n";
    std::cout << "dt_max is estimated from the current fields and recomputed each step.\n";
    std::cout << "Disabled equations report zero; dt_max=inf means no transport restriction.\n";
}

void print_time_step_checks(const Solver &solver, int step) {
    TimeStepLimits limits = solver.time_step_limits();
    std::cout << "step=" << step << '/' << solver.config.steps << "  t=" << solver.time << " s"
              << "  div=" << solver.maximum_divergence() << '\n'
              << "  CFL_x=" << limits.courant_x << "  CFL_y=" << limits.courant_y
              << "  CFL=" << limits.courant << '\n'
              << "  Fo_nu=" << limits.viscous_fourier << "  Fo_alpha=" << limits.thermal_fourier
              << "  S_fluid=" << limits.momentum << "  S_heat=" << limits.thermal
              << "  dt_max=" << limits.maximum_dt << " s" << '\n';
}

int main(int argc, char **argv) {
    try {
        if (argc > 3) {
            throw std::invalid_argument("Usage: simulate [case_name] [new_output_directory]");
        }

        std::string case_name = "air";
        if (argc > 1) {
            case_name = argv[1];
        }
        Config config = select_case(case_name);

        std::filesystem::path output_directory = "results";
        if (argc > 2) {
            output_directory = argv[2];
        }
        if (std::filesystem::exists(output_directory) &&
            !std::filesystem::is_empty(output_directory)) {
            throw std::runtime_error("Output directory must be empty");
        }
        std::filesystem::create_directories(output_directory);
        save_metadata(config, case_name, output_directory);

        Solver solver(config);
        print_simulation_checks(solver, case_name);
        print_time_step_checks(solver, 0);
        save_snapshot(solver, output_directory, 0);

        for (int step = 1; step <= config.steps; ++step) {
            solver.step();

            bool save_this_step = step % config.save_every == 0;
            bool last_step = step == config.steps;
            if (save_this_step || last_step) {
                save_snapshot(solver, output_directory, step);
                print_time_step_checks(solver, step);
            }
        }
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
