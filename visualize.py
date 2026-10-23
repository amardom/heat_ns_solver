import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FuncFormatter, MaxNLocator

parser = argparse.ArgumentParser()
parser.add_argument("directory", nargs="?", default="results", type=Path)
parser.add_argument("--field", choices=["T", "p", "speed", "u"])
parser.add_argument("--save", type=Path)
parser.add_argument("--fps", type=int, default=25)
args = parser.parse_args()
if args.fps < 1:
    parser.error("--fps must be positive")
if args.save and args.save.suffix.lower() not in {".png", ".gif", ".mp4"}:
    parser.error("--save requires .png, .gif or .mp4")

files = sorted(
    args.directory.glob("step_*.csv"), key=lambda path: int(path.stem.split("_")[1])
)
if not files:
    parser.error("No step_*.csv files found")
snapshots = []
for path in files:
    data = np.genfromtxt(path, delimiter=",", names=True)
    snapshots.append(data)
metadata_path = args.directory / "run.json"
metadata = {}
if metadata_path.exists():
    metadata = json.loads(metadata_path.read_text())
case_name = metadata.get("case")
if args.field is None:
    if case_name == "stokes":
        args.field = "u"
    elif case_name in {"buoyancy", "heated_wall", "thermal_bubble"}:
        args.field = "T"
    else:
        args.field = "speed"
x_coordinates = np.unique(snapshots[0]["x"])
y_coordinates = np.unique(snapshots[0]["y"])
times = np.array([data["time"][0] for data in snapshots])
dx = x_coordinates[1] - x_coordinates[0]
dy = y_coordinates[1] - y_coordinates[0]
length = len(x_coordinates) * dx
height = len(y_coordinates) * dy
plotted_fields = []
for data in snapshots:
    if args.field == "speed":
        values = np.hypot(data["u"], data["v"])
    else:
        values = data[args.field]
    plotted_fields.append(values.reshape(len(y_coordinates), len(x_coordinates)))
if case_name in {"buoyancy", "heated_wall", "thermal_bubble"} and args.field == "T":
    for frame in range(len(plotted_fields)):
        plotted_fields[frame] = plotted_fields[frame] - 273.15
field_minimum = min(float(values.min()) for values in plotted_fields)
field_maximum = max(float(values.max()) for values in plotted_fields)
if args.field == "u":
    field_maximum = max(abs(field_minimum), abs(field_maximum))
    field_minimum = -field_maximum
if field_minimum == field_maximum:
    field_maximum = field_minimum + 1e-12

show_velocity_profile = case_name in {"poiseuille", "couette"} and args.field == "speed"
show_stokes = case_name == "stokes" and args.field == "u"
show_velocity_profile = show_velocity_profile or show_stokes
show_temperature_profile = case_name == "heated_wall" and args.field == "T"
plt.style.use("dark_background")
plt.rcParams.update(
    {
        "figure.facecolor": "#101923",
        "axes.facecolor": "#101923",
        "axes.edgecolor": "#667785",
        "text.color": "#e3edf5",
        "axes.labelcolor": "#c0ced9",
        "xtick.color": "#a6bac9",
        "ytick.color": "#a6bac9",
        "font.size": 10,
    }
)
if show_velocity_profile or show_temperature_profile:
    figure, (axes, profile_axes) = plt.subplots(
        1,
        2,
        figsize=(12, 4.4),
        gridspec_kw={"width_ratios": [2.8, 1]},
        layout="constrained",
    )
else:
    if case_name == "thermal_bubble":
        size = (6, 7.5)
    elif case_name == "lid_cavity":
        size = (6.5, 6)
    else:
        size = (10, 5)
    figure, axes = plt.subplots(figsize=size, layout="constrained")

label = {"T": "T [K]", "p": "p [Pa]", "speed": "|u| [m/s]", "u": "u [m/s]"}[args.field]
if case_name in {"buoyancy", "heated_wall", "thermal_bubble"} and args.field == "T":
    label = "T [°C]"
coordinate_scale = 1
unit = "m"
if show_temperature_profile:
    coordinate_scale = 1000
    unit = "mm"
field_colormap = "magma"
if args.field == "u":
    field_colormap = "coolwarm"
if show_stokes:
    field_colormap = LinearSegmentedColormap.from_list(
        "dark_velocity",
        [
            (0.00, "#69d8e5"),
            (0.20, "#3958a7"),
            (0.40, "#171b43"),
            (0.50, "#000004"),
            (0.60, "#24104d"),
            (0.80, "#b63679"),
            (0.93, "#fb8761"),
            (1.00, "#fcfdbf"),
        ],
    )
image = axes.imshow(
    plotted_fields[0],
    origin="lower",
    aspect="equal",
    interpolation="bilinear",
    extent=(0, length * coordinate_scale, 0, height * coordinate_scale),
    vmin=field_minimum,
    vmax=field_maximum,
    cmap=field_colormap,
)
figure.colorbar(image, ax=axes, label=label, shrink=0.75, pad=0.025)
axes.set(xlabel=f"x [{unit}]", ylabel=f"y [{unit}]")
if case_name == "thermal_bubble":
    axes.xaxis.set_major_locator(MaxNLocator(5))
    axes.yaxis.set_major_locator(MaxNLocator(6))
    axes.xaxis.set_major_formatter(
        FuncFormatter(lambda value, position: f"{1000 * value:g}")
    )
    axes.yaxis.set_major_formatter(
        FuncFormatter(lambda value, position: f"{1000 * value:g}")
    )
    axes.set(xlabel="x [mm]", ylabel="y [mm]")
title = figure.suptitle("", fontsize=16, fontweight="bold")

if show_velocity_profile:
    viscosity = metadata["viscosity"]
    profile_heights = np.linspace(0, height, 250)
    if show_stokes:
        amplitude = metadata["oscillation_amplitude_stokes"]
        omega = 2 * np.pi / metadata["oscillation_period_stokes"]
        stokes_wave = np.sqrt(1j * omega / viscosity)
        stokes_shape = np.sinh(stokes_wave * (height - profile_heights)) / np.sinh(
            stokes_wave * height
        )
        steady = np.imag(amplitude * stokes_shape)
    elif case_name == "couette":
        wall_velocity = metadata["top_velocity"]
        steady = wall_velocity * profile_heights / height
    else:
        pressure_drop = metadata["left_pressure"] - metadata["right_pressure"]
        acceleration = pressure_drop / (metadata["density"] * length)
        steady = (
            acceleration
            * profile_heights
            * (height - profile_heights)
            / (2 * viscosity)
        )
    reference_label = "Steady"
    if show_stokes:
        reference_label = "Time-periodic"
    (reference_line,) = profile_axes.plot(
        steady, profile_heights, "--", color="#7a8b98", label=reference_label
    )
    (exact_line,) = profile_axes.plot(
        profile_heights * 0, profile_heights, color="#ffbd69", label="Analytical"
    )
    (numerical_line,) = profile_axes.plot(
        snapshots[0]["u"].reshape(len(y_coordinates), len(x_coordinates))[
            :, len(x_coordinates) // 2
        ],
        y_coordinates,
        "o",
        markersize=3,
        color="#62d8e5",
        label="C++",
    )
    if show_stokes:
        profile_limits = (-amplitude * 1.1, amplitude * 1.1)
    else:
        profile_limits = (0, steady.max() * 1.08)
    profile_axes.set(
        xlim=profile_limits,
        ylim=(0, height),
        xlabel="u [m/s]",
        ylabel="y [m]",
        title="Velocity profile",
    )
    profile_axes.grid(alpha=0.12)
    if show_stokes:
        profile_axes.axvline(0, color="#667785", linewidth=0.8)
    profile_axes.legend(frameon=False, fontsize=8, loc="lower right")

    generator = np.random.default_rng(7)
    particle_y = generator.uniform(dy / 2, height - dy / 2, 100)
    particle_x = generator.uniform(0, length, len(particle_y))
    particle_positions = [np.column_stack((particle_x, particle_y))]
    for frame in range(1, len(snapshots)):
        previous_u = (
            snapshots[frame - 1]["u"]
            .reshape(len(y_coordinates), len(x_coordinates))
            .mean(axis=1)
        )
        current_u = (
            snapshots[frame]["u"]
            .reshape(len(y_coordinates), len(x_coordinates))
            .mean(axis=1)
        )
        velocity = np.interp(particle_y, y_coordinates, 0.5 * (previous_u + current_u))
        particle_x = (
            particle_x + velocity * (times[frame] - times[frame - 1])
        ) % length
        particle_positions.append(np.column_stack((particle_x, particle_y)))
    particle_size = 7
    if show_stokes:
        particle_size = 10
    particles = axes.scatter(
        particle_positions[0][:, 0],
        particle_y,
        s=particle_size,
        c="white",
        alpha=0.8,
        linewidths=0,
    )

if show_temperature_profile:
    heat_positions = np.linspace(0, length, 250)
    cold = metadata["right_temperature"] - 273.15
    hot = metadata["left_temperature"] - 273.15
    heat_difference = hot - cold
    heat_steady = cold + heat_difference * (1 - heat_positions / length)
    profile_axes.plot(
        heat_positions * 1000, heat_steady, "--", color="#7a8b98", label="Steady"
    )
    (heat_exact_line,) = profile_axes.plot(
        heat_positions * 1000,
        heat_positions * 0 + cold,
        color="#ffbd69",
        label="Analytical",
    )
    (heat_numerical_line,) = profile_axes.plot(
        x_coordinates * 1000,
        plotted_fields[0][len(y_coordinates) // 2, :],
        "o",
        markersize=3,
        color="#62d8e5",
        label="C++",
    )
    profile_axes.set(
        xlim=(0, length * 1000),
        ylim=(cold - 0.5, hot + 0.5),
        xlabel="x [mm]",
        ylabel="T [°C]",
        title="Temperature profile",
    )
    profile_axes.grid(alpha=0.12)
    profile_axes.legend(frameon=False, fontsize=8, loc="upper right")
    axes.plot([0, 0], [0, height * 1000], color="#ff7b54", linewidth=4)
    axes.plot(
        [length * 1000, length * 1000], [0, height * 1000], color="#62d8e5", linewidth=4
    )


def analytical_temperature(time):
    if time == 0:
        return np.full_like(heat_positions, cold)
    temperature = heat_steady.copy()
    for n in range(1, 202):
        wave_number = n * np.pi / length
        temperature -= (
            2
            * heat_difference
            / (n * np.pi)
            * np.sin(wave_number * heat_positions)
            * np.exp(-metadata["diffusivity"] * wave_number**2 * time)
        )
    return temperature


show_vectors = case_name in {
    "lid_cavity",
    "buoyancy",
    "thermal_bubble",
} and args.field in {"speed", "T"}
if show_vectors:
    tracer_speed = 0.0
    for data in snapshots:
        speed = np.hypot(data["u"], data["v"])
        tracer_speed = max(tracer_speed, float(speed.max()))
    tracer_speed = max(tracer_speed, metadata["top_velocity"], 1e-12)
    stride = 2
    vector_x, vector_y = np.meshgrid(x_coordinates[::stride], y_coordinates[::stride])
    vectors = axes.quiver(
        vector_x,
        vector_y,
        np.zeros_like(vector_x),
        np.zeros_like(vector_y),
        color="white",
        alpha=0.3,
        angles="xy",
        scale_units="xy",
        scale=tracer_speed / length * 14,
        width=0.003,
    )
    if case_name == "lid_cavity":
        axes.annotate(
            "",
            xy=(0.9 * length, height * 1.025),
            xytext=(0.1 * length, height * 1.025),
            arrowprops={"arrowstyle": "->", "color": "#62d8e5", "lw": 2},
            annotation_clip=False,
        )
    elif case_name == "buoyancy":
        axes.plot([0, length], [0, 0], color="#ff7b54", linewidth=4)
        axes.plot([0, length], [height, height], color="#62d8e5", linewidth=4)

    generator = np.random.default_rng(7)
    particle_x = generator.uniform(dx, length - dx, 180)
    particle_y = generator.uniform(dy, height - dy, 180)
    playback_label = axes.text(
        0.02,
        0.03,
        "",
        transform=axes.transAxes,
        fontsize=9,
        bbox={"facecolor": "#101923", "alpha": 0.8, "edgecolor": "none"},
    )
    cavity_particles = np.column_stack((particle_x, particle_y))
    particle_colors = ["white"] * len(cavity_particles)
    if case_name == "thermal_bubble":
        initial_temperature = snapshots[0]["T"].reshape(
            len(y_coordinates), len(x_coordinates)
        )
        threshold = initial_temperature.min() + 0.5 * np.ptp(initial_temperature)
        warm_cells = np.argwhere(initial_temperature > threshold)
        selected = warm_cells[generator.integers(0, len(warm_cells), 80)]
        warm_particles = np.column_stack(
            (
                x_coordinates[selected[:, 1]]
                + generator.uniform(-dx / 2, dx / 2, len(selected)),
                y_coordinates[selected[:, 0]]
                + generator.uniform(-dy / 2, dy / 2, len(selected)),
            )
        )
        cavity_particles = np.vstack((cavity_particles, warm_particles))
        particle_colors += ["#101923"] * len(warm_particles)
    cavity_positions = [cavity_particles.copy()]
    sample_x = np.concatenate(([0], x_coordinates, [length]))
    sample_y = np.concatenate(([0], y_coordinates, [height]))

    def padded_velocity(data):
        velocity_u = np.pad(
            data["u"].reshape(len(y_coordinates), len(x_coordinates)), 1
        )
        velocity_v = np.pad(
            data["v"].reshape(len(y_coordinates), len(x_coordinates)), 1
        )
        velocity_u[-1, 1:-1] = metadata["top_velocity"]
        return velocity_u, velocity_v

    def interpolate(values, points):
        particle_x = np.clip(points[:, 0], 0, length)
        particle_y = np.clip(points[:, 1], 0, height)
        left_column = np.clip(
            np.searchsorted(sample_x, particle_x) - 1, 0, len(sample_x) - 2
        )
        bottom_row = np.clip(
            np.searchsorted(sample_y, particle_y) - 1, 0, len(sample_y) - 2
        )
        horizontal_weight = (particle_x - sample_x[left_column]) / (
            sample_x[left_column + 1] - sample_x[left_column]
        )
        vertical_weight = (particle_y - sample_y[bottom_row]) / (
            sample_y[bottom_row + 1] - sample_y[bottom_row]
        )
        return (
            (1 - horizontal_weight)
            * (1 - vertical_weight)
            * values[bottom_row, left_column]
            + horizontal_weight
            * (1 - vertical_weight)
            * values[bottom_row, left_column + 1]
            + (1 - horizontal_weight)
            * vertical_weight
            * values[bottom_row + 1, left_column]
            + horizontal_weight
            * vertical_weight
            * values[bottom_row + 1, left_column + 1]
        )

    previous_u, previous_v = padded_velocity(snapshots[0])
    for frame in range(1, len(snapshots)):
        current_u, current_v = padded_velocity(snapshots[frame])
        interval = times[frame] - times[frame - 1]
        substeps = max(1, int(np.ceil(interval * tracer_speed / (0.3 * min(dx, dy)))))
        particle_dt = interval / substeps
        for substep in range(substeps):
            fraction = (substep + 0.5) / substeps
            velocity_u = (1 - fraction) * previous_u + fraction * current_u
            velocity_v = (1 - fraction) * previous_v + fraction * current_v
            velocity = np.column_stack(
                (
                    interpolate(velocity_u, cavity_particles),
                    interpolate(velocity_v, cavity_particles),
                )
            )
            midpoint = cavity_particles + 0.5 * particle_dt * velocity
            midpoint_velocity = np.column_stack(
                (interpolate(velocity_u, midpoint), interpolate(velocity_v, midpoint))
            )
            cavity_particles += particle_dt * midpoint_velocity
            cavity_particles[:, 0] = np.clip(cavity_particles[:, 0], 0, length)
            cavity_particles[:, 1] = np.clip(cavity_particles[:, 1], 0, height)
        cavity_positions.append(cavity_particles.copy())
        previous_u, previous_v = current_u, current_v
    moving_particles = axes.scatter(
        cavity_positions[0][:, 0],
        cavity_positions[0][:, 1],
        s=9,
        color=particle_colors,
        alpha=0.9,
        linewidths=0,
    )


def analytical_profile(time):
    if time == 0:
        return np.zeros_like(profile_heights)
    if show_stokes:
        velocity = amplitude * (1 - profile_heights / height) * np.sin(omega * time)
        for n in range(1, 302):
            wave_number = n * np.pi / height
            decay = viscosity * wave_number**2
            coefficient = (
                -2
                * amplitude
                * omega
                / (n * np.pi)
                * (
                    decay * np.cos(omega * time)
                    + omega * np.sin(omega * time)
                    - decay * np.exp(-decay * time)
                )
                / (decay**2 + omega**2)
            )
            velocity += coefficient * np.sin(wave_number * profile_heights)
        return velocity
    velocity = steady.copy()
    if case_name == "couette":
        for n in range(1, 202):
            wave_number = n * np.pi / height
            velocity += (
                2
                * wall_velocity
                * (-1) ** n
                / (n * np.pi)
                * np.sin(wave_number * profile_heights)
                * np.exp(-viscosity * wave_number**2 * time)
            )
        return velocity
    for n in range(1, 202, 2):
        wave_number = n * np.pi / height
        velocity -= (
            4
            * acceleration
            * height**2
            / (viscosity * np.pi**3 * n**3)
            * np.sin(wave_number * profile_heights)
            * np.exp(-viscosity * wave_number**2 * time)
        )
    return velocity


def update(frame):
    left = int(frame)
    right = min(left + 1, len(files) - 1)
    fraction = frame - left

    def blend(values):
        return (1 - fraction) * values[left] + fraction * values[right]

    image.set_data(blend(plotted_fields))
    name = metadata.get("case", args.directory.name).replace("_", " ").capitalize()
    time = (1 - fraction) * times[left] + fraction * times[right]
    title.set_text(f"{name}   |   t = {time:.2f} s")
    if show_vectors:
        velocity_u = (
            (1 - fraction) * snapshots[left]["u"] + fraction * snapshots[right]["u"]
        ).reshape(len(y_coordinates), len(x_coordinates))
        velocity_v = (
            (1 - fraction) * snapshots[left]["v"] + fraction * snapshots[right]["v"]
        ).reshape(len(y_coordinates), len(x_coordinates))
        vectors.set_UVC(velocity_u[::stride, ::stride], velocity_v[::stride, ::stride])
        moving_particles.set_offsets(blend(cavity_positions))
        if time == 0:
            playback_text = "At rest"
        elif case_name in {"buoyancy", "thermal_bubble"}:
            playback_text = "Slow motion"
        elif time < 5:
            playback_text = "Startup · real time"
        else:
            playback_text = "Evolution · 5x"
        playback_label.set_text(playback_text)
    if show_temperature_profile:
        heat_numerical_line.set_ydata(blend(plotted_fields)[len(y_coordinates) // 2, :])
        heat_exact_line.set_ydata(analytical_temperature(time))
    if show_velocity_profile:
        numerical_line.set_xdata(
            (
                (1 - fraction) * snapshots[left]["u"] + fraction * snapshots[right]["u"]
            ).reshape(len(y_coordinates), len(x_coordinates))[
                :, len(x_coordinates) // 2
            ]
        )
        exact_line.set_xdata(analytical_profile(time))
        if show_stokes:
            reference_line.set_xdata(
                np.imag(amplitude * stokes_shape * np.exp(1j * omega * time))
            )
            particle_frame = particle_positions[left].copy()
            displacement = (
                particle_positions[right][:, 0] - particle_frame[:, 0] + length / 2
            ) % length - length / 2
            particle_frame[:, 0] = (
                particle_frame[:, 0] + fraction * displacement
            ) % length
            particles.set_offsets(particle_frame)
        else:
            particles.set_offsets(particle_positions[left])


if args.save and args.save.suffix.lower() == ".png":
    update(len(files) - 1)
    figure.savefig(args.save, dpi=160)
else:
    update(0)
    figure.canvas.draw()
    figure.set_layout_engine("none")
    animation_frames = range(len(files))
    if show_temperature_profile or show_stokes:
        duration = 15
        if show_stokes:
            duration = 20
        playback_times = np.linspace(0, times[-1], duration * args.fps)
        animation_frames = np.interp(playback_times, times, np.arange(len(files)))
    if show_vectors:
        if case_name in {"buoyancy", "thermal_bubble"}:
            playback_times = np.linspace(0, min(8.0, times[-1]), 15 * args.fps)
        else:
            startup_end = min(5.0, times[-1])
            playback_times = np.concatenate(
                (
                    np.zeros(args.fps),
                    np.arange(0, startup_end, 1 / args.fps),
                    np.arange(startup_end, times[-1], 5 / args.fps),
                    [times[-1]],
                )
            )
        animation_frames = np.interp(playback_times, times, np.arange(len(files)))
    animation = FuncAnimation(
        figure, update, frames=animation_frames, interval=1000 / args.fps, repeat=True
    )
    if args.save:
        if args.save.suffix.lower() == ".mp4":
            import imageio_ffmpeg

            mpl.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
            writer = FFMpegWriter(
                fps=args.fps,
                codec="libx264",
                extra_args=["-pix_fmt", "yuv420p", "-crf", "20"],
            )
        else:
            writer = PillowWriter(fps=args.fps)
        animation.save(args.save, writer=writer, dpi=120)
    else:
        plt.show()
