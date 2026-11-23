from dataclasses import dataclass
from pathlib import Path
from typing import Final

import numpy as np
import numpy.typing as npt
from matplotlib import colormaps
from matplotlib.animation import PillowWriter
from matplotlib.figure import Figure

FloatArray = npt.NDArray[np.float64]
IndexArray = npt.NDArray[np.intp]

LATTICE_VELOCITIES: Final = ((0, 0), (1, 0), (0, 1), (-1, 0), (0, -1), (1, 1), (-1, 1), (-1, -1), (1, -1))
VELOCITY_X: Final = np.array([step_x for step_x, _ in LATTICE_VELOCITIES], dtype=np.float64)
VELOCITY_Y: Final = np.array([step_y for _, step_y in LATTICE_VELOCITIES], dtype=np.float64)
WEIGHTS: Final = np.array([4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36], dtype=np.float64)
OPPOSITE: Final = np.array([LATTICE_VELOCITIES.index((-step_x, -step_y)) for step_x, step_y in LATTICE_VELOCITIES])
REST: Final = 0
MOVING_EAST: Final = np.flatnonzero(VELOCITY_X > 0)
MOVING_WEST: Final = np.flatnonzero(VELOCITY_X < 0)
MOVING_UP: Final = np.flatnonzero(VELOCITY_Y > 0)
MOVING_DOWN: Final = np.flatnonzero(VELOCITY_Y < 0)
MOVING_ALONG_LID: Final = np.flatnonzero(VELOCITY_Y == 0)
COLOR_LEVELS: Final = 32


@dataclass(frozen=True)
class Parameters:
    cavity_size: int = 100
    lid_velocity: float = 0.2
    viscosity: float = 0.02
    prandtl_number: float = 0.71
    lid_temperature: float = 1.0
    side_wall_temperature: float = 0.0
    initial_density: float = 5.0
    steps: int = 20000
    snapshot_interval: int = 200

    @property
    def thermal_diffusivity(self) -> float:
        return self.viscosity / self.prandtl_number

    @property
    def flow_relaxation_rate(self) -> float:
        return 1 / (3 * self.viscosity + 0.5)

    @property
    def heat_relaxation_rate(self) -> float:
        return 1 / (3 * self.thermal_diffusivity + 0.5)

    @property
    def reynolds_number(self) -> float:
        return self.lid_velocity * self.cavity_size / self.viscosity

    @property
    def peclet_number(self) -> float:
        return self.lid_velocity * self.cavity_size / self.thermal_diffusivity


@dataclass(frozen=True)
class Snapshot:
    step: int
    speed: FloatArray
    temperature: FloatArray


def projected_velocity(velocity_x: FloatArray, velocity_y: FloatArray) -> FloatArray:
    return VELOCITY_X[:, None, None] * velocity_x + VELOCITY_Y[:, None, None] * velocity_y


def flow_equilibrium(density: FloatArray, velocity_x: FloatArray, velocity_y: FloatArray) -> FloatArray:
    projected = projected_velocity(velocity_x, velocity_y)
    speed_squared = velocity_x**2 + velocity_y**2
    return WEIGHTS[:, None, None] * density * (1 + 3 * projected + 4.5 * projected**2 - 1.5 * speed_squared)


def heat_equilibrium(temperature: FloatArray, velocity_x: FloatArray, velocity_y: FloatArray) -> FloatArray:
    return WEIGHTS[:, None, None] * temperature * (1 + 3 * projected_velocity(velocity_x, velocity_y))


def shifted_slices(shift: int) -> tuple[slice, slice]:
    if shift > 0:
        return slice(shift, None), slice(None, -shift)
    if shift < 0:
        return slice(None, shift), slice(-shift, None)
    return slice(None), slice(None)


def stream(populations: FloatArray) -> FloatArray:
    streamed = populations.copy()
    for direction, (step_x, step_y) in enumerate(LATTICE_VELOCITIES):
        target_rows, source_rows = shifted_slices(-step_y)
        target_columns, source_columns = shifted_slices(step_x)
        streamed[direction, target_rows, target_columns] = populations[direction, source_rows, source_columns]
    return streamed


def apply_flow_walls(populations: FloatArray, lid_velocity: float) -> None:
    populations[MOVING_EAST, :, 0] = populations[OPPOSITE[MOVING_EAST], :, 0]
    populations[MOVING_WEST, :, -1] = populations[OPPOSITE[MOVING_WEST], :, -1]
    populations[MOVING_UP, -1, :] = populations[OPPOSITE[MOVING_UP], -1, :]
    lid = populations[:, 0, :]
    lid_density = lid[MOVING_ALONG_LID].sum(axis=0) + 2 * lid[MOVING_UP].sum(axis=0)
    lid_push = 6 * (WEIGHTS * VELOCITY_X)[MOVING_DOWN, None] * lid_velocity * lid_density
    lid[MOVING_DOWN] = lid[OPPOSITE[MOVING_DOWN]] + lid_push


def set_wall_temperature(wall: FloatArray, incoming: IndexArray, temperature: float) -> None:
    wall[incoming] = 2 * WEIGHTS[incoming, None] * temperature - wall[OPPOSITE[incoming]]
    wall[REST] = WEIGHTS[REST] * temperature


def apply_heat_walls(populations: FloatArray, lid_temperature: float, side_wall_temperature: float) -> None:
    set_wall_temperature(populations[:, :, 0], MOVING_EAST, side_wall_temperature)
    set_wall_temperature(populations[:, :, -1], MOVING_WEST, side_wall_temperature)
    set_wall_temperature(populations[:, 0, :], MOVING_DOWN, lid_temperature)
    populations[:, -1, :] = populations[:, -2, :]


def status_line(step: int, density: FloatArray, initial_mass: float, latest: Snapshot, previous: Snapshot) -> str:
    return (
        f"step {step:6d}  density {density.min():.4f} to {density.max():.4f}  "
        f"mass change {density.sum() / initial_mass - 1:+.3%}  max speed {latest.speed.max():.4f}  "
        f"temperature {latest.temperature.min():.4f} to {latest.temperature.max():.4f}  "
        f"largest speed change since the last frame {np.abs(latest.speed - previous.speed).max():.1e}"
    )


def simulate(parameters: Parameters) -> list[Snapshot]:
    size = parameters.cavity_size + 1
    density = np.full((size, size), parameters.initial_density, dtype=np.float64)
    velocity_x = np.zeros((size, size), dtype=np.float64)
    velocity_y = np.zeros((size, size), dtype=np.float64)
    velocity_x[0, :] = parameters.lid_velocity
    temperature = np.zeros((size, size), dtype=np.float64)
    flow = WEIGHTS[:, None, None] * density
    heat = WEIGHTS[:, None, None] * temperature
    initial_mass = float(density.sum())
    snapshots = [Snapshot(0, np.hypot(velocity_x, velocity_y), temperature.copy())]

    for step in range(1, parameters.steps + 1):
        flow_rate = parameters.flow_relaxation_rate
        flow = (1 - flow_rate) * flow + flow_rate * flow_equilibrium(density, velocity_x, velocity_y)
        flow = stream(flow)
        apply_flow_walls(flow, parameters.lid_velocity)
        density = flow.sum(axis=0)
        velocity_x = (VELOCITY_X[:, None, None] * flow).sum(axis=0) / density
        velocity_y = (VELOCITY_Y[:, None, None] * flow).sum(axis=0) / density

        heat_rate = parameters.heat_relaxation_rate
        heat = (1 - heat_rate) * heat + heat_rate * heat_equilibrium(temperature, velocity_x, velocity_y)
        heat = stream(heat)
        apply_heat_walls(heat, parameters.lid_temperature, parameters.side_wall_temperature)
        temperature = heat.sum(axis=0)

        if step % parameters.snapshot_interval == 0:
            if not (np.isfinite(density).all() and np.isfinite(temperature).all()):
                raise FloatingPointError(f"the simulation diverged before step {step}")
            snapshots.append(Snapshot(step, np.hypot(velocity_x, velocity_y), temperature.copy()))
            if step % (10 * parameters.snapshot_interval) == 0:
                print(status_line(step, density, initial_mass, snapshots[-1], snapshots[-2]))

    return snapshots


def render(snapshots: list[Snapshot], parameters: Parameters, path: Path) -> None:
    figure = Figure(figsize=(8, 3.9), layout="constrained")
    speed_axis = figure.add_subplot(1, 2, 1)
    temperature_axis = figure.add_subplot(1, 2, 2)
    speed_image = speed_axis.imshow(
        snapshots[0].speed / parameters.lid_velocity,
        cmap=colormaps["Blues"].resampled(COLOR_LEVELS),
        vmin=0,
        vmax=1,
    )
    temperature_image = temperature_axis.imshow(
        snapshots[0].temperature,
        cmap=colormaps["inferno"].resampled(COLOR_LEVELS),
        vmin=parameters.side_wall_temperature,
        vmax=parameters.lid_temperature,
    )
    figure.colorbar(speed_image, ax=speed_axis, label="speed / lid speed")
    figure.colorbar(temperature_image, ax=temperature_axis, label="temperature")
    speed_axis.set_title("Flow speed")
    temperature_axis.set_title("Temperature")
    for axis in (speed_axis, temperature_axis):
        axis.set_xticks([])
        axis.set_yticks([])
    title = figure.suptitle("")

    writer = PillowWriter(fps=10)
    with writer.saving(figure, str(path), dpi=80):
        for snapshot in snapshots:
            speed_image.set_data(snapshot.speed / parameters.lid_velocity)
            temperature_image.set_data(snapshot.temperature)
            title.set_text(f"Heated lid-driven cavity, Re {parameters.reynolds_number:.0f}, step {snapshot.step}")
            writer.grab_frame()


def main() -> None:
    parameters = Parameters()
    size = parameters.cavity_size + 1
    print(
        f"Re {parameters.reynolds_number:.0f}, Pr {parameters.prandtl_number}, Pe {parameters.peclet_number:.0f}, "
        f"grid {size}x{size}, {parameters.steps} steps"
    )
    snapshots = simulate(parameters)
    output = Path("cavity.gif")
    render(snapshots, parameters, output)
    print(f"wrote {output}: {len(snapshots)} frames, {output.stat().st_size / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
