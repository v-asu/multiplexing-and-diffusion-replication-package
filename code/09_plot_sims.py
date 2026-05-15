################# Bootstrap ###############################
import importlib
import os
from pathlib import Path
import re
import site
import subprocess
import sys


REQUIRED_PACKAGES = {
    "numpy": "numpy",
    "pandas": "pandas",
    "matplotlib": "matplotlib",
    "scipy": "scipy",
}


def ensure_python_packages():
    user_site = site.getusersitepackages()
    if user_site not in sys.path:
        site.addsitedir(user_site)

    missing = []

    for module_name, package_name in REQUIRED_PACKAGES.items():
        try:
            importlib.import_module(module_name)
        except ModuleNotFoundError:
            missing.append(package_name)

    if not missing:
        return

    print(
        "Installing missing Python packages: " + ", ".join(missing),
        file=sys.stderr,
    )

    install_commands = [
        [sys.executable, "-m", "pip", "install", "--user", *missing],
        [
            sys.executable,
            "-m",
            "pip",
            "install",
            "--user",
            "--break-system-packages",
            *missing,
        ],
    ]

    last_error = None

    for command in install_commands:
        try:
            subprocess.check_call(command)
            site.addsitedir(user_site)
            return
        except subprocess.CalledProcessError as exc:
            last_error = exc

    raise RuntimeError(
        "Failed to install required Python packages: " + ", ".join(missing)
    ) from last_error


def find_package_root():
    script_root = Path(__file__).resolve().parent.parent
    candidates = [
        Path.cwd(),
        Path.cwd() / "replication_package",
        script_root,
    ]

    for candidate in candidates:
        if (candidate / "data" / "raw" / "sims").is_dir():
            return candidate

    raise FileNotFoundError("Could not locate replication_package root.")


ensure_python_packages()

import matplotlib  # noqa: E402

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.interpolate import griddata  # noqa: E402


PACKAGE_ROOT = find_package_root()
SIMS_DIR = PACKAGE_ROOT / "data" / "raw" / "sims"
FIGURES_DIR = PACKAGE_ROOT / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


################# Functions ###############################

def get_interpolation(df, outcome):
    df_pivot = (
        df.groupby(["q", "delta"])
        .mean(numeric_only=True)
        .reset_index()
        .pivot(index="q", columns="delta", values=outcome)
    )

    q_un = df_pivot.index.values
    delta_un = df_pivot.columns.values
    q, delta = np.meshgrid(q_un, delta_un)

    y_mean = df_pivot.values.T

    q_grid, delta_grid = np.meshgrid(
        np.linspace(q_un.min(), q_un.max(), 100),
        np.linspace(delta_un.min(), delta_un.max(), 100),
    )

    y_mean_interp = griddata(
        (q.flatten(), delta.flatten()),
        y_mean.flatten(),
        (q_grid, delta_grid),
        method="cubic",
    )

    return q_grid, delta_grid, y_mean_interp


def mapping_level_curves(grid_1, grid_2, num_levels=10):
    valid_grid_2 = grid_2[2][np.isfinite(grid_2[2])]
    if valid_grid_2.size == 0:
        raise ValueError("Interpolation grid for infection is entirely missing.")

    min_inf = valid_grid_2.min()
    max_inf = valid_grid_2.max()

    equidistant_values = np.linspace(min_inf, max_inf, num=num_levels + 2)[1:-1]
    sorted_values = np.sort(valid_grid_2.flatten())
    equidistant_indices = np.searchsorted(sorted_values, equidistant_values)
    equidistant_indices = np.clip(equidistant_indices, 0, len(sorted_values) - 1)
    equidistant_points = sorted_values[equidistant_indices]

    out = np.full((num_levels, 2), np.nan)

    for k, level in enumerate(equidistant_points):
        contour_set = plt.contour(
            grid_2[0], grid_2[1], grid_2[2], levels=[level], colors="r"
        )

        contour_points = None
        if hasattr(contour_set, "collections"):
            for path_collection in contour_set.collections:
                paths = path_collection.get_paths()
                if paths:
                    contour_points = paths[-1].vertices
        elif hasattr(contour_set, "get_paths"):
            paths = contour_set.get_paths()
            if paths:
                contour_points = paths[-1].vertices

        plt.close()

        if contour_points is None or len(contour_points) == 0:
            continue

        interpolated_values_at_level = griddata(
            (grid_1[0].flatten(), grid_1[1].flatten()),
            grid_1[2].flatten(),
            contour_points,
            method="cubic",
        )

        if interpolated_values_at_level is None:
            continue

        out[k, 0] = level
        out[k, 1] = np.nanmean(interpolated_values_at_level)

    return out[np.isfinite(out).all(axis=1)]


def master_function(data):
    """
    data: a DataFrame with columns village, q, delta, mpex1, inf1, mpex2, inf2
    """
    mpex_grid = get_interpolation(data, "binary_inf")
    infection_grid = get_interpolation(data, "inf2")
    return mapping_level_curves(mpex_grid, infection_grid, num_levels=10)


def bootstrap_group(group):
    return group.sample(frac=1, replace=True)


def generate_bootstrap_sample(data):
    df = data.groupby("village", group_keys=False).apply(bootstrap_group)
    df.reset_index(drop=True, inplace=True)
    return df


def load_simulation_dataframe(files):
    arrays = [
        np.hstack(
            (
                np.full((data.shape[0], 1), int(re.search(r"\d+", file.name).group())),
                data,
            )
        )
        for file, data in ((file, np.load(file)) for file in files)
    ]

    return pd.DataFrame(
        np.vstack(arrays),
        columns=["village", "q", "delta", "mpex1", "inf1", "mpex2", "inf2"],
    )


def add_outcome_columns(df):
    df = df.copy()
    df["mpex_diff"] = df["mpex2"] - df["mpex1"]
    df["inf_diff"] = df["inf2"] - df["inf1"]
    df["binary_inf"] = np.where((df["mpex_diff"] * df["inf_diff"]) > 0, 1, 0)
    return df


def save_scatter(points, output_path, ymin=None, ymax=None):
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(points[:, 0], points[:, 1], marker="o")
    ax.axhline(y=0.5, color="r", linestyle="--")
    ax.set_ylabel("Fraction of simulations where \n more multiplexing leads to more diffusion")
    ax.set_xlabel("Extent of Diffusion ($p$)")

    if ymin is not None or ymax is not None:
        ax.set_ylim(ymin, ymax)

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


################# Data ####################################

file_list = sorted(
    SIMS_DIR.iterdir(),
    key=lambda path: (
        0 if "simple" in path.name else 1,
        int(re.search(r"\d+", path.name).group()),
    ),
)

simple_files = [file for file in file_list if "simple" in file.name]
complex_files = [file for file in file_list if "simple" not in file.name]

df_simple = add_outcome_columns(load_simulation_dataframe(simple_files))
df_complex = add_outcome_columns(load_simulation_dataframe(complex_files))


################## Interpolation ##########################

full_sample_simple = master_function(df_simple)
save_scatter(full_sample_simple, FIGURES_DIR / "fig_main_03a_simple_contagion_simulation.pdf")

full_sample_complex = master_function(df_complex)
save_scatter(
    full_sample_complex,
    FIGURES_DIR / "fig_main_03b_complex_contagion_simulation.pdf",
    ymin=0.42,
    ymax=0.6,
)
