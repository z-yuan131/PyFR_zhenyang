from functools import cache
from io import StringIO
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator


_DATA_DIR = Path(__file__).resolve().parent


@cache
def _polar_interpolators():
    with np.load(_DATA_DIR / 'polarCLCD.npz') as mesh:
        re = mesh['Re'].copy()
        alpha = mesh['alpha'].copy()
        data = mesh['data'].copy()

    cl_grid, cd_grid = data[0], data[1]

    return (
        RegularGridInterpolator((re, alpha), cl_grid,
                                bounds_error=False, fill_value=None),
        RegularGridInterpolator((re, alpha), cd_grid,
                                bounds_error=False, fill_value=None)
    )


@cache
def _blade_data():
    data_lines = []
    with open(_DATA_DIR / '13x6-PERF.txt') as f:
        for line in f:
            stripped = line.strip()
            if not stripped or not stripped[0].isdigit():
                continue

            data_lines.append(line)

    data = np.loadtxt(StringIO(''.join(data_lines)))

    return (
        np.asarray(data[:, 0], dtype=float).ravel(),
        np.asarray(data[:, 7]*np.pi/180, dtype=float).ravel(),
        np.asarray(data[:, 1], dtype=float).ravel()
    )


def clcd(alphas, reynolds):
    interp_cl, interp_cd = _polar_interpolators()
    points = np.column_stack((reynolds, alphas))

    return interp_cl(points), interp_cd(points)


def pitch(r):
    radius, pitch_angle, chord = _blade_data()

    return np.interp(r, radius, pitch_angle), np.interp(r, radius, chord)
