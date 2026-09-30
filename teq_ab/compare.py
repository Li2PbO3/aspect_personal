#!/usr/bin/env python3
"""Compare the two A/B runs in teq_ab/ (fluid pressure vs adiabatic/lithostatic
pressure as the pressure of the thermodynamic equilibrium calculation).

Usage:  python3 teq_ab/compare.py [step]
        step defaults to the last available output step.
"""
import glob
import os
import sys

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

FIELDS = ("p_f", "p_c", "porosity", "T", "density", "viscosity",
          "permeability", "melting_rate", "debug_1", "debug_2", "debug_3")


def load(directory, step):
    acc = {}
    for f in sorted(glob.glob(f"{directory}/solution/{step}.0*.vtu")):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(f)
        reader.Update()
        grid = reader.GetOutput()
        point_data = grid.GetPointData()
        for name in FIELDS:
            arr = point_data.GetArray(name)
            if arr is not None:
                acc.setdefault(name, []).append(vtk_to_numpy(arr))
        pts = vtk_to_numpy(grid.GetPoints().GetData())
        acc.setdefault("x", []).append(pts[:, 0])
        acc.setdefault("y", []).append(pts[:, 1])
    return {k: np.concatenate(v).astype(np.float64) for k, v in acc.items()}


def main():
    step = sys.argv[1] if len(sys.argv) > 1 else None
    steps = sorted({os.path.basename(f).split(".")[0]
                    for f in glob.glob("teq_ab/out_fluid/solution/solution-*.vtu")})
    if not steps:
        sys.exit("no output found in teq_ab/out_fluid; run the two decks first")
    step = step or steps[-1]
    print(f"comparing step {step}\n")

    a = load("teq_ab/out_fluid", step)
    b = load("teq_ab/out_adiabatic", step)

    print(f"{'field':14s} {'max|diff|':>12s} {'rms diff':>12s} {'max|value|':>12s} "
          f"{'rel (rms/scale)':>16s}")
    for name in FIELDS:
        if name not in a or name not in b:
            continue
        d = b[name] - a[name]
        scale = max(np.abs(a[name]).max(), 1e-300)
        print(f"{name:14s} {np.abs(d).max():12.4e} {np.sqrt((d**2).mean()):12.4e} "
              f"{scale:12.4e} {np.sqrt((d**2).mean())/scale:16.3e}")

    # vertical profile along the axis of the column (x = 0)
    order = np.argsort(a["y"], kind="stable")
    x = a["x"][order]
    sel = order[np.abs(x) < 1e-9]
    depth = a["y"].max() - a["y"][sel]
    print(f"\n{'depth [m]':>10s} {'phi fluid':>13s} {'phi adiab':>13s} "
          f"{'dphi':>12s} {'p_c [Pa]':>12s} {'T [C]':>9s}")
    for k in range(0, depth.size, max(1, depth.size // 25)):
        i = sel[k]
        print(f"{depth[k]:10.0f} {a['porosity'][i]:13.5e} {b['porosity'][i]:13.5e} "
              f"{b['porosity'][i]-a['porosity'][i]:12.3e} {a['p_c'][i]:12.1f} "
              f"{a['T'][i]-273.15:9.2f}")

    print(f"\nporosity: max={a['porosity'].max():.4e} "
          f"rel_rms_diff={np.sqrt(((b['porosity']-a['porosity'])**2).mean())/np.sqrt((a['porosity']**2).mean()):.3e}")
    print(f"p_c: max|.|={np.abs(a['p_c']).max():.4e} Pa  "
          f"({100*np.abs(a['p_c']).max()/a['p_f'].mean():.3e} % of p_f)")


if __name__ == "__main__":
    main()
