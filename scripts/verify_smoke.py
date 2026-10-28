"""Validate serial/MPI conservation and Exodus I/O for the environment test."""
from pathlib import Path
import json
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import netCDF4
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parents[1]
out = root / "results" / "environment_smoke"


def read_case(stem):
    frame = pd.read_csv(out / f"{stem}.csv")
    if len(frame) != 11 or not np.all(np.isfinite(frame[["time", "mean_c"]])):
        raise RuntimeError(f"{stem}: expected 11 finite output records")
    if not np.isclose(frame.time.iloc[-1], 0.5):
        raise RuntimeError(f"{stem}: simulation did not reach its final time")
    drift = float(np.max(np.abs(frame.mean_c - frame.mean_c.iloc[0])))
    if drift > 1e-8:
        raise RuntimeError(f"{stem}: mean concentration drift {drift} exceeds 1e-8")
    with netCDF4.Dataset(out / f"{stem}.e") as ds:
        names = [str(n).strip() for n in netCDF4.chartostring(ds["name_nod_var"][:])]
        c_index = names.index("c") + 1
        c = np.asarray(ds[f"vals_nod_var{c_index}"][:])
        x = np.asarray(ds["coordx"][:])
        y = np.asarray(ds["coordy"][:])
        order = np.lexsort((x, y))
        c, x, y = c[:, order], x[order], y[order]
    if not np.all(np.isfinite(c)) or np.max(np.abs(c[-1] - c[0])) < 1e-8:
        raise RuntimeError(f"{stem}: concentration field is invalid or did not evolve")
    return frame, c, x, y, drift


serial, cs, xs, ys, ds = read_case("serial")
parallel, cp, xp, yp, dp = read_case("mpi2")
if not (np.allclose(xs, xp) and np.allclose(ys, yp)):
    raise RuntimeError("Serial and MPI meshes differ")
field_difference = float(np.max(np.abs(cs - cp)))
if field_difference > 1e-8:
    raise RuntimeError(f"Serial/MPI field discrepancy: {field_difference}")

kks = pd.read_csv(out / "kks.csv")
if len(kks) != 4 or not np.isclose(kks.time.iloc[-1], 0.3):
    raise RuntimeError("KKS test did not reach the end of all 3 steps")
if not np.all(np.isfinite(kks[["time", "mean_c"]])):
    raise RuntimeError("KKS output contains nonfinite values")
kks_drift = float(np.max(np.abs(kks.mean_c - kks.mean_c.iloc[0])))
if kks_drift > 1e-8:
    raise RuntimeError(f"KKS mean concentration drift: {kks_drift}")
with netCDF4.Dataset(out / "kks.e") as kds:
    knames = [str(n).strip() for n in netCDF4.chartostring(kds["name_nod_var"][:])]
    for required in ["c", "eta", "w", "cm", "cd"]:
        values = np.asarray(kds[f"vals_nod_var{knames.index(required) + 1}"][:])
        if not np.all(np.isfinite(values)):
            raise RuntimeError(f"KKS: nonfinite field {required}")

fig, axes = plt.subplots(1, 3, figsize=(12, 3.6), layout="constrained")
nx, ny = len(np.unique(xs)), len(np.unique(ys))
for ax, field, title in zip(axes[:2], [cs[0], cs[-1]], ["Initial field", "After 10 time steps"]):
    mesh = ax.pcolormesh(xs.reshape(ny, nx), ys.reshape(ny, nx), field.reshape(ny, nx),
                         shading="auto", vmin=0.1, vmax=0.9, cmap="viridis")
    ax.set(aspect="equal", title=title, xlabel="x (dimensionless)", ylabel="y (dimensionless)")
fig.colorbar(mesh, ax=list(axes[:2]), label="c")
axes[2].plot(serial.time, serial.mean_c - serial.mean_c.iloc[0], "o-", label="Serial")
axes[2].plot(parallel.time, parallel.mean_c - parallel.mean_c.iloc[0], "x--", label="2 MPI ranks")
axes[2].set(title="Conservation check", xlabel="Time (dimensionless)", ylabel="Mean c - initial mean c")
axes[2].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
axes[2].legend()
limit = max(1e-12, 1.2 * max(ds, dp))
axes[2].set_ylim(-limit, limit)
fig.suptitle("MOOSE environment verification - generic model, not Ni-Al thermodynamics")
fig.savefig(out / "verification.png", dpi=160)
plt.close(fig)

report = {
    "verified_at_utc": datetime.now(timezone.utc).isoformat(),
    "status": "PASS",
    "test": "Dimensionless Cahn-Hilliard environment smoke test",
    "steps": 10,
    "mesh_elements": 576,
    "serial_mean_concentration_drift": ds,
    "mpi2_mean_concentration_drift": dp,
    "maximum_serial_mpi_nodal_difference": field_difference,
    "absolute_tolerance": 1e-8,
    "exodus_readback": "PASS",
    "kks_mpi2": "PASS",
    "kks_steps": 3,
    "kks_mean_concentration_drift": kks_drift,
    "conservation_note": "Drifts computed from MOOSE CSV output (14 significant digits).",
    "note": "This does not validate a Ni-Al thermodynamic or kinetic model."
}
(out / "verification.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
print(json.dumps(report, indent=2))
