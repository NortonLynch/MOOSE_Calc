"""Convert vendor TABULATE output, validate thermodynamics, and plot results.

Run with the existing WSL moose-post Python. No Thermo-Calc API is required.
Raw vendor files are never modified. All compositions are atomic/mole fractions.
"""
from pathlib import Path
import hashlib
import json
import re

import numpy as np
import pandas as pd
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import root
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
RAW = BASE / "raw"
OUT = BASE / "ni_al_1000K"
OUT.mkdir(exist_ok=True)
INTERCEPT = -47425.823
SLOPE = -135156.45


def table(name, columns):
    """Retain phase-region labels; merge duplicate endpoints only if consistent."""
    rows, phases, in_header = [], [], False
    for line in (RAW / name).read_text().splitlines():
        s = line.strip()
        if s.startswith("Phase Region for:"):
            phases, in_header = [], True
        elif s.startswith("col-1="):
            in_header = False
        elif re.match(r"^[+-]?\d\.\d+E[+-]\d+", s):
            values = [float(x) for x in s.split()]
            assert len(values) == len(columns), (name, line)
            rows.append(dict(zip(columns, values), phase_region=" + ".join(phases)))
        elif in_header and s:
            phases.append(s)
    df = pd.DataFrame(rows)
    assert len(df) > 0, name
    for _, group in df.groupby("x_Al"):
        if len(group) > 1:
            # Near phase-boundary duplicate rows may differ in absent-phase values.
            if name.startswith("gamma"):
                assert np.allclose(group[columns], group[columns].iloc[0], atol=1e-9, rtol=1e-5)
    return df.drop_duplicates("x_Al").sort_values("x_Al").reset_index(drop=True)


def write_csv(df, name):
    df.to_csv(OUT / name, index=False, float_format="%.10g")


bulk_cols = ["x_Al", "Gm_raw_J_mol", "mu_Al_raw_J_mol", "mu_Ni_raw_J_mol", "Vm_m3_mol"]
reduced_cols = ["x_Al", "dGm_dx_raw_J_mol", "Gm_reduced_J_mol", "dGm_reduced_dx_J_mol"]
checks = {}
single = {}
splines = {}
for name, reduced_file, limits in [
    ("gamma", "gamma_shifted_1000K.txt", (0.05, 0.30)),
    ("gamma_prime", "gamma_prime_reduced_1000K.txt", (0.18, 0.30)),
]:
    df = table(f"{name}_bulk_1000K.txt", bulk_cols)
    red = table(reduced_file, reduced_cols).drop(columns="phase_region")
    df = df.merge(red, on="x_Al", validate="one_to_one")
    df["Gm_J_mol"] = df.Gm_reduced_J_mol + INTERCEPT + SLOPE * df.x_Al
    df["dGm_dx_J_mol"] = df.dGm_reduced_dx_J_mol + SLOPE
    # Euler-reconstructed potentials avoid cancellation in the six-digit raw export.
    df["mu_Ni_J_mol"] = df.Gm_J_mol - df.x_Al * df.dGm_dx_J_mol
    df["mu_Al_J_mol"] = df.mu_Ni_J_mol + df.dGm_dx_J_mol
    if name == "gamma_prime":
        order = table("gamma_prime_order_1000K.txt", ["x_Al", "phase_mole_fraction", "y_Al_sublattice_1", "y_Al_sublattice_2"])
        df = df.merge(order.drop(columns="phase_region"), on="x_Al", validate="one_to_one")
        df["order_indicator_y2_minus_y1"] = df.y_Al_sublattice_2 - df.y_Al_sublattice_1
        assert np.allclose(df.phase_mole_fraction, 1)
        assert (df.order_indicator_y2_minus_y1 > 0.6).all(), "Ordered branch lost"
    assert np.isfinite(df.select_dtypes("number")).all().all()
    assert np.allclose([df.x_Al.min(), df.x_Al.max()], limits, atol=1e-8)
    assert df.x_Al.diff().dropna().gt(0).all()
    assert df.x_Al.diff().max() <= 0.000501
    assert (df.Vm_m3_mol > 0).all()
    euler_error = float(np.max(np.abs(df.Gm_raw_J_mol - (df.x_Al*df.mu_Al_raw_J_mol + (1-df.x_Al)*df.mu_Ni_raw_J_mol))))
    g_error = float(np.max(np.abs(df.Gm_raw_J_mol - df.Gm_J_mol)))
    assert euler_error < 0.5, (name, euler_error)
    assert g_error < 0.1, (name, g_error)
    checks[name] = {"rows": len(df), "x_Al_range": list(limits), "max_euler_residual_raw_J_mol": euler_error, "max_G_raw_vs_reconstructed_J_mol": g_error}
    write_csv(df, f"{name}_1000K.csv")
    single[name] = df
    splines[name] = CubicHermiteSpline(df.x_Al, df.Gm_reduced_J_mol, df.dGm_reduced_dx_J_mol, extrapolate=False)

eq = table("equilibrium_bulk_1000K.txt", bulk_cols)
ph = table("equilibrium_phases_1000K.txt", ["x_Al", "mole_fraction_fcc1", "x_Al_fcc1", "mole_fraction_fcc2", "x_Al_fcc2"])
vol = table("equilibrium_volumes_1000K.txt", ["x_Al", "Vm_fcc1_m3_mol", "Vm_fcc2_m3_mol", "volume_fraction_fcc1", "volume_fraction_fcc2"])
order = table("equilibrium_order_1000K.txt", ["x_Al", "y_Al_fcc1_sub1", "y_Al_fcc1_sub2", "y_Al_fcc2_sub1", "y_Al_fcc2_sub2"])
for extra in [ph, vol, order]:
    eq = eq.merge(extra.drop(columns="phase_region"), on="x_Al", validate="one_to_one")
eq["other_phase_mole_fraction"] = (1-eq.mole_fraction_fcc1-eq.mole_fraction_fcc2).clip(lower=0)
eq["other_phase_volume_fraction"] = (1-eq.volume_fraction_fcc1-eq.volume_fraction_fcc2).clip(lower=0)
for i in (1, 2):
    absent = eq[f"mole_fraction_fcc{i}"] < 1e-8
    # The solver retains arbitrary starting constitutions for absent phases.
    for col in [f"x_Al_fcc{i}", f"Vm_fcc{i}_m3_mol", f"y_Al_fcc{i}_sub1", f"y_Al_fcc{i}_sub2"]:
        eq.loc[absent, col] = np.nan
write_csv(eq, "equilibrium_composition_scan_1000K.csv")
target_rows = eq.loc[np.isclose(eq.x_Al, 0.14, atol=1e-10, rtol=0)]
assert len(target_rows) == 1
t = target_rows.iloc[0]
assert t.phase_region == "FCC_L12#1 + FCC_L12#2"
mass_error = float(abs(t.mole_fraction_fcc1*t.x_Al_fcc1+t.mole_fraction_fcc2*t.x_Al_fcc2-0.14))
vol_error = float(abs(t.Vm_m3_mol-(t.mole_fraction_fcc1*t.Vm_fcc1_m3_mol+t.mole_fraction_fcc2*t.Vm_fcc2_m3_mol)))
assert mass_error < 1e-6
assert vol_error < 1e-10
sg, sp = splines["gamma"], splines["gamma_prime"]


def common_tangent(c):
    a, b = c
    da, db = sg(a, 1), sp(b, 1)
    return [(da-db)/1e5, (sg(a)-a*da-sp(b)+b*db)/1e4]


ct = root(common_tangent, [t.x_Al_fcc1, t.x_Al_fcc2], tol=1e-10)
assert ct.success, ct.message
tie_error = float(np.max(np.abs(ct.x-[t.x_Al_fcc1, t.x_Al_fcc2])))
assert tie_error < 2e-6, tie_error
checks["equilibrium"] = {"rows": len(eq), "target_mass_balance_abs_error": mass_error, "target_molar_volume_balance_abs_error_m3_mol": vol_error, "common_tangent_x_Al_from_interpolation": ct.x.tolist(), "max_common_tangent_vs_TC_x_Al_error": tie_error}

phases = pd.DataFrame([
    {"phase": "gamma", "thermo_calc_phase": "FCC_L12#1 (DISORD)", "x_Al": t.x_Al_fcc1, "Al_at_percent": 100*t.x_Al_fcc1, "mole_fraction": t.mole_fraction_fcc1, "volume_fraction": t.volume_fraction_fcc1, "Vm_m3_mol": t.Vm_fcc1_m3_mol},
    {"phase": "gamma_prime", "thermo_calc_phase": "FCC_L12#2 (ORD)", "x_Al": t.x_Al_fcc2, "Al_at_percent": 100*t.x_Al_fcc2, "mole_fraction": t.mole_fraction_fcc2, "volume_fraction": t.volume_fraction_fcc2, "Vm_m3_mol": t.Vm_fcc2_m3_mol},
])
write_csv(phases, "equilibrium_target_1000K_xAl014.csv")
metadata = {
    "software": "Thermo-Calc Educational 2026b", "database": "NIDEMO", "database_version": "2.0", "export_date": "2026-10-19", "method": "Licensed Console Mode, POLY TABULATE; not TC-Python",
    "temperature_K": 1000, "pressure_Pa": 101325, "bulk_x_Al": 0.14,
    "composition_basis": "mole fraction of atoms; Ni + Al = 1",
    "chemical_potential_reference": "SER (vendor default)",
    "affine_reference_J_mol": {"intercept": INTERCEPT, "slope_per_x_Al": SLOPE, "definition": "G_reduced = G - (intercept + slope*x_Al); same reference for both phases"},
    "raw_numeric_precision": "TABULATE: normally six significant digits. Reconstruct G and diffusion potential from separately exported reduced values to preserve small energy differences. This is not extra database accuracy.",
    "gamma_branch": "Only DIS_FCC_A1 retrieved, disordered homogeneous fcc.",
    "gamma_prime_branch": "Only FCC_L12#2 entered, all other phases/composition sets suspended; ordinary minimization; site fractions relaxed and order checked.",
    "scope": "Bulk chemical thermodynamics only; no elastic coherency, curvature, interface energy or diffusion mobility. Single-phase branches may be metastable against other phases.",
    "target_phases": phases.to_dict("records"),
    "validation": checks,
    "input_sha256": {str(p.relative_to(BASE)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(RAW.glob("*.txt"))},
}
(OUT / "metadata.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
(OUT / "validation.json").write_text(json.dumps({"passed": True, **checks}, indent=2)+"\n")

plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
fig, axs = plt.subplots(1, 3, figsize=(14, 4.4), constrained_layout=True)
colors = {"gamma": "#2166ac", "gamma_prime": "#d6604d"}
labels = {"gamma": r"$\gamma$ (disordered)", "gamma_prime": r"$\gamma^\prime$ (ordered)"}
for name, df in single.items():
    axs[0].plot(100*df.x_Al, df.Gm_reduced_J_mol/1000, color=colors[name], label=labels[name], lw=2)
    axs[1].plot(100*df.x_Al, df.dGm_reduced_dx_J_mol/1000, color=colors[name], lw=2)
for ax in axs[:2]:
    ax.axhline(0, color="#777", ls="--", lw=1)
    ax.set_xlim(7, 28)
    ax.set_xlabel("Al (at.%)")
    ax.grid(alpha=0.15)
axs[0].set_ylim(-0.05, 3)
axs[0].set_ylabel("Reduced molar Gibbs energy (kJ/mol)")
axs[0].set_title("Homogeneous free-energy branches")
axs[0].legend(frameon=False, fontsize=10)
axs[1].set_ylim(-30, 50)
axs[1].set_ylabel(r"$\mu_{Al}-\mu_{Ni}-$ reference (kJ/mol)")
axs[1].set_title("Diffusion potential")
axs[0].scatter(100*ct.x, [sg(ct.x[0])/1000, sp(ct.x[1])/1000], color=[colors["gamma"], colors["gamma_prime"]], zorder=4)
v = phases.volume_fraction.to_numpy()*100
bars = axs[2].bar([r"$\gamma$", r"$\gamma^\prime$"], v, color=list(colors.values()), width=0.55)
axs[2].bar_label(bars, labels=[f"{x:.2f}%" for x in v], padding=4)
axs[2].set_ylim(0, 100)
axs[2].set_ylabel("Equilibrium volume fraction (%)")
axs[2].set_title("Ni–14 at.% Al at 1000 K")
axs[2].set_xlabel("Phase")
fig.suptitle("Thermo-Calc Educational 2026b / NIDEMO 2.0 / 1000 K / 101325 Pa", fontsize=13)
fig.savefig(OUT / "thermodynamics_overview.png", dpi=180)
fig.savefig(OUT / "thermodynamics_overview.pdf")
plt.close(fig)
print(json.dumps({"output": str(OUT), "validation": checks}, indent=2))
