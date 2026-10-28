"""Use the installed, licensed TC-Python entry point without modifying it."""
from pathlib import Path
import json
import traceback
from tc_python import TCPython, ThermodynamicQuantity

out = Path(__file__).resolve().parent / "probe"
out.mkdir(exist_ok=True)
try:
    with TCPython() as session:
        system = session.select_database_and_elements("NIDEMO", ["NI", "AL"]).get_system()
        calc = (system.with_single_equilibrium_calculation()
                .set_condition(ThermodynamicQuantity.temperature(), 1000.0)
                .set_condition(ThermodynamicQuantity.pressure(), 101325.0)
                .set_condition(ThermodynamicQuantity.mole_fraction_of_a_component("AL"), 0.14))
        result = calc.calculate()
        info = {"status": "success", "stable_phases": result.get_stable_phases()}
        print(json.dumps(info, indent=2))
        (out / "tc_python_status.json").write_text(json.dumps(info, indent=2), encoding="utf-8")
except Exception as exc:
    info = {"status": "failed", "error_type": type(exc).__name__, "error": str(exc)}
    (out / "tc_python_status.json").write_text(json.dumps(info, indent=2), encoding="utf-8")
    traceback.print_exc()
    raise SystemExit(1)
