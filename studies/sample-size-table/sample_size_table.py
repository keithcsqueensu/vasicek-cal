# SPDX-License-Identifier: Apache-2.0
"""S-2 sample-size planning table (studies/sample-size-table/PREDICTION.md).

Reads the committed recovery summary (tests/golden/recovery/summary.csv; no new fits) and, for
each cell (PD, rho, n) and each parameter, fits RMSE ~ C / sqrt(T) across T in {20, 40, 100},
checks the scaling, and solves for the T at which 1.96 x RMSE meets the target: +-0.05 absolute
for rho, +-25% relative for PD. Standard library only.

    python studies/sample-size-table/sample_size_table.py > table.md
"""
import csv
import math
import pathlib
from statistics import NormalDist

ROOT = pathlib.Path(__file__).resolve().parents[2]
SUMMARY = ROOT / "tests" / "golden" / "recovery" / "summary.csv"

Z = 1.959963984540054          # z_0.975
RHO_TARGET = 0.05              # absolute
PD_TARGET_REL = 0.25           # relative
EDGE_MAJORITY = 0.5            # "mostly at a bound": more than half the replicates on the grid edge
SLOPE_RANGE = (-0.75, -0.30)   # the scaling check (step 2)
RESIDUAL_MAX = 0.15
T_STUDIED = (20, 100)
T_FIT = (20, 40, 100)           # the fitted points, exactly (PREDICTION.md, step 1)


def load():
    with open(SUMMARY, newline="") as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    cells = {}
    for r in rows:
        if int(r["periods"]) not in T_FIT:
            continue
        key = (float(r["pd"]), float(r["rho"]), int(r["obligors"]))
        cells.setdefault(key, []).append({
            "T": int(r["periods"]),
            "edge": int(r["edge"]) / int(r["replicates"]),
            "pd": float.fromhex(r["pd_rmse_hex"]),
            "rho": float.fromhex(r["rho_rmse_hex"]),
            "scenario": int(r["scenario"]),
        })
    for key, v in cells.items():
        v.sort(key=lambda p: p["T"])
        if tuple(p["T"] for p in v) != T_FIT:
            raise SystemExit(f"cell {key}: the summary holds T = {[p['T'] for p in v]}, not exactly {list(T_FIT)}")
    return cells


def ols(xs, ys):
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    a = my - b * mx
    return a, b, max(abs(y - (a + b * x)) for x, y in zip(xs, ys))


def analyse(points, param, target):
    used = [p for p in points if p["edge"] <= EDGE_MAJORITY]
    out = {"used": len(used)}
    if not used:
        out["status"] = "not estimable"
        return out
    logs = [(math.log(p["T"]), math.log(p[param])) for p in used]
    log_c = sum(y + 0.5 * x for x, y in logs) / len(logs)
    out["C"] = math.exp(log_c)
    out["T_fixed"] = (Z * out["C"] / target) ** 2
    if len(used) == 3:
        a, b, res = ols([x for x, _ in logs], [y for _, y in logs])
        out["slope"], out["residual"] = b, res
        out["scaling"] = "holds" if SLOPE_RANGE[0] <= b <= SLOPE_RANGE[1] and res <= RESIDUAL_MAX else "fails"
        if out["scaling"] == "holds":
            # RMSE(T) = exp(a) T^b = Z-scaled target  =>  T = (target / (Z exp(a)))^(1/b)
            out["T_slope"] = (target / (Z * math.exp(a))) ** (1.0 / b)
    else:
        out["scaling"] = "unchecked"
    out["status"] = "ok"
    return out


def fmt_t(t):
    if t is None:
        return "-"
    whole = math.ceil(t - 1e-9)
    mark = "" if T_STUDIED[0] <= whole <= T_STUDIED[1] else " (extrap.)"
    return f"{whole:,}{mark}"


def benchmark(pd, rho):
    nd = NormalDist()
    c = nd.inv_cdf(pd)
    lam = nd.pdf(c) / nd.cdf(c)
    return math.sqrt(2.0) * rho * (1.0 - rho), lam * math.sqrt(rho)


def main():
    cells = load()
    print("| PD | rho | n | n PD | param | points used | slope | max resid | scaling | C | large-n C | T* (slope -1/2) | T* (fitted slope) |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    results = {}
    for key in sorted(cells):
        pd, rho, n = key
        bench_rho, bench_pd = benchmark(pd, rho)
        for param, target, bench in (("rho", RHO_TARGET, bench_rho), ("pd", PD_TARGET_REL * pd, bench_pd * pd)):
            r = analyse(cells[key], param, target)
            results[(key, param)] = r
            if r["status"] == "not estimable":
                print(f"| {pd:g} | {rho:g} | {n:,} | {n * pd:g} | {param} | 0 | | | | | | not estimable at any T studied | |")
                continue
            slope = f"{r['slope']:.3f}" if "slope" in r else ""
            resid = f"{r['residual']:.3f}" if "residual" in r else ""
            c_str = f"{r['C'] / pd:.4g} x PD" if param == "pd" else f"{r['C']:.4g}"
            b_str = f"{bench / pd:.4g} x PD" if param == "pd" else f"{bench:.4g}"
            print(f"| {pd:g} | {rho:g} | {n:,} | {n * pd:g} | {param} | {r['used']} | {slope} | {resid} | {r['scaling']} "
                  f"| {c_str} | {b_str} | {fmt_t(r['T_fixed'])} | {fmt_t(r.get('T_slope'))} |")
    # The planning grids: fixed-slope T* per (PD, n) row and rho column, one grid per target.
    pds = sorted({k[0] for k in cells})
    rhos = sorted({k[1] for k in cells})
    ns = sorted({k[2] for k in cells})
    for param, title in (("rho", "Years for rho within +-0.05 (95%)"), ("pd", "Years for PD within +-25% (95%)")):
        print(f"\n{title}\n")
        print("| PD | n | n PD | " + " | ".join(f"rho = {r:g}" for r in rhos) + " |")
        print("|---|---|---|" + "---|" * len(rhos))
        for pd in pds:
            for n in ns:
                cols = []
                for rho in rhos:
                    r = results[((pd, rho, n), param)]
                    if r["status"] == "not estimable":
                        cols.append("not estimable")
                    else:
                        cols.append(fmt_t(r["T_fixed"]) + ("" if r["scaling"] == "holds" else f" ({r['scaling']})"))
                print(f"| {pd:g} | {n:,} | {n * pd:g} | " + " | ".join(cols) + " |")
    return results


if __name__ == "__main__":
    main()
