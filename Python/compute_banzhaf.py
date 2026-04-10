#!/usr/bin/env python3
"""
compute_banzhaf.py — Post-hoc Banzhaf values from GP-SEIR ablation results.

Reads per-coalition validation metrics produced by the R analysis pipeline
and computes exact Shapley and Banzhaf values for all four data sources
across all available characteristic functions.

Usage
-----
    python compute_banzhaf.py                          # process FALSE_ and TRUE_ defaults
    python compute_banzhaf.py input.csv output.csv     # custom paths

Output columns
--------------
metric, metric_name, source, shapley, banzhaf, norm_banzhaf, banzhaf_stable,
singleton, v_grand

norm_banzhaf = "NA" when the normalisation is unreliable because either
  |sum(beta_i)| < 0.05, or sign(sum) != sign(v(N)).

Dependencies: Python standard library only (no external packages).
"""

from __future__ import annotations
import csv, sys
from itertools import combinations
from math import factorial

PLAYERS = ["cases", "hosp", "icu", "radar"]
N = len(PLAYERS)

# All characteristic functions produced by gp_seir_complete_v2.R
METRICS = [
    ("rt_corr",       "Full-period Rt Pearson"),
    ("rt_spearman",   "Full-period Rt Spearman"),
    ("rt_skill",      "Full-period Rt skill score"),
    ("prev_corr",     "Full-period Prev Pearson"),
    ("prev_spearman", "Full-period Prev Spearman"),
    ("prev_skill",    "Full-period Prev skill score"),
    ("rt_corr_pre",   "Pre-June-12 Rt Pearson"),
    ("rt_corr_post",  "Post-June-12 Rt Pearson"),
    ("rt_wave1",      "Wave-1 Rt Pearson"),
    ("rt_inter",      "Inter-wave Rt Pearson"),
    ("rt_wave2",      "Wave-2 Rt Pearson"),
    ("prev_wave1",    "Wave-1 Prev Pearson"),
    ("prev_inter",    "Inter-wave Prev Pearson"),
    ("prev_wave2",    "Wave-2 Prev Pearson"),
    ("rt_sk_w1",      "Wave-1 Rt skill score"),
    ("rt_sk_w2",      "Wave-2 Rt skill score"),
]

BANZHAF_STABLE_THRESHOLD = 0.05


def read_ablation(path: str) -> dict[str, dict[str, str]]:
    """Read ablation CSV keyed by sorted source-set label."""
    with open(path, newline="") as f:
        return {row["sources"]: row for row in csv.DictReader(f)}


def get_v(abl: dict, src_set: list[str], col: str) -> float | None:
    """Return v(S) for coalition src_set; v(empty) = 0 by convention."""
    if not src_set:
        return 0.0
    key = "+".join(sorted(src_set))
    val = abl.get(key, {}).get(col, "")
    if val in (None, "", "NA"):
        return None
    try:
        return float(val)
    except ValueError:
        return None


def shapley(abl: dict, col: str) -> dict[str, float]:
    """Exact Shapley values for n=4 players."""
    phi: dict[str, float] = {}
    for i in PLAYERS:
        others = [p for p in PLAYERS if p != i]
        total = 0.0
        for s_size in range(N):
            coals = [()] if s_size == 0 else list(combinations(others, s_size))
            w = factorial(s_size) * factorial(N - s_size - 1) / factorial(N)
            for S in coals:
                S = list(S)
                vw, vnw = get_v(abl, S + [i], col), get_v(abl, S, col)
                if vw is not None and vnw is not None:
                    total += w * (vw - vnw)
        phi[i] = total
    return phi


def banzhaf(abl: dict, col: str) -> dict[str, float]:
    """Banzhaf values: equal weight 1/2^(n-1) on every coalition."""
    norm = 1.0 / (2 ** (N - 1))
    beta: dict[str, float] = {}
    for i in PLAYERS:
        others = [p for p in PLAYERS if p != i]
        total = 0.0
        for s_size in range(N):
            coals = [()] if s_size == 0 else list(combinations(others, s_size))
            for S in coals:
                S = list(S)
                vw, vnw = get_v(abl, S + [i], col), get_v(abl, S, col)
                if vw is not None and vnw is not None:
                    total += vw - vnw
        beta[i] = norm * total
    return beta


def normalise_banzhaf(
    beta_dict: dict[str, float], vN: float
) -> tuple[dict[str, float | None], bool]:
    """
    Rescale Banzhaf so sum = v(N).

    Returns (values_dict, is_stable).
    is_stable = False when |sum(beta)| < threshold or sign mismatch.
    """
    s = sum(beta_dict.values())
    sign_mismatch = abs(s) > 1e-10 and abs(vN) > 0.01 and (s * vN < 0)
    if abs(s) < BANZHAF_STABLE_THRESHOLD or sign_mismatch:
        return {p: None for p in PLAYERS}, False
    return {p: beta_dict[p] * vN / s for p in PLAYERS}, True


def process(label: str, in_path: str, out_path: str) -> None:
    print(f"\n{'=' * 72}\n{label}\nInput : {in_path}\nOutput: {out_path}\n{'=' * 72}")
    try:
        abl = read_ablation(in_path)
    except FileNotFoundError:
        print(f"  File not found — skipping.")
        return

    rows: list[dict] = []
    for col, name in METRICS:
        vN = get_v(abl, PLAYERS, col)
        if vN is None:
            print(f"  {name}: v(N) missing — skipping")
            continue
        phi  = shapley(abl, col)
        beta = banzhaf(abl, col)
        nb, stable = normalise_banzhaf(beta, vN)

        phi_sum  = sum(phi.values())
        beta_sum = sum(beta.values())
        eff_err  = abs(phi_sum - vN)
        flag = "  *** EFFICIENCY VIOLATED ***" if eff_err > 1e-3 else ""
        print(
            f"  {name:<40s}  v(N)={vN:+.4f}  "
            f"sum(phi)={phi_sum:+.4f} (err={eff_err:.1e})  "
            f"sum(beta)={beta_sum:+.4f}  stable={stable}{flag}"
        )

        for p in PLAYERS:
            singleton = get_v(abl, [p], col)
            rows.append({
                "metric":         col,
                "metric_name":    name,
                "source":         p,
                "shapley":        f"{phi[p]:.8f}",
                "banzhaf":        f"{beta[p]:.8f}",
                "norm_banzhaf":   f"{nb[p]:.8f}" if nb[p] is not None else "NA",
                "banzhaf_stable": str(stable).upper(),
                "singleton":      f"{singleton:.8f}" if singleton is not None else "NA",
                "v_grand":        f"{vN:.8f}",
            })

    if rows:
        with open(out_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"\n  Saved {len(rows)} rows to: {out_path}")
    else:
        print("  No rows produced.")


def main() -> None:
    if len(sys.argv) >= 3:
        process("Custom", sys.argv[1], sys.argv[2])
        return
    # Default: both configurations
    process(
        "FALSE — no delay kernels (primary results)",
        "FALSE_ablation_results.csv",
        "FALSE_banzhaf_values.csv",
    )
    process(
        "TRUE  — with log-normal delay kernels",
        "TRUE_ablation_results.csv",
        "TRUE_banzhaf_values.csv",
    )


if __name__ == "__main__":
    main()
