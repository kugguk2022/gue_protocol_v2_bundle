#!/usr/bin/env python3
'''
Run the three-model experiment suite:

  1) Independent primes baseline (partial Euler product)
  2) Riemann-Siegel Z(t) approximation (functional equation symmetry)
  3) Full zeta ground truth (mpmath.zeta) mapped into Z(t)

Also runs two negative controls:
  - jittered "fake primes"
  - phase-randomized RS sum
'''
from __future__ import annotations

import argparse
import json
from pathlib import Path

from guesuite.models import (
    Z_full, Z_riemann_siegel, Z_euler_partial,
    jitter_primes, Z_euler_partial_direct_float_primes, phase_randomized_rs
)
from guesuite.zeros import ZeroScanConfig, find_zeros
from guesuite.stats import unfold_spacings, summarize_spacings


def primes_up_to(n: int) -> list[int]:
    sieve = [True] * (n + 1)
    if n >= 0:
        sieve[0] = False
    if n >= 1:
        sieve[1] = False
    r = int(n ** 0.5)
    for p in range(2, r + 1):
        if sieve[p]:
            for q in range(p * p, n + 1, p):
                sieve[q] = False
    return [i for i in range(2, n + 1) if sieve[i]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--t-min", type=float, default=10.0)
    ap.add_argument("--t-max", type=float, default=200.0)
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--max-zeros", type=int, default=120)

    ap.add_argument("--p-max", type=int, default=2000, help="max prime for partial Euler product baseline")
    ap.add_argument("--k-max", type=int, default=1, help="Euler product power truncation (k>=1)")
    ap.add_argument("--seed", type=int, default=0)

    ap.add_argument("--outdir", type=str, default="out_gue")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    primes = primes_up_to(args.p_max)

    cfg = ZeroScanConfig(t_min=args.t_min, t_max=args.t_max, dt=args.dt, max_zeros=args.max_zeros)

    fake_primes = jitter_primes(primes, width=0.5, seed=args.seed)

    models = {
        "independent_primes": lambda t: float(Z_euler_partial(t, primes=primes, k_max=args.k_max, use_log=False)),
        "riemann_siegel": lambda t: float(Z_riemann_siegel(t)),
        "full_zeta": lambda t: float(Z_full(t)),
        "fake_primes_jitter": lambda t: float(Z_euler_partial_direct_float_primes(t, primes_like=fake_primes)),
        "rs_phase_randomized": lambda t: float(phase_randomized_rs(t, seed=args.seed)),
    }

    report = {}

    for name, f in models.items():
        zeros = find_zeros(f, cfg)
        spacings = unfold_spacings(zeros)
        summary = summarize_spacings(spacings)

        (outdir / f"zeros_{name}.csv").write_text(
            "t\\n" + "\\n".join(f"{z:.12f}" for z in zeros) + "\\n",
            encoding="utf-8"
        )

        report[name] = {
            "scan": {"t_min": args.t_min, "t_max": args.t_max, "dt": args.dt, "max_zeros": args.max_zeros},
            "p_max": args.p_max,
            "k_max": args.k_max,
            "zeros_found": len(zeros),
            "spacing_summary": {
                "n": summary["n"],
                "mean": summary["mean"],
                "var": summary["var"],
                "r_mean": summary["r_mean"],
                "ks_poisson_D": summary["ks_poisson"][0],
                "ks_poisson_p": summary["ks_poisson"][1],
                "ks_gue_D": summary["ks_gue"][0],
                "ks_gue_p": summary["ks_gue"][1],
            },
        }

    (outdir / "summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
