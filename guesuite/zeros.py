'''
Zero finding utilities for real-valued signals.
'''
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, List, Tuple, Optional

import mpmath as mp


@dataclass
class ZeroScanConfig:
    t_min: float = 10.0
    t_max: float = 200.0
    dt: float = 0.02
    refine_max_iter: int = 80
    refine_tol: float = 1e-10
    max_zeros: Optional[int] = 200


def bracket_zeros(f: Callable[[float], float], cfg: ZeroScanConfig) -> List[Tuple[float, float]]:
    t = cfg.t_min
    fa = mp.mpf(f(t))
    brackets: List[Tuple[float, float]] = []
    while t < cfg.t_max:
        t2 = min(cfg.t_max, t + cfg.dt)
        fb = mp.mpf(f(t2))
        if fa == 0:
            brackets.append((t, t))
        elif fb == 0:
            brackets.append((t2, t2))
        elif fa * fb < 0:
            brackets.append((t, t2))
        t, fa = t2, fb
        if cfg.max_zeros is not None and len(brackets) >= cfg.max_zeros:
            break
    return brackets


def refine_root_bisect(f: Callable[[float], float], a: float, b: float, cfg: ZeroScanConfig) -> float:
    if a == b:
        return a
    fa = mp.mpf(f(a))
    fb = mp.mpf(f(b))
    if fa == 0:
        return a
    if fb == 0:
        return b
    if fa * fb > 0:
        raise ValueError("Interval does not bracket a root.")

    lo, hi = mp.mpf(a), mp.mpf(b)
    for _ in range(cfg.refine_max_iter):
        mid = (lo + hi) / 2
        fm = mp.mpf(f(mid))
        if abs(fm) < cfg.refine_tol or abs(hi - lo) < cfg.refine_tol:
            return float(mid)
        if fa * fm < 0:
            hi, fb = mid, fm
        else:
            lo, fa = mid, fm
    return float((lo + hi) / 2)


def find_zeros(f: Callable[[float], float], cfg: ZeroScanConfig) -> List[float]:
    brackets = bracket_zeros(f, cfg)
    zeros: List[float] = []
    for a, b in brackets:
        try:
            z = refine_root_bisect(f, a, b, cfg)
            if not zeros or abs(z - zeros[-1]) > max(cfg.dt, 1e-6):
                zeros.append(z)
        except Exception:
            continue
        if cfg.max_zeros is not None and len(zeros) >= cfg.max_zeros:
            break
    return zeros
