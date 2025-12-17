'''
Statistics: unfolding, spacing distributions, KS tests, r-statistic.
'''
from __future__ import annotations

import math
from typing import Callable, Sequence, Dict, Any, Tuple

import numpy as np


def unfold_spacings(t_zeros: Sequence[float]) -> np.ndarray:
    t = np.asarray(t_zeros, dtype=float)
    if t.size < 2:
        return np.array([], dtype=float)
    dt = np.diff(t)
    density = np.log(t[:-1] / (2 * math.pi)) / (2 * math.pi)
    s = dt * density
    m = float(s.mean()) if s.size else 1.0
    return s / m


def cdf_poisson(s: np.ndarray) -> np.ndarray:
    return 1.0 - np.exp(-s)


def pdf_gue_wigner(s: np.ndarray) -> np.ndarray:
    # GUE Wigner surmise (beta=2):
    # P(s) = (32/pi^2) s^2 exp(-4 s^2 / pi)
    return (32.0 / (math.pi ** 2)) * (s ** 2) * np.exp(-4.0 * (s ** 2) / math.pi)


def cdf_gue_wigner(s: np.ndarray) -> np.ndarray:
    # F(s) = erf(2 s / sqrt(pi)) - (4/pi) s exp(-4 s^2 / pi)
    s = np.asarray(s, dtype=float)
    erf_vec = np.vectorize(lambda x: math.erf(2.0 * x / math.sqrt(math.pi)))
    return erf_vec(s) - (4.0 / math.pi) * s * np.exp(-4.0 * (s ** 2) / math.pi)


def ks_statistic(sample: np.ndarray, cdf: Callable[[np.ndarray], np.ndarray]) -> float:
    x = np.sort(np.asarray(sample, dtype=float))
    n = x.size
    if n == 0:
        return float("nan")
    F = cdf(x)
    i = np.arange(1, n + 1)
    d_plus = np.max(i / n - F)
    d_minus = np.max(F - (i - 1) / n)
    return float(max(d_plus, d_minus))


def ks_pvalue(d: float, n: int) -> float:
    if n <= 0 or not np.isfinite(d):
        return float("nan")
    en = math.sqrt(n)
    lam = (en + 0.12 + 0.11 / en) * d
    s = 0.0
    for k in range(1, 200):
        term = 2.0 * ((-1) ** (k - 1)) * math.exp(-2.0 * (k * k) * (lam * lam))
        s += term
        if abs(term) < 1e-12:
            break
    return max(0.0, min(1.0, s))


def ks_test(sample: np.ndarray, cdf: Callable[[np.ndarray], np.ndarray]) -> Tuple[float, float]:
    sample = np.asarray(sample, dtype=float)
    d = ks_statistic(sample, cdf)
    p = ks_pvalue(d, sample.size)
    return d, p


def r_statistic(spacings: np.ndarray) -> float:
    s = np.asarray(spacings, dtype=float)
    if s.size < 2:
        return float("nan")
    a = s[1:]
    b = s[:-1]
    r = np.minimum(a, b) / np.maximum(a, b)
    return float(np.mean(r))


def summarize_spacings(spacings: np.ndarray) -> Dict[str, Any]:
    s = np.asarray(spacings, dtype=float)
    ks_poi = ks_test(s, cdf_poisson)
    ks_gue = ks_test(s, cdf_gue_wigner)
    return {
        "n": int(s.size),
        "mean": float(np.mean(s)) if s.size else float("nan"),
        "var": float(np.var(s)) if s.size else float("nan"),
        "r_mean": r_statistic(s),
        "ks_poisson": ks_poi,
        "ks_gue": ks_gue,
    }
