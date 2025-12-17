'''
Spectral models for zeta / prime-block experiments.

- Independent Euler-product model (baseline, tends to Poisson unless extra constraints added)
- Riemann-Siegel Z(t) approximation (encodes functional equation symmetry)
- Full zeta ground truth (mpmath.zeta), mapped into Z(t) via theta(t)

Notes:
  * Z(t) is defined so that zeta(1/2+it) = Z(t) * exp(-i*theta(t)), hence Z(t) is (approximately) real.
'''
from __future__ import annotations

import math
from typing import Sequence, Optional, List

import mpmath as mp

mp.mp.dps = 50  # default precision; override in caller if needed


def theta(t: float) -> mp.mpf:
    '''
    Riemann-Siegel theta function:
        theta(t) = Im(log Gamma(1/4 + i t/2)) - (t/2) log(pi)
    '''
    t = mp.mpf(t)
    return mp.im(mp.log(mp.gamma(mp.mpf("0.25") + 0.5j * t))) - (t / 2) * mp.log(mp.pi)


def Z_full(t: float) -> mp.mpf:
    '''
    Ground truth Z(t) computed from mpmath.zeta using:
        Z(t) = zeta(1/2 + i t) * exp(i theta(t))
    '''
    t = mp.mpf(t)
    return mp.re(mp.zeta(mp.mpf("0.5") + 1j * t) * mp.e ** (1j * theta(t)))


def Z_riemann_siegel(t: float, n_terms: Optional[int] = None) -> mp.mpf:
    '''
    Basic Riemann-Siegel Z(t) approximation:
        Z(t) approx 2 * sum_{n=1}^{N} cos(t log n - theta(t)) / sqrt(n)
    where N = floor(sqrt(t/(2*pi))) if n_terms is not provided.

    This is not the full RS formula with remainder terms; it's enough to test whether
    "global mirror symmetry" changes spacing statistics compared to independent prime models.
    '''
    t = mp.mpf(t)
    th = theta(t)
    if n_terms is None:
        if t <= 0:
            return mp.mpf("0.0")
        n_terms = int(mp.floor(mp.sqrt(t / (2 * mp.pi))))
        n_terms = max(1, n_terms)

    s = mp.mpf("0.0")
    for n in range(1, n_terms + 1):
        s += mp.cos(t * mp.log(n) - th) / mp.sqrt(n)
    return 2 * s


def Z_euler_partial(
    t: float,
    primes: Sequence[int],
    k_max: int = 1,
    use_log: bool = True,
) -> mp.mpf:
    '''
    Baseline "independent primes" model from a partial Euler product on the critical line.

    Option A (use_log=True): uses log Euler product amplitude:
        log A(t) approx sum_{p in primes} sum_{k=1..k_max} (1/k) p^{-k/2} cos(k t log p)
    and returns exp(log A(t)) as a crude real-positive signal.

    Option B (use_log=False): returns Re(prod_{p} (1 - p^{-1/2 - i t})^{-1})
    (more oscillatory, slower).

    WARNING: this model does NOT encode the functional equation, and typically fails to show
    GUE level repulsion (Poisson is the expected baseline).
    '''
    t = mp.mpf(t)

    if use_log:
        re_log = mp.mpf("0.0")
        for p in primes:
            lp = mp.log(p)
            for k in range(1, k_max + 1):
                re_log += (mp.cos(k * t * lp) / k) * (mp.power(p, -mp.mpf(k) / 2))
        return mp.e ** re_log

    prod = 1
    for p in primes:
        prod *= 1 / (1 - mp.power(p, -mp.mpf("0.5") - 1j * t))
    return mp.re(prod)


def jitter_primes(primes: Sequence[int], width: float = 0.5, seed: int = 0) -> List[float]:
    '''
    Negative control: pseudo-primes r_k ~ Uniform(p_k - width, p_k + width).
    '''
    import random
    rng = random.Random(seed)
    return [p + rng.uniform(-width, width) for p in primes]


def Z_euler_partial_float_primes(
    t: float,
    primes_like: Sequence[float],
    k_max: int = 1,
) -> mp.mpf:
    '''
    Log Euler product amplitude using non-integer 'primes' for negative controls.
    '''
    t = mp.mpf(t)
    re_log = mp.mpf("0.0")
    for p in primes_like:
        lp = mp.log(p)
        for k in range(1, k_max + 1):
            re_log += (mp.cos(k * t * lp) / k) * (mp.power(p, -mp.mpf(k) / 2))
    return mp.e ** re_log


def phase_randomized_rs(t: float, seed: int = 0, n_terms: Optional[int] = None) -> mp.mpf:
    '''
    Negative control: RS-like sum with randomized phases.
    Preserves amplitude scaling 1/sqrt(n) but destroys coherent structure.
    '''
    import random
    rng = random.Random(seed)
    t = mp.mpf(t)
    if n_terms is None:
        n_terms = int(mp.floor(mp.sqrt(t / (2 * mp.pi))))
        n_terms = max(1, n_terms)

    s = mp.mpf("0.0")
    for n in range(1, n_terms + 1):
        phi = rng.random() * 2 * math.pi
        s += mp.cos(t * mp.log(n) + phi) / mp.sqrt(n)
    return 2 * s

def Z_euler_partial_direct_float_primes(t: float, primes_like: Sequence[float]) -> mp.mpf:
    '''
    Direct Euler product real part using non-integer 'primes' (negative control).
    Returns Re(prod_{p} (1 - p^{-1/2 - i t})^{-1}).
    '''
    t = mp.mpf(t)
    prod = 1
    for p in primes_like:
        prod *= 1 / (1 - mp.power(p, -mp.mpf("0.5") - 1j * t))
    return mp.re(prod)
