# Gaussian Unitary Ensemble / Poisson spacing experiment

# GUE Protocol v1 — zeta/prime **spectral rigidity** test suite

A small, reproducible experiment suite for one question:

> **What “global constraints” turn a prime-based signal from Poisson-like spacings into GUE-like level repulsion?**

This repo implements **three model channels** (baseline → symmetry → ground truth) plus **two mandatory negative controls**, then runs the same zero-finding + unfolding + spacing statistics pipeline across all of them.

If you’re looking for a proof of RH: this is not that.  
If you’re looking for a **clean numerical microscope** for “rigidity vs collapse” in zeta-like spectra: this *is* that.

---

## Inputs

### Model channels
| Channel | What it is | What it’s testing |
|---|---|---|
| `independent_primes` | Partial Euler product baseline (no functional equation) | “Independent primes” ≈ **Poisson** baseline |
| `riemann_siegel` | Riemann–Siegel **Z(t)** approximation | Effect of **global symmetry** (functional equation–type structure) |
| `full_zeta` | `mpmath.zeta(1/2 + it)` mapped to **Z(t)** | “Ground truth” reference |

### Negative controls (non‑negotiable)
| Control | Why it exists |
|---|---|
| `fake_primes_jitter` | Same pipeline, but primes are “jittered” → breaks arithmetic coherence |
| `rs_phase_randomized` | Same RS amplitude scaling, but phases randomized → destroys coherence |

---

## Quickstart

```bash
pip install -r requirements.txt
python3 run_full_experiment.py --t-min 10 --t-max 200 --dt 0.02 --max-zeros 120 --p-max 2000
```

Outputs go to `out_gue/` by default:
- `zeros_<model>.csv` — zero locations `t`
- `summary.json` — metrics per model

### CLI options (most useful)
- `--t-min`, `--t-max`, `--dt`, `--max-zeros` — scan window + resolution
- `--p-max` — maximum prime for the Euler-product baseline
- `--k-max` — Euler product power truncation (k ≥ 1)
- `--seed` — determinism for negative controls
- `--outdir` — output folder

---

## What the pipeline actually does

### 1) Zero finding (robust, boring, correct)
- Scan for sign changes on a grid
- Refine each bracket by **bisection** (no free Newton wandering)

Implementation: `guesuite/zeros.py`.

### 2) Unfolding (density normalization)
Given zeros at heights \(t_n\), spacings are unfolded via the standard density factor:
$$
s_n = (t_{n+1}-t_n)\,\frac{\log(t_n/2\pi)}{2\pi},
$$
then normalized so \(\mathbb{E}[s]\approx 1\).

Implementation: `guesuite/stats.py`.

### 3) Metrics
We report:
- **KS vs Poisson** (exponential spacing CDF)
- **KS vs GUE Wigner surmise** (β=2):
 $$
  P_{\mathrm{GUE}}(s)=\frac{32}{\pi^2}s^2\exp\!\left(-\frac{4}{\pi}s^2\right)
 $$
- Mean **r-statistic**:
 $$
  r_n=\frac{\min(s_n,s_{n-1})}{\max(s_n,s_{n-1})}.
 $$

---

## How to read the results (don’t overfit)

- With too few zeros, **everything is noise**.  
  If `spacing_summary.n < ~100`, treat KS p-values as “suggestive at best.”
- The “interesting” signal is not a single metric — it’s the **pattern across channels**:
  - baseline vs symmetry vs ground truth
  - and whether controls collapse as expected

---

## “Landscape vs Swampland” (metaphor, not a physics claim)

If you want a spicy framing (use it responsibly):

- **Rigid phase (“landscape”)**: coherent global constraints → level repulsion emerges
- **Collapsed phase (“swampland”)**: perturb the coherence (jitter phases/primes) → repulsion dies, clustering returns

This repo is a **rigidity probe**. It does *not* prove a swampland conjecture.

---

## Repo layout

```
guesuite/
  models.py        # Z(t) models + negative controls
  zeros.py         # sign-change scan + bisection refinement
  stats.py         # unfolding + KS + r-statistics
run_full_experiment.py
requirements.txt
protocol_v2_addendum.tex / .pdf
PATCH_NOTES_v2.md
```

---

## Protocol docs (recommended)

- **`protocol_v2_addendum.pdf` / `.tex`** — the lab protocol: what must be present for a run to be considered valid
- **`PATCH_NOTES_v2.md`** — minimal fixes that led to this v2 bundle (formula typos, non-negotiables)

---

## Common pitfalls

- **dt too large** → you miss sign changes → garbage zeros.
- **t_max too small** → too few zeros → meaningless KS.
- **p_max too small** → baseline degenerates (you’re not testing anything).
- **Comparing raw spacings without unfolding** → invalid.

---

## Roadmap (if you want this to look “paper-ready”)
- Add plots: spacing histogram + empirical CDF vs Poisson/GUE
- Add CI: fast smoke test with small `t_max`, verify JSON schema + no regressions
- Add precision controls: `mp.mp.dps` knob + timing logs
- Add more deformation knobs: jitter width as CLI arg; controlled phase noise level

---

## License
No license is included yet. If you want external contributions, add one (MIT/BSD/Apache-2.0).


- 
