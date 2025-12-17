# GUE / Poisson spacing experiment

This implements the "three-model" comparison plus two negative controls.

## Run

```bash
pip install -r requirements.txt
python run_full_experiment.py --t-min 10 --t-max 200 --dt 0.02 --max-zeros 120 --p-max 2000
```

Outputs land in `out_gue/`:
- `zeros_<model>.csv`
- `summary.json`

Notes:
- Increase `--t-max` to gather more zeros and stabilize the KS tests.
