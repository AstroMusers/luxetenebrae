# Luxetenebrae

We use COMPAS to simulate stellar binaries for population studies and detailed analyses of selected systems.

Our `paths.py` module reads `LUXETENEBRAE_PATH`, which defaults to this checkout. It exposes local `data/` and `visuals/` directories. Git ignores both.

Our checkout-level scripts use TDpy's shared repository-path utility. The conda environment includes it, or we can install the workflow dependencies from this checkout.

```bash
python -m pip install -e ".[workflow]"
```
