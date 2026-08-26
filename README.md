# xas-batch-processor

Batch processing of X-ray absorption spectroscopy (XAS) data with [xraylarch](https://github.com/xraypy/xraylarch) — merge scans, calibrate energy against references, and generate publication-quality XANES plots without clicking through the Athena GUI.

Developed for beamline data collected at BNL NSLS-II (BMM and ISS), where each sample produces many repeated scans that need to be merged, calibrated, and compared.

- Documentation: https://xraypy.github.io/xraylarch
- Code: https://github.com/xraypy/xraylarch

> Citation: M. Newville, *Larch: An Analysis Package For XAFS And Related Spectroscopies*. Journal of Physics: Conference Series, 430:012007 (2013).

## What it does

The script `Larch_XAS.py` runs one of three workflows, selected by the `FILE_TYPE` constant at the top of the file:

| `FILE_TYPE` | Typical source | Workflow |
|---|---|---|
| `''` or `.dat` | BMM single scans (`''`), ISS transmission scans (`.dat`) | **Raw scans → Athena project.** Reads raw beamline ASCII scans, computes transmission (`ln(I0/It)`) or fluorescence (`ΣIf/I0`) μ(E) — auto-detected from the scan header (`Scan.plot_hint`) or forced via `TRANSMISSION_MODE` — merges repeated scans of each sample, extracts the reference channel (`ln(It/Ir)`), and bundles everything into a single Athena project (`Created_transmission_group.prj`), followed by energy calibration. |
| `.prj` | BMM fluorescence scans | **Merge & calibrate Athena projects.** Reads existing `.prj` files, merges the scans within each project, tags each merge with its element/edge (via `find_e0` and `guess_edge`), collects the merges into `Created_group.prj`, and then energy-calibrates. |
| `.txt` | Normalized data exported from Athena | **Plot XANES.** Reads a normalized μ(E) table and produces a stacked, publication-quality comparison plot, styled through a `.ini` config file. |

### Energy calibration

For every group whose label contains `foil` or `reference`, the script compares the measured E0 (`find_e0`) with the tabulated edge energy (`xraydb.xray_edge`) and applies the resulting shift to the matching sample scans. The calibrated dataset is saved as `Created_group_with_calibration.prj`. A warning is printed when the shift exceeds 3 eV — inspect those manually in Athena.

## Files

- **`Larch_XAS.py`** — main script; all switches live in the constants block at the top.
- **`athena_project.py`** — patched local copy of larch's Athena project reader/writer, imported instead of `larch.io`'s version (see below).
- **`larch_plot_config.ini`** — template config for the plotting workflow.

### Why a local `athena_project.py`?

The stock `larch.io.athena_project` writer assigns every group a **random hash key** when saving a project, so the sample names are lost and the groups show up with meaningless labels when the `.prj` is opened in Athena. The `add_group()` method in the local copy was rewritten to accept an explicit `scan_name` and use it as the group key, so saved projects keep human-readable sample names:

```python
new_merge_project.add_group(merges, scan_name)   # stock larch add_group() has no name argument
```

The local copy is pinned to an older larch API (`_athena_groups`, `_larch=` keyword) that `Larch_XAS.py` relies on, so keep importing this file rather than swapping in the installed package's version.

## Requirements

- Python 3.9+
- [xraylarch](https://xraypy.github.io/xraylarch/), `xraydb`, `numpy`, `matplotlib`, `palettable`, `colorama`, `asteval`
- wxPython (the script sets matplotlib's backend to `wxAgg`)

**Recommended: install from conda-forge.** xraylarch and especially wxPython have compiled dependencies that conda-forge ships prebuilt, which avoids most installation headaches. You don't need the full Anaconda distribution for this — a minimal [Miniforge](https://github.com/conda-forge/miniforge) install (conda preconfigured for conda-forge) is enough:

```
conda create -n xas python=3.12
conda activate xas
conda install -c conda-forge xraylarch wxpython palettable colorama
```

Or with [pixi](https://pixi.sh/) (also uses conda-forge packages):

```
pixi init && pixi add python=3.12 xraylarch wxpython palettable colorama
```

Installing with pip / [uv](https://docs.astral.sh/uv/) works too, but wxPython may need to compile from source on Linux:

```
uv pip install xraylarch wxpython palettable colorama
```

(`xraydb`, `numpy`, `matplotlib`, and `asteval` are pulled in automatically as xraylarch dependencies.)

## Usage

1. Open `Larch_XAS.py` and edit the constants block:
   - `FILE_TYPE` — choose the workflow (see table above)
   - `INPUT_PATH` — folder containing your data files
   - `SKIP_SCANS` — names of bad scans to exclude from merging
   - `CONFIG_FILE` — path to the plotting `.ini` (only used for `.txt` mode)
   - `IF_SAVE` — save figures as 300-dpi PNGs. Applies to all three workflows: the per-sample merge plots (saved in `INPUT_PATH`) and the `.txt` comparison plot (saved next to the `.ini`). Set it to `False` to skip writing PNGs entirely — the `.prj` output is written either way. (The merge workflows close their figures immediately, so with `IF_SAVE = False` they produce no visible plot; only the `.txt` workflow displays one.)
2. Run the script with the environment you installed into (see [Requirements](#requirements)):
   ```
   # conda / Miniforge
   conda activate xas
   python Larch_XAS.py

   # pixi (from the project folder, no activation needed)
   pixi run python Larch_XAS.py

   # uv
   uv run python Larch_XAS.py
   ```
   If you run it from an IDE (PyCharm, VS Code, Spyder), point the project's Python interpreter at that environment instead.
3. Processed output (`*_merged`, `*_reference`, `Created_*.prj`) lands in an `Output_files/` subfolder of `INPUT_PATH`.

### Plot configuration (`.ini`)

For the `.txt` plotting workflow, copy `larch_plot_config.ini` next to your data and adjust:

- **`[samples]`** — `file_index` picks which file in the folder to plot; `sample_list` selects and orders the data columns; indices in `standard_list` are drawn as dashed lines (e.g., reference foils).
- **`[legends]`** — `sample_label` overrides the column names in the legend (LaTeX like `Cu$_2$O` and unicode like `°` work).
- **`[format]`** — figure size, [palettable](https://jiffyclub.github.io/palettable/) color palette, line widths, per-curve y-offset for stacked plots, energy/y-axis ranges, tick interval, and output filename.

Run the script once with an empty `sample_list` to print every column index and name, then use those indices to set up the plot you want.
