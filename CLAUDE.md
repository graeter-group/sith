# CLAUDE.md

Guidance for working in the SITH repository.

## What SITH is

**SITH** (Splitting Intramolecular Tension via stretcHing) is a quantum-chemistry
analysis package. It splits the change in electronic energy of a mechanically
stretched molecule into contributions from individual internal degrees of
freedom (DOFs): bond lengths, angles and dihedrals. It does this by numerically
integrating the work-energy theorem along a COGEF-like stretching path:

    ΔE_i = -∫ F_i dq_i  ≈  -Σ_j (F_i^j + F_i^{j+1})/2 · Δq_i^{j→j+1}
    F_i^j = -∂E_DFT/∂q_i  evaluated at configuration j

Paper: *SITH: A quantum-chemical framework for predicting bond destabilization in
stretched molecules*, J. Chem. Phys. 164 (23), 2026, doi:10.1063/5.0325167.

The package also implements **JEDI** (harmonic approximation, which needs the
Hessian of the relaxed structure) as an alternative energy-distribution method.

## End-to-end workflow (from the docs and tutorials)

1. **Stretch** (bash, Gaussian): `sith stretching -c -i 1,16 -m G.pdb`
   - Optimize, then repeatedly move two atoms apart, constrain that distance and
     re-optimize until a bond breaks.
   - Produces `<mol>-stretched<n>.{com,log,chk,fchk,xyz}` (n=0 is the
     unconstrained optimum), `frozen_dofs.dat` (line 1: pulled atoms; line 2:
     atoms of the broken bond, each ending in `F`) and a `rupture/` directory.
   - An alternative path starts from the broken or extreme end:
     `workflow_from_extreme.sh`, with outputs in `from_extreme/`, which needs sklearn.
   - `workflow.sh` wraps stretching and force extraction and can resubmit itself
     on slurm.
2. **Extract forces** (bash, Gaussian):
   `sith find_forces -c -f G-stretched<n>.chk -p stretched -v`, written to `forces/`.
   This writes `.fchk` files in the "sith format" that `G09Reader` expects.
   Internally it runs `formchk -3` and merges in the forces in internal
   coordinates from the log.
3. **Analyse** (Python):
   ```python
   from sith import SITH
   sith = SITH('./forces')          # dir of *.fchk, sorted alphabetically = path order
   sith.killer(killElements='H')    # optional: drop DOFs, e.g. those involving H
   sith.sith_analysis()             # or sith.jedi_analysis()
   ```
4. **Visualize**: `sith.visualize.vmol.EnergiesVMol` (VPython plus VMol) or
   `sith.visualize.ngl.EnergiesNGL` (nglview), plus plotting helpers in
   `sith.utils.sith_plots`.

Tutorial data, including the glycine G example: https://github.com/Sucerquia/Gdata4tut

## Gaussian / cluster setup (required by the bash tools)

The bash tools call a shell function named `gaussian` that the user must define.
It takes `gaussian <input.com> [output.log]`. You can define it in either of two ways:
- Define it in the shell or job script, then `source $(sith <tool> -path) ...`.
- Put it, together with any `module load` lines, in `$HOME/sw/load_modules.sh`
  and pass `-c` to the tool (`sith stretching -c ...`).

The scripts carry `#SBATCH` headers and are designed for slurm:
`sbatch -n 8 -J name $(sith stretching -path) -c -i 1,37 -m mol.pdb`.
You cannot run these tools for real without Gaussian (g09 or newer).

## Code layout (`src/sith/`)

| Path | Role |
|---|---|
| `SITH.py` | `SITH` class, the main entry point. Loads structures through a reader, validates that all structures share atoms and DOFs, and holds the result arrays. Methods: `sith_analysis(integration_method=...)` (`trapezoid_integration` by default, also `simpson_integration` and `rectangle_integration`), `jedi_analysis()`, `killer()`, `rem_first_last()`, `reduce_path()`, `energies_error()`. |
| `Utilities.py` | `Geometry`, the per-structure container that readers must fill (`n_atoms`, `scf_energy`, `dims`, `dim_indices`, `dof`, `internal_forces`, `hessian`, `atoms`), plus colour and colourbar helpers. |
| `SithWriter.py` | `WriteSITH` / `write_sith_data`, which write SITH data back out to fchk-like files. |
| `readers/g09_reader.py` | `FileReader` reads a single `.fchk` into a `Geometry`; it doubles as the template for new readers. `G09Reader` reads a list of files or a directory. Readers are registered in `SITH._get_reader` (`readers_objs` dict). |
| `energy_analysis/` | `SithAnalysis` (numerical integration; `delta_q` holds consecutive differences with angles wrapped into [-π, π]) and `JediAnalysis`. |
| `g09_stretching/` | Bash pipeline (`stretching.sh`, `workflow.sh`, `find_forces.sh`, `extract_forces.sh`, `restart_g09.sh`, `proline_mod.sh`, ...) plus Python helpers (`g09_xyz`, `protonate`, `sith_tools`, `change_dofs`). `from_extreme/` holds the reverse-path workflow. |
| `utils/` | `basics.sh`, a shared bash library that other scripts source (labels, `verbose`, `warning`, `fail`, `finish`, `create_bck`, `load_modules`, ...). Also `tools.py` (CLI-exposed geometry utilities), `molecules.py` (`MoleculeSetter`, alignment, Gaussian input creation), `peptides.py` (`PepSetter`), `analysis.py` and `sith_plots.py` (peptide dataset statistics and plots), and `compare_siths.py`. |
| `visualize/` | `ngl.py` (nglview) and `vmol.py` (VPython). |
| `cli/main.py` | The `sith` console script (see below). |

### Data conventions
- All per-structure arrays are indexed as `[structure, dof]`. DOFs are ordered
  distances, then angles, then dihedrals. `dims = [total, n_dist, n_ang, n_dih]`,
  so slice with `dims[1]` etc.
- `dim_indices` has shape `(#DOFs, 4)`. Indices are **1-based** (the Gaussian
  convention), and 0 pads DOFs that have fewer than 4 atoms.
- Units: Hartree, Ångström, radians. Energies are relative to `reference`
  (default structure 0).
- The order of structures equals the sorted filename order of the `.fchk` files.

## The `sith` CLI

`cli/main.py` dispatches `sith <name> ...` through two dicts:
- `pymodules`: name → Python module. It runs the function of that name, with
  `--kw value` parsed through `eval`. `-h` prints the docstring.
- `sh_executers`: name → bash script path. It runs the script, and `-path` prints
  the absolute path so the script can be sourced or `sbatch`ed.
- `sith -h` lists all tools. `sith tests` runs pytest in `tests/`.

**To add a new CLI tool, register it in one of these dicts.** Each bash script
documents its flags in a `print_help()` heredoc, reachable with `-h`.

## Documentation (`doc/`)

- Sphinx with `sphinx_rtd_theme`, `napoleon`, `autodoc` and `sphinxcontrib.mermaid`.
  `conf.py` adds `../src` to `sys.path` and defines a custom `:bashscript:` role
  (for example `:bashscript:`sith.g09_stretching.stretching`` links to the
  script's doc page).
- Build: `cd doc && make html`. `BUILDDIR=.`, so the output goes to `doc/html/`
  and `doc/doctrees/`, which are **committed**. The GitHub Pages workflow
  (`.github/workflows/jekyll-gh-pages.yml`) rebuilds on push to `main`. It also
  installs `PlottingTool` and `vmol` from GitHub, plus `nglview`, `ipywidgets`
  and `vpython`.
- `doc/modules/bash_rsts_doc/*.rst` mirror each script's `print_help()` text, and
  `bash_rsts_scripts/*.rst` `literalinclude` the script source. If you change a
  script's flags, update the matching `bash_rsts_doc` page; no generator
  script lives in this repo.
- `doc/workflow/*.rst` (extraction, analysis, output, visualization) are mostly
  empty stubs. The real user docs are `install.rst` and the two tutorials in
  `doc/tutorials/`.
- Docstrings use numpy-ish sections with `=====` underlines (`Parameters`,
  `Return`, `Attributes`). Follow that style.

## Dev setup and testing

- Install: `pip install -e .`. `pyproject.toml` requires Python >=3.9 and lists
  only `ase`, `numpy`, `matplotlib` and `sphinx-rtd-theme` as dependencies.
- Optional imports that are *not* declared: `scipy`, `sklearn`, `nglview`,
  `ipywidgets`, `vpython`, `vmol`, `PlottingTool`, `pytest`.
- Run tests with `pytest tests` (or `sith tests`). Coverage is minimal:
  - `tests/test_importall.py` builds module names with `str(path).split('/')[9:]`,
    which depends on how deep the checkout sits. In this checkout it produces
    `sith.sith.Utilities` and **fails**. Fix the path logic, for example with
    `relative_to`, rather than trusting the result.
  - `tests/visualize/test_ngl.py` is a placeholder.
  - Test data: `tests/resources/G/` (9 `structure_00N.fchk` files and `G.pdb`
    for capped glycine). The helper `sith.variables4tests.G` points there.
- Code style: PEP8 (recent commits were PEP8 passes).

## Gotchas

- Version strings disagree: `pyproject.toml` and `doc/conf.py` say 1.0.0, the
  README describes "Version 2.0" changes, and `dist/` holds stale 0.0.5 and
  0.1.0 builds. `package.json` exists only for `gh-pages` deploys, and its JSON
  is invalid (a comma is missing).
- The only reader so far is Gaussian (`G09Reader`). The architecture is designed
  so readers for QChem, ORCA or GPAW can be added. A new reader must fill every
  "mandatory" `Geometry` attribute and be registered in `SITH._get_reader`.
- `sith_analysis` needs small steps between structures for the integration to be
  accurate. Large jumps in dihedrals, such as methyl hydrogens flipping, cause
  numerical artefacts, which is why the tutorials use
  `killer(killElements='H')`.
- `utils/analysis.py` and `utils/sith_plots.py` are research scripts for
  analysing peptide datasets and change often. They have heavy plotting
  dependencies (`PlottingTool`).
