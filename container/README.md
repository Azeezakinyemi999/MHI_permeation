# Hydrogen Permeation Model — Offline Analysis Environment

A self-contained JupyterLab environment for the MHI hydrogen-permeation model.
Everything needed to run the analysis is inside the container image. **No Python
packages are installed, and no internet access is required, at any point.**

Verified to run with networking fully disabled.

---

## What you received

```
Hydrogen_Model_Environment_1.1.0/
├── image/
│   ├── hydrogen-model-1.1.0.tar      the container image
│   └── hydrogen-model-1.1.0.sha256   checksum for the above
├── scripts/
│   ├── verify.sh                     run the two acceptance gates
│   └── start.sh                      launch JupyterLab
├── workspace/                         your working files (see below)
├── requirements.in                    what we deliberately depend on
├── requirements.lock.txt              the exact 106 packages the image installs
├── installed-packages.txt             pip freeze taken from inside the built
│                                      image; byte-identical to the lock, which
│                                      is the evidence the build honoured it
├── Dockerfile                         how the image was built
└── README.md
```

## Requirements

- A Docker-compatible container runtime (Docker Engine, Docker Desktop, or
  Podman in Docker-compatible mode) able to run **Linux `amd64`** containers.
- ~3 GB free disk for the loaded image.
- Permission to bind-mount a local directory and to bind a port on `127.0.0.1`.

The container runs as a **non-root** user and needs no elevated privileges.

## First-time setup

**1. Check the archive is intact.** From the package root:

```bash
shasum -a 256 -c image/hydrogen-model-1.1.0.sha256
```
Expect `image/hydrogen-model-1.1.0.tar: OK`. On Linux without `shasum`, use
`sha256sum -c`. If this fails, the transfer was corrupted — do not proceed.

**2. Load the image.**

```bash
docker load -i image/hydrogen-model-1.1.0.tar
docker images hydrogen-model
```

**3. Verify it.**

```bash
./scripts/verify.sh
```
This runs with `--network none` and must end with:

```
VERIFICATION PASSED — hydrogen-model:1.1.0 is the validated environment.
```

It checks two separate things: that all 106 pinned packages are exactly the
validated versions, and that the model still reproduces known reference values
through its real solvers. If either fails, do not use the image for analysis.

## Daily use

```bash
./scripts/start.sh
```

Then **copy the full URL printed in the terminal, including the `?token=...`
part**, into your browser. Without the token the page will not load. The token
is new each time the server starts.

Press `Ctrl-C` in the terminal to stop.

## The workspace

```
workspace/
├── calculations/          the model. Python package; do not rename.
└── Application/
    ├── Proposal.ipynb                  levels 1-5, self-contained
    ├── Surface_proposal.ipynb          surface kinetics (level 6), self-contained
    ├── sensitivity_regime_L5.ipynb     regime sensitivity analysis
    ├── sensitivity_regime_L5L6.ipynb   regime sensitivity incl. surface
    ├── regime_parallel_coords.ipynb    visualises sensitivity_regime_L5L6 output
    ├── regime_parallel_coords_L5.ipynb visualises sensitivity_regime_L5 output
    ├── sa_results/                     results land here (starts empty)
    └── sa_results_L5/                  results land here (starts empty)
```

Anything you put in `workspace/` is visible inside JupyterLab, and anything the
notebooks write appears on your own filesystem. Results are written **next to the
notebook that produces them**, not to a central output directory.

### Run order matters

`Proposal.ipynb` and `Surface_proposal.ipynb` are self-contained — run them in
any order.

The other four are a pipeline. The two `regime_parallel_coords*` notebooks only
*plot* results; they read CSVs produced by the sensitivity notebooks and cannot
regenerate them. **No precomputed results are shipped**, so run the producer
first:

```
sensitivity_regime_L5.ipynb    ──>  regime_parallel_coords_L5.ipynb
sensitivity_regime_L5L6.ipynb  ──>  regime_parallel_coords.ipynb
```

Running a `parallel_coords` notebook first fails with
`FileNotFoundError: sa_results_L5/master_clusters.csv`. That is expected, not a
fault — run its producer.

The sensitivity notebooks run a large parameter scan and take a while. Once they
have run, their CSVs are cached in `sa_results*/scans/`, and re-running the plots
is fast.

## Adding a Python package

Don't `pip install` inside a running container — the change is lost when the
container stops, and it silently breaks the guarantee that `verify.sh` checks.

Instead, tell us what you need. We add it to `requirements.in`, regenerate the
lock, rebuild, re-run both gates, and issue a new version. That is the whole
point of the pinned lock: the environment that produced a result can be
reconstructed exactly.

## Reporting a problem

Please include:

- the output of `./scripts/verify.sh` in full
- `docker --version` and your OS
- `docker image inspect hydrogen-model:1.1.0 --format '{{.Id}} {{.Architecture}}'`
- which notebook and which cell, with the error text

## Version

`hydrogen-model:1.1.0` — Python 3.12, linux/amd64. Cite this tag alongside any
results, so they can be tied back to a known environment.

### Changes since 1.0.1

Adds `kaleido` so plotly figures export to vector PDF/SVG offline — previously
the parallel-coordinates figures existed only as HTML and could not be placed in
a manuscript. `kaleido` is pinned to 0.2.1 deliberately: 1.x drives a real Chrome
browser, which is neither present in this image nor installable without network.
`plotly` is pinned alongside it for the same reason.

The model's figure styling also changed: figures are now drawn at Elsevier
column widths with fonts that stay above 7 pt at final printed size, vector
output, and colour-vision-safe palettes. Sensitivity figures now carry 95%
bootstrap intervals and the noise floor. No computed value changed — the
numerically significant packages are identical to 1.0.0.

### Changes in 1.0.1

1.0.0 required `--user "$(id -u):0"`; running it as `--user "$(id -u):$(id -g)"`
— the more natural form — failed with
`PermissionError: [Errno 13] Permission denied: '/home/appuser/.ipython'`.
1.0.1 redirects every per-user state directory (IPython, matplotlib, XDG) to
`/tmp`, so any `uid:gid` works and your files are written with your own group
rather than group 0. `scripts/start.sh` uses the natural form. Nothing about the
scientific environment changed — the 106 locked package versions and every
reference value are identical.
