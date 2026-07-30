# Reproducibility: what runs where, and how to rebuild it

This is the authoritative account of the software environments this analysis runs
in, how they are recorded, and how a reader rebuilds one. Everything below was
verified on 2026-07-30 rather than inferred; where something is inferred it says so.

## The three environments

The analysis runs in three places, and they are not the same. Pretending otherwise
is the single most likely way for a result to become irreproducible, so the
differences are recorded rather than wished away.

| | R | Library | Purpose |
|---|---|---|---|
| **ARC** (environment of record) | 4.4.2, module `R/4.4.2-gfbf-2024a` | 190 packages in `$DATA/PROJECT_GROUP/R/library_4.4` | Produces every reported result |
| **CI** | 4.4.2, GitHub Actions | Installed per run from Posit's binary repository | Static checks and unit tests |
| **Workstation** | whatever is installed | whatever is installed | Editing, and the parts of the suite that need no Stan toolchain |

Only ARC fits models. A workstation typically lacks brms and the Stan toolchain
entirely, which is expected and supported, but it means a green local test run
exercises materially less than CI does. See "Watertight checks" below.

## What records the environment

Three artefacts, with distinct jobs. Confusing them is what let the environment
drift unnoticed for as long as it did.

**`R/12_environment.R` — the requirement.** Declares, in one place, which packages
the analysis needs. This is a specification of what *must* be present, not a
description of what is. The groups matter:

| Group | Contents | Needed for |
|---|---|---|
| `core` | testthat, here, withr, dplyr, tibble, forcats, readr, assertr, yaml, purrr, stringr, tidyr | Loading and unit-testing the deterministic logic |
| `modelling` | brms, posterior | Loading the analysis modules, hence running the suite at all |
| `fitting` | cmdstanr, MASS, ordinal, bridgesampling | Actually sampling a model |
| `pipeline` | optparse, targets, tarchetypes | The runnable scripts and `{targets}` |
| `reporting` | quarto, bib2df, httr2, ggplot2 | Rendering the report, auditing the bibliography |

`modelling` is kept deliberately narrow because it is the group the test gate
enforces; demanding more would block environments that can run every test. `fitting`
is separate because **cmdstanr is not on CRAN** — it comes from the Stan repository
at <https://mc-stan.org/r-packages/> — so CI reasonably has brms without it, and the
unit tests never sample anything. An ARC node about to spend days fitting should
check everything with `--groups all`.

**`renv.lock` — the pinned versions.** Records exact versions and content hashes
for **176 packages** under R 4.4.2. Regenerated on ARC on 2026-07-30 with
`renv::snapshot(type = "all")` against the live library, so every version in it
matches `config/arc_environment_recorded.csv` exactly: 176 locked, 176 matching,
none differing. 176 rather than 190 because renv excludes the ~14 base packages
that ship with R itself, which the module pin already fixes.

**`config/arc_environment_recorded.csv` — the observation.** The 190 packages
actually present on ARC, with the library each resolves from, written by
`scripts/record_environment.R`. This is evidence, not intent: it is what the
reported results were produced under.

### Regenerating the lockfile

It must be regenerated **on ARC**, because that is where the library it describes
lives. Doing it anywhere else records the wrong machine.

```bash
ssh <arc-host>
cd ~/design_analysis
module load R/4.4.2-gfbf-2024a
export R_LIBS_USER=$DATA/PROJECT_GROUP/R/library_4.4
export RENV_PATHS_CACHE=$DATA/PROJECT_GROUP/renv/cache
Rscript -e 'renv::snapshot(project = ".", library = .libPaths(), \
                           lockfile = "~/renv_new.lock", type = "all", \
                           prompt = FALSE, force = TRUE)'
```

Three arguments are doing real work and are not optional:

- `library = .libPaths()` snapshots **both** libraries. Restricting it to the
  project library alone fails pre-flight validation, because roughly a hundred
  dependencies live in the R module's site library and renv correctly reports them
  as unsatisfied.
- `type = "all"` records everything installed rather than only what renv infers
  from the source, which is what makes the lockfile able to rebuild the
  environment rather than merely describe its direct dependencies.
- `force = TRUE` bypasses one genuine pre-flight failure: `pkgdown` requires
  httr2 (>= 1.0.2) while the project library pins **httr2 1.0.1**, which shadows
  the site library's 1.0.6. Neither `pkgdown` nor `devtools` is referenced anywhere
  in the analysis — they are incidental tooling — and httr2 1.0.1 is the version
  the reference audit has always used. The inconsistency is therefore recorded
  rather than repaired: upgrading httr2 would change a package the analysis
  actually depends on, to satisfy one it does not.

Writing to a new path and reviewing before replacing `renv.lock` is worth the extra
step; validate with `scripts/check_environment.R --lockfile <path>`.

### Known limits of the lockfile

- **`cmdstanr` is not on CRAN.** It is locked from an r-universe mirror
  (`https://bbsbayes.r-universe.dev`) with the upstream `stan-dev/cmdstanr` remote
  and commit SHA recorded alongside. `renv::restore()` uses the per-package
  `Repository` field, so this resolves, but that mirror must be reachable. If it is
  not, install CmdStan the ordinary way with
  `cmdstanr::install_cmdstan(cores = 4)`.
- **The lockfile pins R packages, not the Stan toolchain.** CmdStan 2.39.0 and GCC
  13.3.0 are recorded in `config/arc_modules.yaml` and in the environment record,
  not in `renv.lock`. Stan models are compiled C++, so the compiler is part of what
  determines whether a fit reproduces exactly.
- renv is **not activated**: there is no `.Rprofile` and no `renv/` directory, so
  nothing forces a session to obey the lockfile. It is a specification that is
  *checked* rather than one *enforced* by library redirection.

That last point is deliberate. Activating renv would override `R_LIBS_USER`, which
every `hpc/` script sets explicitly to the project library on `$DATA`, and would
break the working cluster setup. On a cluster where the library is managed outside
the project, verification is the appropriate mechanism; library hijacking is not.

## Why not a container

Containers would be the stronger answer, and they are not available here.
Checked on ARC on 2026-07-30: `apptainer`, `singularity`, `podman` and `docker` are
absent from `PATH` and from the module tree on the login node. Nothing in this
repository can therefore be containerised for execution on this cluster today.

Reproducibility instead rests on three things together: the module pin in
`config/arc_modules.yaml`, the version pins in `renv.lock`, and the observed record
in `config/arc_environment_recorded.csv`. If ARC gains Apptainer, a container built
from the recorded package set would be a strict improvement and worth doing.

## Watertight checks

The test suite used to report success in an environment where it could not run.
Three test files source modules that call `library(brms)`; without brms those files
fail to load, and testthat reports a load failure as a **skip**. The suite printed

```
[ FAIL 0 | WARN 0 | SKIP 4 | PASS 37 ]
```

while 29 of its checks had not executed. Nothing in that line distinguishes "all
passed" from "most never ran".

Two changes fixed it:

1. **`tests/testthat/test-environment.R` fails the run** when a required `core` or
   `modelling` package is missing. It sorts first, so it reports before anything
   else, and it prints the environment and the exact missing packages.
2. The three affected files now carry an explicit `skip_if_not_installed("brms")`
   before their `source()`, so their skip is a stated dependency rather than an
   accident with a misleading reason.

A partial run is still possible, and still useful while editing:

```bash
CLAPS_ALLOW_PARTIAL_TESTS=1 Rscript tests/testthat.R
```

That downgrades the gate to a skip and prints a prominent `!! PARTIAL TEST RUN !!`
banner. CI sets nothing, so CI cannot pass with an incomplete environment.

## Linting

The lint standard lives in `.lintr.R`, next to the code, so that

```bash
Rscript -e 'lintr::lint_dir("R")'
```

reproduces CI exactly. It used to be passed inline in the workflow YAML, which meant
running lintr locally applied different rules, and because the step only warned, 259
findings had accumulated unread. The step now **fails** on any finding, and `R/` is
clean against the configuration.

Of those 259, some 60% were `object_usage_linter` reporting dplyr's NSE column names
and functions defined in `source()`d files as undefined. Those are artefacts of this
project's architecture rather than defects, so the linter is disabled with that
reasoning recorded in `.lintr.R`; the remaining 87 were real formatting issues and
were fixed. Every deviation from lintr's defaults is argued in that file. If a new
rule needs relaxing, argue it there rather than restoring a warn-only step.

## Checking an environment

```bash
Rscript scripts/check_environment.R                                   # core + modelling
Rscript scripts/check_environment.R --groups all                      # everything
Rscript scripts/check_environment.R --compare config/arc_environment_recorded.csv \
                                    --lockfile renv.lock              # full report
```

Exits 0 if the requested groups are satisfied and 1 otherwise, so it can gate a
shell pipeline or the head of a SLURM job. Drift against the recorded environment
and disagreement with the lockfile are **reported, not fatal**: a workstation is
legitimately not the cluster, and the value is in the difference being visible in
the log rather than in forcing the two to be identical.

Both scripts use base R only, deliberately. They must run in exactly the broken
environments they exist to diagnose, so they cannot depend on the packages they are
checking for.

## Recording an environment

```bash
Rscript scripts/record_environment.R                        # -> outputs/
Rscript scripts/record_environment.R --out config/arc_environment_recorded.csv
```

Writes a CSV (package table, diffable in review, read by the checker) and a JSON
(the same table plus R version, platform, OS, compiler, CmdStan version, and the
SLURM job identifiers). Run it on ARC, inside or alongside a job, so the record
describes the machine that fits the models.

One subtlety it handles: `installed.packages()` returns one row per package **per
library**, so a package present in both the project library and the R module's site
library appears twice with different versions. Seven do on ARC (`cpp11`, `httr2`,
`knitr`, `ps`, `rlang`, `systemfonts`, `yaml`). R resolves these by search order, so
the record keeps the copy that actually loads, names its library, and lists the
shadowed copies in the preamble. A shadowed newer version is a common reason for a
result differing between two accounts on the same cluster, so the information is
kept rather than discarded.

## Rebuilding the ARC environment from scratch

```bash
ssh <arc-host>
srun --partition=devel --cpus-per-task=4 --mem=16G --time=01:00:00 --pty bash
module purge && module load R/4.4.2-gfbf-2024a

export R_LIBS_USER=$DATA/PROJECT_GROUP/R/library_4.4
export RENV_PATHS_CACHE=$DATA/PROJECT_GROUP/renv/cache
mkdir -p "$R_LIBS_USER" "$RENV_PATHS_CACHE"

cd ~/design_analysis
Rscript -e 'renv::restore(prompt = FALSE)'          # pinned versions from renv.lock
Rscript -e 'cmdstanr::install_cmdstan(cores = 4)'   # CmdStan 2.39.0
Rscript scripts/check_environment.R --groups all
bash hpc/submit_devel_smoke_test.sh
```

Never build on a login node: compiling CmdStan there is both slow and against
cluster policy. The `check_environment.R` call is the point at which you find out
whether the rebuild worked, before committing days of compute to it.

## Provenance of a result

Every run writes `outputs/manifest.csv` via `write_manifest()` in
`R/10_job_status.R`, recording the git SHA, timestamp, and the R, brms and cmdstanr
versions. `scripts/record_environment.R` captures the fuller picture when a complete
record is wanted. Between the git SHA, the manifest and the recorded environment, a
result can be traced to the code and the software that produced it.
