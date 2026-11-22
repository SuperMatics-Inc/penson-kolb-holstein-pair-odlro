# Positive-temperature pair long-range order in the three-dimensional pi-flux Penson-Kolb-Hubbard-Holstein model with unbounded phonons

Davis Rens, SuperMatics Inc., [ORCID 0000-0001-5839-5197](https://orcid.org/0000-0001-5839-5197)

| | |
| --- | --- |
| Preprint | [10.5281/zenodo.24433375](https://doi.org/10.5281/zenodo.24433375) |
| Replay supplement | [10.5281/zenodo.24433377](https://doi.org/10.5281/zenodo.24433377) |
| Version | 1.0.0, 21 November 2026 |
| Manuscript | [`paper/main.pdf`](paper/main.pdf), 28 pages |

This repository is the tagged release tree of that preprint. It holds the reviewed manuscript
source, the released PDF and the complete replay supplement, pinned byte for byte through
[`MANIFEST.sha256`](MANIFEST.sha256) and [`replay/SHA256SUMS`](replay/SHA256SUMS).

## What the paper does

The paper extends a released theorem on explicit, constructive positive-temperature onsite pair
off-diagonal long-range order in three-dimensional Penson-Kolb-Hubbard models to on-site Holstein
phonons with an unbounded oscillator Fock space. Each site carries a quantum harmonic oscillator
of frequency `omega` that couples to the local charge. The pair-hopping amplitude is `J > 0`; the
static polaronic shift is `U_ph = g^2/omega`, with `g` the on-site electron-phonon coupling.

1. **Reflection positivity survives the phonons.** Graded reflection positivity, Gaussian
   domination and the infrared bound hold for every real `g`, every `omega > 0`, every real `U`
   and every real hopping amplitude, on the full unbounded fermion-phonon Hilbert space. The
   odd-sublattice particle-hole transformation used for the fermions must be accompanied by the
   odd-sublattice phonon parity, and the oscillator cutoff is then removed by a trace-norm bridge.
2. **The sum rule acquires a static polaronic term.** The double-commutator sum rule for the bare
   pair operator gains an exact term `+4 U_ph (1 - r)`, with `r` the vacancy density, together
   with a density of displaced-oscillator excitations bounded by three elementary estimates.
3. **Fifteen certified parameter points.** The sum rule is closed with a Gibbs variational
   polaron-dimer bound, an energy-entropy decomposition at a split phonon frequency and the
   vacancy partition function. The fifteen rows of the manuscript's table cover adiabatic to
   anti-adiabatic phonons, on-site attraction from weak to strong, weak on-site repulsion shifted
   into the admissible window by the phonons, and the pi-flux lattice with `|kappa|/J <= 0.50`.
   Each row carries a temperature ceiling `T_0/J` and a positive lower bound on the staggered pair
   order parameter, certified by 50-digit outward-rounded interval arithmetic and valid at every
   temperature `0 < T <= T_0`. The oscillator cutoff is taken to infinity first at fixed even
   lattice size, then the lattice size to infinity.
4. **Where the bound stops.** For `U = 0`, `kappa = 0`, `alpha = 1/2` and `T_0 = 0.06 J` the
   certified floor is positive only for `g/J <= 0.3169`, `0.4568` and `0.7499` at `omega/J = 1`,
   `2` and `5`, that is for `U_ph/J <= 0.1004`, `0.1043` and `0.1125`. At the next coupling in the
   bisection the method returns no floor. The loss is not the static term: at `omega = J` the
   static part contributes `0.0271` and the displaced-oscillator part `0.0447` to the cost `c`,
   the latter being the piece that grows with the coupling. `replay/probe_xover.py` reproduces the
   boundary and `replay/probe_static.py` reproduces the decomposition.

## What the paper does not claim

The bound is a lower bound on a staggered pair correlation at explicit temperatures, and nothing
more. No critical temperature, no phase transition and no upper bound on a critical temperature is
proved or excluded, and no statement is made about a measured transition temperature anywhere. No
material is claimed, no identification of the model with any superconductor is claimed, and the
order parameter bounded here is not a measured quantity. The interaction is a pair hopping, not an
electronic Coulomb interaction.

The sum rule is written for the bare pair operator. Whether the polaron-dressed pair operator has
long-range order at positive temperature is an open problem, and no statement about a dressed
operator is made. Only the pi-flux lattice with `|kappa|/J <= 0.50` is covered in the hopping
regime, together with the zero-hopping case; the zero-flux lattice with `kappa != 0` stays open,
because the reflection-positive decomposition fails there and no other route to a spatial pair
infrared bound is available. Polaronic and charge-density-wave physics at strong coupling is not
treated.

The fifteen rows are the certified parameter points, not a continuum: a parameter point absent
from the table is a point at which this method has not been certified, not a point at which the
order parameter vanishes. The interval enclosures and the finite-matrix checks verify finite
arithmetic and small-system identities; they do not formalise the analytic or thermodynamic
arguments that carry the theorems.

## Contents

| Path | What it is |
| --- | --- |
| `paper/main.pdf` | The released manuscript, 28 pages. |
| `paper/main.tex` | The LaTeX source that produces that PDF. |
| `paper/refs.bib` | The bibliography. |
| `BUILD.md` | How the PDF is built, with the compiler version, the fixed epoch and the observed log counts. |
| `replay/README.md` | What each replay script checks, how to run it and its byte custody. |
| `replay/replay.py` | The Brillouin-zone interval enclosures, the fifteen certified rows, the exact finite-matrix operator identities, the reflection-positivity Gram matrices at oscillator cutoffs, the two obstruction checks and the seven in-script mutants. |
| `replay/probe_xover.py` | The positivity boundary of the certified floor in the Holstein coupling, by bisection. |
| `replay/probe_static.py` | The decomposition of the sum-rule cost `c`, rectangle by rectangle. |
| `replay/mutants.py` | The independent mutation audit, one process per mutant, with the expected diagnostic. |
| `replay/run_all.py` | Standard-library driver that runs the four commands and reports each exit code. |
| `replay/probe_xover.txt`, `replay/probe_static.txt` | The saved output of the two probes. |
| `replay/SHA256SUMS` | Pins for the eight other files of the replay directory. |
| `MANIFEST.sha256` | Pins for every other tracked file of this repository. |
| `CITATION.cff` | Citation metadata for the preprint and the supplement. |
| `LICENSE` | Apache-2.0 for the code, CC-BY-4.0 for the manuscript, as stated in its first lines. |
| `.github/workflows/ci.yml` | The two hosted jobs described below. |

`paper/main.tex` differs from the admitted reviewed source in one place only, the data
availability paragraph, which now carries the two archival identifiers of this release and the
repository URL. Nothing else differs from the admitted reviewed source; the staging record
compares the remainder line for line. `paper/refs.bib` is unchanged.
`replay/` is unchanged, and its own `SHA256SUMS` still verifies.

## Reproduce

### The replay

Python 3.12 or later. The replay scripts import only NumPy, SciPy and mpmath, at the verified
versions:

```sh
python -m pip install numpy==2.5.5 scipy==1.18.1 mpmath==1.4.1
```

From the `replay/` directory of an unpacked release:

```sh
cd replay
python -B run_all.py
```

`run_all.py` runs `replay.py`, `replay.py --mutants`, `probe_xover.py` and `mutants.py`, prints a
summary line per step with its exit status and wall time, and exits nonzero unless all four pass,
so a partial run cannot look like a pass. The individual commands are:

```sh
python -B replay.py            # interval certificates, finite checks, seven in-script mutants
python -B probe_xover.py       # positivity boundary in the Holstein coupling
python -B probe_static.py      # decomposition of the cost c at the boundary coupling
python -B mutants.py           # independent mutation audit
```

Timings vary with the machine and its load. The probe scripts check their own printed values
against the manuscript and exit nonzero on disagreement. [`replay/README.md`](replay/README.md)
maps every check to the statement it backs.

### The paper

```sh
cd paper
SOURCE_DATE_EPOCH=1795132800 tectonic -X compile --keep-logs main.tex
```

Tectonic 0.17.0 with that fixed epoch reproduces the released PDF byte for byte (SHA-256
`a651bdf1bf846b4db89c5922c7b5e93ab79cd69af16cf5eb71d815ca315d212a`, 265507 bytes, 28 pages). The
log has no overfull box, no underfull box, no undefined reference and no warning; see
[`BUILD.md`](BUILD.md).

### Continuous integration

[`.github/workflows/ci.yml`](.github/workflows/ci.yml) has two jobs. The `replay` job verifies
`MANIFEST.sha256` and `replay/SHA256SUMS`, installs the pinned NumPy, SciPy and mpmath, and runs
`run_all.py`, that is the certificates, the in-script mutants, the positivity-boundary probe and the
independent mutation audit. The `paper` job builds the manuscript with the epoch above, requires zero overfull
boxes, zero undefined references and zero warnings, requires exactly 28 pages, and requires both
DOIs in the built PDF.

## Byte custody

Every released file is a copy of reviewed bytes, and `MANIFEST.sha256` plus `replay/SHA256SUMS`
pin all of them except the two manifests themselves:

```sh
sha256sum -c MANIFEST.sha256
(cd replay && sha256sum -c SHA256SUMS)
```

The deposited supplement archive is `git archive` of the `v1.0.0` tag, so its members are the
files listed by `git ls-files` at that commit, unchanged.

## Licences

The manuscript (`paper/main.tex`, `paper/refs.bib`, `paper/main.pdf`) is licensed under Creative
Commons Attribution 4.0 International (CC-BY-4.0). The replay scripts and all remaining files are
licensed under the Apache License 2.0. The first lines of [`LICENSE`](LICENSE) state that split;
the full Apache-2.0 text follows.

## Review and AI assistance

AI agents generated proofs, certificates, code, and manuscript drafts and edits under D. Rens'
direction and review. Numerical enclosures used mpmath interval arithmetic, and the linear algebra
checks used NumPy and SciPy. The analytic statements are presented with complete arguments, but no
hand verification of every algebraic step was performed. Two independent AI reviews of the
repaired draft returned PASS on the released bytes. This preprint has not been peer reviewed.
