# Replay supplement: pair long-range order in the 3D pi-flux Penson-Kolb-Hubbard-Holstein model

This directory contains the replay script for the preprint on positive-temperature pair
off-diagonal long-range order in the three-dimensional pi-flux Penson-Kolb-Hubbard
model with on-site Holstein phonons.

| Item | Value |
| --- | --- |
| Script | `replay.py` |
| SHA-256 of `replay.py` | `1395c7457086db77efcb9c8e0478830b11491e395a6a9849250d741390d6e53b` |
| Companion scripts | `probe_xover.py`, `probe_static.py`, `mutants.py`, `run_all.py` |
| Saved output | `probe_xover.txt`, `probe_static.txt` |
| Digest list | `SHA256SUMS`, covering all eight files in this directory |
| Language | Python 3 |
| Dependencies | numpy, scipy, mpmath, sympy |
| Runtime | about 30 s for the certificates, about 2 min for everything |
| Exit status | 0 when every check passes; 1 on any failed check |

## How to run

```
PYTHONDONTWRITEBYTECODE=1 uv run --no-project --with numpy --with scipy \
  --with mpmath --with sympy python replay.py
PYTHONDONTWRITEBYTECODE=1 uv run --no-project --with numpy --with scipy \
  --with mpmath --with sympy python replay.py --mutants
```

Notation differs between the paper and the script. `replay.py` calls the pair-hopping
amplitude of the paper `g` and fixes it to 1, and it calls the on-site Holstein coupling
of the paper `g_ph`, printed as `gph`. The printed ceilings `T/g` are therefore the
paper's `T_0/J`, and a printed row such as `gph=0.25` is a Holstein coupling `g=0.25J`.

The first command prints the enclosures of the three lattice integrals, the finite
operator and reflection-positivity checks, and, for each of the fifteen parameter rows of
the paper, the temperature ceiling `T/g`, the weight `alpha`, the floor `m_0^2` and the
exact lower endpoint of the 50-digit enclosure from which the floor was truncated.
The second command runs the adversarial mutants; each mutant must exit with status 1,
and the suite reports success only if all of them are rejected.

The third script locates the positivity boundary of the certified floor that the paper
quotes in its discussion of strong coupling. It imports the unchanged `replay.py` and
bisects the Holstein coupling, asking only whether the certificate still returns a
positive floor, with the same 300 by 16 rectangle cover used for the fifteen rows:

```
python probe_xover.py
```

It prints the largest Holstein coupling `g_ph/J` with a positive floor, and the
corresponding `U_ph/J = g_ph^2/omega`, at `omega/J = 1, 2, 5` for zero hopping, `U = 0`,
`kappa = 0`, `alpha = 1/2` and `T_0 = 0.06J`, and then the fraction of that boundary at
which rows 6 and 7 of the paper sit. The values it prints are
`g_ph/J <= 0.3169, 0.4568, 0.7499`, that is `U_ph/J <= 0.1004, 0.1043, 0.1125`, and
93.3 and 87.0 per cent for the two rows. It exits 0 when the printed values agree with
the paper and 1 otherwise; `probe_xover.txt` is its saved output.

The fourth script decomposes the cost `c` of the sum rule at the minimising rectangle of
the largest certified coupling, for the three frequencies quoted in the discussion of
strong coupling. It repeats the interval certificate of `replay.py` rectangle by
rectangle, checks that the two floors agree exactly at every parameter point, and prints
the static and Cauchy-Schwarz fluctuation parts of `B_ph/(3J)`, the temperature piece and
the total `c` at that rectangle (the hopping piece vanishes at `kappa = 0`):

```
python probe_static.py
```

It prints `r ~ 0.190, 0.186, 0.178`, static parts `0.0271, 0.0283, 0.0308`, fluctuation
parts `0.0447, 0.0467, 0.0502`, temperature pieces `6.1e-4` and total `c` of
`0.0724, 0.0756, 0.0817` at `omega/J = 1, 2, 5`. It exits 0 when the mirrored floor
matches `replay.py` exactly and every printed number agrees with the paper to four
significant digits, and 1 otherwise; `probe_static.txt` is its saved output.

The remaining companion scripts sit next to `replay.py`; `probe_xover.py`,
`probe_static.py` and `replay.py` need numpy and mpmath, `mutants.py` and `run_all.py`
need nothing beyond a Python interpreter:

```
python run_all.py          # replay.py, replay.py --mutants, probe_xover.py, mutants.py
python probe_xover.py      # positivity boundary of the certified floor
python probe_static.py     # decomposition of c at the boundary coupling
python mutants.py          # mutation audit on its own
python mutants.py --list   # what each of the seven mutants corrupts
```

`run_all.py` prints the exit status and the wall time of each step and exits 0 only if
all four steps pass, so a partial run cannot look like a pass. `mutants.py` runs one
mutant per process and checks not only that the run exits with status 1 but also that the
rejection reports the quantity and the sign the paper says is broken, so that a mutant
cannot be rejected by an unrelated numerical accident.

| File | Contents |
| --- | --- |
| `replay.py` | interval certificates for the fifteen rows, the lattice integrals, the finite reflection-positivity checks, and the seven mutants |
| `probe_xover.py` | positivity boundary of the certified floor in the Holstein coupling, for the three frequencies and the two rows the paper quotes |
| `probe_xover.txt` | saved output of `probe_xover.py` |
| `probe_static.py` | decomposition of `c` at the minimising rectangle of the boundary coupling, for the three frequencies the paper quotes |
| `probe_static.txt` | saved output of `probe_static.py` |
| `mutants.py` | independent audit of the seven mutants, one process each, with the expected diagnostic |
| `run_all.py` | runs all four commands and reports every exit status |
| `SHA256SUMS` | digests of the eight files in this directory |
| `README.md` | this file |

## What the paper does

The paper proves, for the pi-flux Penson-Kolb-Hubbard model with on-site Holstein
phonons of arbitrary coupling and frequency, that a staggered pair correlation is
bounded below by an explicit positive number at every temperature below an explicit
ceiling. The proof has four inputs: graded reflection positivity with an
odd-sublattice phonon parity, Gaussian domination, a bound that removes the oscillator
cutoff at fixed lattice size, and a sum rule with a static polaronic term. The two
limits are taken in the order oscillator cutoff to infinity at fixed even lattice size,
then lattice size to infinity along even sizes.

`replay.py` reproduces the numbers of the table and the finite checks. It recomputes
each floor with 50-digit outward-rounded interval arithmetic in `mpmath.iv`, so each
printed interval is enclosed between the printed endpoints and the floor quoted in the
paper lies below the true value. `probe_xover.py` reproduces the coupling boundary quoted
in the discussion of strong coupling, from the same certificate.

## What the paper does not claim

The paper does not claim an exact critical temperature, does not claim that any material
is a superconductor, and makes no statement about phases. The bound is a lower bound on
an order parameter at explicit temperatures, nothing more. In particular:

- The bare-pair sum rule that the argument uses gives a positive floor only up to a
  finite Holstein coupling. For zero hopping, `U = 0`, `kappa = 0`, `alpha = 1/2` and
  `T_0 = 0.06J` the certificate is positive only for `U_ph/J <= 0.1004, 0.1043, 0.1125`
  at `omega/J = 1, 2, 5`, and above those couplings the paper gives no bound. That limit
  belongs to the method, not to the model.
- Polaronic and charge-density-wave physics in the strong-coupling regime is not treated.
- Dressed-pair long-range order, as opposed to the pair correlation analysed here, is an
  open problem.
- Only the pi-flux lattice is covered; the zero-flux model with nonzero single-electron
  hopping is not.
- The fifteen table rows are the parameter points for which the argument is carried out
  in full; they are not a continuum.
- The paper does not claim that every step was checked by hand.
