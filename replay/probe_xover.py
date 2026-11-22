# pyright: reportMissingImports=false
"""Locate the positivity boundary of the certified floor in the Holstein coupling.

This script imports the unchanged replay.py that sits next to it and bisects the
Holstein coupling g_ph for fixed U, kappa, omega, alpha and T_0, asking only whether
the certified floor of equation (5.14) of the paper is positive. The rectangle cover is
the one used for the fifteen certified rows: 300 cells in the vacancy direction and 16
in the pair-energy direction.

Two things are reported.

  1. At zero hopping, U = 0, kappa = 0, alpha = 1/2 and T_0 = 0.06J, the boundary in
     g_ph/J at omega = J, 2J and 5J, together with U_ph/J = g_ph^2/omega there.
  2. For rows 6 and 7 of the paper, how far the row's own g_ph lies below the boundary
     computed at that row's U, kappa, omega, alpha and T_0.

Exit status is 0 if the printed values are the ones the paper quotes, 1 otherwise.

    python probe_xover.py
"""

from __future__ import annotations

import sys
import time

sys.path.insert(0, ".")

import mpmath as mp

import replay

mp.mp.dps = 50

I_W_iv, _W_3, Jsq_half, Jsq_three_quarters, Jsq_one, I_sin_iv = replay.compute_3d_integrals()
ALPHA_MAP = {(1, 2): Jsq_half, (3, 4): Jsq_three_quarters, (1, 1): Jsq_one}


def floor(row: dict, g_ph_num: int, g_ph_den: int = 10**4):
    return replay.certify_holstein_all_T_le_T0(
        row["T0n"], row["T0d"], row["Un"], row["Ud"], row["kn"], row["kd"],
        g_ph_num, g_ph_den, row["omn"], row["omd"], row["an"], row["ad"],
        ALPHA_MAP[(row["an"], row["ad"])], I_W_iv, I_sin_iv, nr=300, ne=16,
    )


def boundary(row: dict, lo: float, hi: float, iterations: int) -> tuple[float, float]:
    if not floor(row, int(round(lo * 10**4))) > 0:
        raise SystemExit(f"starting point g_ph={lo} is not certified")
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        if floor(row, int(round(mid * 10**4))) > 0:
            lo = mid
        else:
            hi = mid
    return lo, hi


ZERO_HOPPING = {
    "T0n": 6, "T0d": 100, "Un": 0, "Ud": 1, "kn": 0, "kd": 1,
    "an": 1, "ad": 2,
}

ROWS = {
    6: dict(T0n=60, T0d=1000, Un=0, Ud=1, kn=0, kd=1, gph=0.70, omn=5, omd=1, an=1, ad=2),
    7: dict(T0n=60, T0d=1000, Un=0, Ud=1, kn=0, kd=1, gph=1.60, omn=25, omd=1, an=1, ad=2),
}

EXPECTED = {
    ("omega", 1): (0.3169, 0.1004),
    ("omega", 2): (0.4568, 0.1043),
    ("omega", 5): (0.7499, 0.1125),
    ("row", 6): 0.933,
    ("row", 7): 0.870,
}


def main() -> int:
    start = time.time()
    print("Positivity boundary of the certified floor, zero hopping, U = 0, kappa = 0,")
    print("alpha = 1/2, T_0 = 0.06J, 300 x 16 rectangle cover.")
    print()
    worst = 0.0
    for om_num in (1, 2, 5):
        row = dict(ZERO_HOPPING, gph=0.05, omn=om_num, omd=1)
        lo, hi = boundary(row, 0.05, 1.5, iterations=12)
        u_ph = lo * lo / om_num
        ref = EXPECTED[("omega", om_num)]
        worst = max(worst, abs(lo - ref[0]), abs(u_ph - ref[1]))
        print(
            f"  omega/J = {om_num:>2}: positive floor for g_ph/J <= {lo:.4f} "
            f"(next trial {hi:.4f} gives no floor), U_ph/J = g_ph^2/omega = {u_ph:.4f} "
            f"[paper: {ref[0]:.4f}, {ref[1]:.4f}]"
        )
    print()
    print("Distance of the certified rows to that boundary at their own parameters:")
    for n, row in ROWS.items():
        lo, hi = boundary(row, row["gph"], 1.4 * row["gph"], iterations=12)
        frac = row["gph"] / lo
        ref = EXPECTED[("row", n)]
        worst = max(worst, abs(frac - ref))
        print(
            f"  row {n}: g_ph/J = {row['gph']:.2f} lies at {100 * frac:.1f}% of the boundary "
            f"g_ph/J <= {lo:.4f} (omega/J = {row['omn']}, boundary U_ph/J = "
            f"{lo * lo / row['omn']:.4f}) [paper: {100 * ref:.1f}%]"
        )
    print()
    elapsed = time.time() - start
    if worst > 0.002:
        print(f"MISMATCH: largest deviation from the quoted values is {worst:.5f} (in {elapsed:.1f}s)")
        return 1
    print(f"PROBE_XOVER_PASSED: all values agree within {worst:.5f} ({elapsed:.1f}s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
