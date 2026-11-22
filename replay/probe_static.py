# pyright: reportMissingImports=false, reportOperatorIssue=false, reportArgumentType=false, reportAttributeAccessIssue=false
"""Decomposition of the certified cost c at the boundary coupling.

This script imports the unchanged replay.py that sits next to it and repeats the
interval certificate of replay.certify_holstein_all_T_le_T0 rectangle by rectangle,
recording which rectangle attains the minimum of the certified floor.  At that
rectangle (the "minimising rectangle" of Section 7.3 of the paper) it reports the
three pieces of

    c = e_hop/(6J) + B_ph/(3J) + I_W/(9 beta^2 J^2),

namely the hopping piece, the two phonon pieces of B_ph/(3J) (the static term
U_ph (1-r)/(3J) and the Cauchy-Schwarz fluctuation piece) and the temperature piece,
for the boundary coupling of Section 7.3: U = 0, kappa = 0, alpha = 1/2,
T_0 = 0.06J, g_ph/J at the largest certified value at omega/J = 1, 2 and 5.

The mirrored loop is checked against replay.certify_holstein_all_T_le_T0 itself: the
two floors must agree exactly at every parameter point, otherwise the script exits 1.

    python probe_static.py
"""

from __future__ import annotations

import sys
import time

sys.path.insert(0, ".")

import mpmath as mp
from mpmath import iv

import replay

mp.mp.dps = 50
iv.dps = 50

iv_lo = replay.iv_lo
iv_hi = replay.iv_hi
H2_iv = replay.H2_iv

I_W_iv, _W_3, Jsq_half, _Jsq_3q, _Jsq_1, I_sin_iv = replay.compute_3d_integrals()

PARAMS = {1: dict(gph_n=3169, gph_d=10**4), 2: dict(gph_n=4568, gph_d=10**4), 5: dict(gph_n=7499, gph_d=10**4)}

T0 = iv.mpf(6) / 100
EPS = iv.mpf(1) / 2


def certify_tracked(gph_n: int, gph_d: int, om_n: int, om_d: int, nr: int = 300, ne: int = 16):
    """Mirror of replay.certify_holstein_all_T_le_T0 with U = 0, kappa = 0, alpha = 1/2.

    Returns (floor, record) where record describes the minimising rectangle.
    """
    beta0 = 1 / T0
    kap = iv.mpf(0)
    gph = abs(iv.mpf(gph_n) / gph_d)
    omega = iv.mpf(om_n) / om_d
    alpha = iv.mpf(1) / 2
    eps = EPS
    Uph = (gph * gph) / omega
    Ueff = Uph

    F_up = replay.F_dimer_holstein_upper(beta0, 0, 1, gph, omega, 1, 1)
    for l_n in range(0, 6):
        for hn in range(0, 16):
            cand = replay.F_dimer_holstein_upper(beta0, hn, 20, gph, omega, l_n, 5)
            if iv_hi(cand) < iv_hi(F_up):
                F_up = cand
    beta_ent = beta0
    s_osc_iv = -iv.log(1 - iv.exp(-beta_ent * (1 - eps) * omega)) / beta_ent
    r_max_iv = 1 / (1 + iv.exp(beta0 * Ueff / 2))  # kappa = 0 and Ueff >= 0 here

    e_hop_iv = iv.mpf(0)
    comm_norm = 12 + 2 * iv.sqrt(6) * kap
    worst = mp.mpf("inf")
    record = None

    for k in range(nr):
        r0_iv = (r_max_iv * k) / nr
        r1_iv = (r_max_iv * (k + 1)) / nr
        r_iv = iv.mpf([iv_lo(r0_iv), iv_hi(r1_iv)])
        S_el = iv.log(2) + H2_iv(r_iv)

        rhs_r = F_up + s_osc_iv + S_el / beta_ent - (Ueff * r_iv) / 2 + e_hop_iv
        E_lo_raw = -rhs_r / 3
        E_hi_iv = (1 - r_iv) / iv.sqrt(3)
        if iv_lo(E_lo_raw) > iv_hi(E_hi_iv):
            continue

        E_min_iv = iv.mpf(max(mp.mpf(0), iv_lo(E_lo_raw)))
        E_max_iv = iv.mpf(max(iv_lo(E_min_iv), iv_hi(E_hi_iv)))
        E_span_iv = E_max_iv - E_min_iv
        for j in range(ne):
            e0_iv = E_min_iv + (E_span_iv * j) / ne
            e1_iv = E_min_iv + (E_span_iv * (j + 1)) / ne
            E_lo = iv.mpf(iv_lo(e0_iv))
            E_hi = iv.mpf(iv_hi(e1_iv))
            rhs_nu = rhs_r + 3 * E_hi
            nu_max_iv = iv.mpf([mp.mpf(0), max(mp.mpf(0), iv_hi(rhs_nu / (eps * omega)))])

            static_uph = Uph * (1 - r_iv)
            neg_eint_direct = static_uph + iv.sqrt(2) * gph * iv.sqrt((1 - r_iv) * nu_max_iv)
            neg_eint_comm = (
                static_uph
                + (iv.sqrt(2) * gph / omega) * comm_norm * iv.sqrt(nu_max_iv)
                + (Uph / omega) * comm_norm * iv.sqrt(1 - r_iv)
            )
            direct_wins = iv_hi(neg_eint_direct) <= iv_hi(neg_eint_comm)
            neg_eint_hi = min(iv_hi(neg_eint_direct), iv_hi(neg_eint_comm))
            err_ph_iv = iv.mpf([mp.mpf(0), neg_eint_hi]) / 3
            c_tot = e_hop_iv / 6 + err_ph_iv + I_W_iv / (9 * beta_ent * beta_ent)

            E_star_iv = Jsq_half / 4 - c_tot
            f_at_lo = E_lo - iv.sqrt(Jsq_half * (E_lo + c_tot))
            f_at_hi = E_hi - iv.sqrt(Jsq_half * (E_hi + c_tot))
            if iv_hi(E_hi) <= iv_lo(E_star_iv):
                f_lo = iv_lo(f_at_hi)
            elif iv_lo(E_lo) >= iv_hi(E_star_iv):
                f_lo = iv_lo(f_at_lo)
            else:
                f_star = -Jsq_half / 4 - c_tot
                f_lo = min(iv_lo(f_at_lo), iv_lo(f_at_hi), iv_lo(f_star))
            f_iv = iv.mpf([f_lo, max(iv_hi(f_at_lo), iv_hi(f_at_hi))])

            m2_iv = (alpha * (1 - r_iv) + f_iv) / (alpha + 1)
            if iv_lo(m2_iv) < worst:
                worst = iv_lo(m2_iv)
                record = dict(
                    k=k,
                    j=j,
                    r_lo=iv_lo(r_iv),
                    r_hi=iv_hi(r_iv),
                    static=iv_hi(static_uph) / 3,
                    direct_wins=direct_wins,
                    cs_direct=(iv_hi(neg_eint_direct) - iv_lo(static_uph)) / 3,
                    cs_comm=(iv_hi(neg_eint_comm) - iv_lo(static_uph)) / 3,
                    temp=iv_hi(I_W_iv / (9 * beta_ent * beta_ent)),
                    c=iv_hi(c_tot),
                )
    return worst, record


EXPECTED = {
    # omega/J: (r at the minimising rectangle, static piece, fluctuation piece,
    #           temperature piece, upper end of c)
    1: (0.190, 0.0271, 0.0447, 6.1e-4, 0.0724),
    2: (0.186, 0.0283, 0.0467, 6.1e-4, 0.0756),
    5: (0.178, 0.0308, 0.0502, 6.1e-4, 0.0817),
}


def main() -> int:
    start = time.time()
    print("Decomposition of c at the minimising rectangle, zero hopping, U = 0, kappa = 0,")
    print("alpha = 1/2, T_0 = 0.06J, g_ph at the largest certified value, 300 x 16 cover.")
    print()
    worst_dev = 0.0
    bad = 0
    for om in (1, 2, 5):
        p = PARAMS[om]
        floor, rec = certify_tracked(p["gph_n"], p["gph_d"], om, 1)
        ref = replay.certify_holstein_all_T_le_T0(
            6, 100, 0, 1, 0, 1, p["gph_n"], p["gph_d"], om, 1, 1, 2,
            Jsq_half, I_W_iv, I_sin_iv, nr=300, ne=16,
        )
        if floor != ref:
            print(f"  MIRROR MISMATCH at omega/J = {om}: mirrored {mp.nstr(floor, 20)}, "
                  f"replay {mp.nstr(ref, 20)}")
            bad += 1
            continue
        cs = rec["cs_direct"] if rec["direct_wins"] else rec["cs_comm"]
        form = "Cauchy-Schwarz" if rec["direct_wins"] else "commutator"
        u_ph = p["gph_n"] / p["gph_d"] * p["gph_n"] / p["gph_d"] / om
        print(
            f"  omega/J = {om}, g_ph/J = {p['gph_n'] / p['gph_d']:.4f}, U_ph/J = {u_ph:.4f}:"
            f" floor {mp.nstr(floor, 8)}"
        )
        print(
            f"    minimising rectangle k={rec['k']}, j={rec['j']}, r in "
            f"[{mp.nstr(rec['r_lo'], 6)}, {mp.nstr(rec['r_hi'], 6)}]; "
            f"static {mp.nstr(rec['static'], 4)}, {form} fluctuation {mp.nstr(cs, 4)}, "
            f"temperature {mp.nstr(rec['temp'], 2)}, total c {mp.nstr(rec['c'], 4)}"
        )
        r_ref, s_ref, f_ref, t_ref, c_ref = EXPECTED[om]
        dev = max(
            abs(float(rec["static"]) - s_ref),
            abs(float(cs) - f_ref),
            abs(float(rec["temp"]) - t_ref),
            abs(float(rec["c"]) - c_ref),
            0.0
            if float(rec["r_lo"]) <= r_ref <= float(rec["r_hi"])
            else min(abs(float(rec["r_lo"]) - r_ref), abs(float(rec["r_hi"]) - r_ref)),
        )
        worst_dev = max(worst_dev, dev)
        print(
            f"    paper: r ~ {r_ref:.3f}, static {s_ref:.4f}, fluctuation {f_ref:.4f}, "
            f"temperature {t_ref:.1e}, total c {c_ref:.4f}"
        )
    print()
    elapsed = time.time() - start
    if bad:
        print(f"MIRROR FAILED at {bad} parameter point(s) ({elapsed:.1f}s)")
        return 1
    if worst_dev > 5e-4:
        print(f"MISMATCH: largest deviation from the quoted values is {worst_dev:.5f} ({elapsed:.1f}s)")
        return 1
    print(f"PROBE_STATIC_PASSED: mirror exact, values agree within {worst_dev:.5f} ({elapsed:.1f}s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
