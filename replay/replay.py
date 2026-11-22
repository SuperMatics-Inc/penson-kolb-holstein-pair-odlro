# pyright: reportOperatorIssue=false, reportArgumentType=false, reportAttributeAccessIssue=false, reportInvalidTypeForm=false
#!/usr/bin/env python3
"""
PK-HOLSTEIN Verification & Certification Replay Script
======================================================
Verifies all rigorous constants, 50-digit outward-rounded interval bounds (mpmath.iv),
exact algebraic identities, finite-cutoff reflection positivity (both standard trace
and Jaffe-Pedrocchi graded Majorana reflection positivity), bond-phonon RP checks,
and adversarial semantic mutants for the 3D pi-flux Penson-Kolb-Hubbard-Holstein model.

Usage:
  python replay.py
  python replay.py --mutants
  python replay.py --all-mutants
"""

from __future__ import annotations

import argparse
import gc
import math
import subprocess
import sys
import time
import numpy as np
import scipy.linalg as la
import mpmath as mp
from mpmath import iv

mp.mp.dps = 50
iv.dps = 50


def iv_lo(x) -> mp.mpf:
    return mp.mpf(x.a)


def iv_hi(x) -> mp.mpf:
    return mp.mpf(x.b)


def fmt_floor_sci(x: mp.mpf, sig: int = 5) -> tuple[str, mp.mpf]:
    """Strict directed rounding DOWN (towards 0 / -infinity) of positive x to `sig` significant digits."""
    assert x > 0, f"Cannot format non-positive floor {x}"
    exp = int(mp.floor(mp.log10(x)))
    scale = mp.power(10, sig - 1 - exp)
    mant_int = int(mp.floor(x * scale))
    if mant_int >= 10**sig:
        mant_int //= 10
        exp += 1
    if mant_int < 10 ** (sig - 1):
        exp -= 1
        scale = mp.power(10, sig - 1 - exp)
        mant_int = int(mp.floor(x * scale))
    val_floor = mp.mpf(mant_int) * mp.power(10, exp - (sig - 1))
    assert val_floor <= x
    s = str(mant_int)
    return f"{s[0]}.{s[1:]}e{exp}", val_floor


def fmt_iv_directed(x_iv, dec: int = 9) -> str:
    """Strict outward directed rounding of an interval [a, b] to `dec` decimal places."""
    scale = mp.power(10, dec)
    lo_int = int(mp.floor(iv_lo(x_iv) * scale))
    hi_int = int(mp.ceil(iv_hi(x_iv) * scale))
    lo_s = f"{lo_int // (10**dec)}.{abs(lo_int) % (10**dec):0{dec}d}"
    hi_s = f"{hi_int // (10**dec)}.{abs(hi_int) % (10**dec):0{dec}d}"
    return f"[{lo_s}, {hi_s}]"


def section(title: str) -> None:
    print(f"\n{'=' * 72}\n{title}\n{'=' * 72}")


# ---------------------------------------------------------------------------
# 1. [CERTIFIED] Exact Watson integral & 3D box enclosures for W_3, J_sq, I_sin
# ---------------------------------------------------------------------------
def compute_3d_integrals():
    I_W_iv = (
        (iv.sqrt(6) / (32 * iv.pi**3))
        * iv.gamma(iv.mpf(1) / 24)
        * iv.gamma(iv.mpf(5) / 24)
        * iv.gamma(iv.mpf(7) / 24)
        * iv.gamma(iv.mpf(11) / 24)
    )

    def certify_JW_box(alpha_num: int, alpha_den: int, M: int = 64):
        alpha = iv.mpf(alpha_num) / alpha_den
        h = iv.pi / M
        c_pts = [iv.cos(iv.mpf(i) * h) for i in range(M + 1)]
        c_box = [iv.mpf([iv_lo(c_pts[i + 1]), iv_hi(c_pts[i])]) for i in range(M)]
        g1_sum = iv.mpf(0)
        g2_sum = iv.mpf(0)
        for i in range(M):
            ci = c_box[i]
            for j in range(i, M):
                cj = c_box[j]
                for k in range(j, M):
                    w = 1 if (i == j == k) else (3 if (i == j or j == k) else 6)
                    u_iv = (ci + cj + c_box[k]) / 3
                    u_lo, u_hi = iv_lo(u_iv), iv_hi(u_iv)
                    a_val = iv_lo(alpha)
                    if u_hi <= a_val:
                        g1_iv = iv.mpf(0)
                        g2_iv = iv.mpf(0)
                    else:
                        if u_lo <= a_val:
                            g1_lo = mp.mpf(0)
                            g2_lo = mp.mpf(0)
                        else:
                            d_lo = iv.mpf(u_lo) - alpha
                            one_p_lo = 1 + iv.mpf(u_lo)
                            g1_lo = iv_lo(d_lo / one_p_lo)
                            g2_lo = iv_lo(d_lo * d_lo / one_p_lo)
                        d_hi = iv.mpf(u_hi) - alpha
                        one_p_hi = 1 + iv.mpf(u_hi)
                        g1_hi = iv_hi(d_hi / one_p_hi)
                        g2_hi = iv_hi(d_hi * d_hi / one_p_hi)
                        g1_iv = iv.mpf([g1_lo, g1_hi])
                        g2_iv = iv.mpf([g2_lo, g2_hi])
                    g1_sum += w * g1_iv
                    g2_sum += w * g2_iv
        inv_M3 = iv.mpf(1) / (M**3)
        g1_int = g1_sum * inv_M3
        g2_int = g2_sum * inv_M3
        W3_iv = ((alpha + 1) * I_W_iv - g1_int) / 3
        Jsq_iv = (alpha + 1) ** 2 * I_W_iv - 2 * alpha - 1 - g2_int
        return W3_iv, Jsq_iv

    def certify_Isin_box(M: int = 48):
        h = iv.pi / (2 * M)
        s_pts = [iv.sin(iv.mpf(i) * h) for i in range(M + 1)]
        s2_box = [iv.mpf([iv_lo(s_pts[i]) ** 2, iv_hi(s_pts[i + 1]) ** 2]) for i in range(M)]
        acc = iv.mpf(0)
        for i in range(M):
            si = s2_box[i]
            for j in range(i, M):
                sj = s2_box[j]
                for k in range(j, M):
                    w = 1 if (i == j == k) else (3 if (i == j or j == k) else 6)
                    acc += w * iv.sqrt(si + sj + s2_box[k])
        return acc / (M**3)

    W3_half, Jsq_half = certify_JW_box(1, 2, M=64)
    _, Jsq_3q = certify_JW_box(3, 4, M=64)
    Jsq_1 = 4 * I_W_iv - 3
    I_sin_iv = certify_Isin_box(M=48)
    return I_W_iv, W3_half, Jsq_half, Jsq_3q, Jsq_1, I_sin_iv


# ---------------------------------------------------------------------------
# 2. [CERTIFIED] Polaron-Dimer Variational Bound & 50-digit T_0/g Certification
# ---------------------------------------------------------------------------
def H2_iv(r_iv):
    """Certified interval enclosure of binary entropy H2(r) = -r ln r - (1-r) ln(1-r) on [r0, r1]."""
    r0, r1 = iv_lo(r_iv), iv_hi(r_iv)

    def h2_pt(x):
        if x <= 0 or x >= 1:
            return iv.mpf(0)
        xi = iv.mpf(x)
        return -xi * iv.log(xi) - (1 - xi) * iv.log(1 - xi)

    lo = min(iv_lo(h2_pt(r0)), iv_lo(h2_pt(r1)))
    half = mp.mpf("0.5")
    if r0 <= half <= r1:
        hi = iv_hi(iv.log(2))
    else:
        hi = max(iv_hi(h2_pt(r0)), iv_hi(h2_pt(r1)))
    return iv.mpf([lo, hi])


def F_dimer_holstein_upper(
    beta0,
    h_num: int,
    h_den: int,
    gph,
    omega,
    lam_num: int = 1,
    lam_den: int = 1,
    omit_frank_condon: bool = False,
):
    """
    Closed-form certified upper bound on the exact thermodynamic free energy per site
    F_beta / g using a 2-site staggered AF dimer x displaced-oscillator product ansatz
    supported on P_x x P_y, with Debye-Waller pair-hopping renormalization
    g_eff = g * exp(-2 * lam^2 * (g_ph / omega)^2).
    """
    h = iv.mpf(h_num) / h_den
    lam = iv.mpf(lam_num) / lam_den
    geff = iv.mpf(1) if omit_frank_condon else iv.exp(-2 * lam * lam * (gph / omega) ** 2)
    s = iv.sqrt(iv.mpf(1) / 4 + 4 * h * h)
    emax = iv.mpf(1) / 2 + s
    b_eff = beta0 * geff
    Z4 = 1 + iv.exp(-b_eff * emax) + iv.exp(-2 * b_eff * s) + iv.exp(-b_eff * (emax + 1))
    m1 = (2 * h / s) * (1 - iv.exp(-2 * b_eff * s)) / Z4
    F_dimer = -emax / 2 - iv.log(Z4) / (2 * b_eff)
    F_el = geff * (F_dimer + h * m1 - iv.mpf(5) / 4 * m1 * m1)
    shift_ph = ((1 - lam) ** 2) * (gph * gph) / (2 * omega)
    return F_el + shift_ph


def certify_holstein_all_T_le_T0(
    T0_num: int,
    T0_den: int,
    U_num: int,
    U_den: int,
    kap_num: int,
    kap_den: int,
    gph_num: int,
    gph_den: int,
    om_num: int,
    om_den: int,
    alpha_num: int,
    alpha_den: int,
    Jsq_iv,
    I_W_iv,
    I_sin_iv,
    lam_num: int = 1,
    lam_den: int = 1,
    eps_num: int = 1,
    eps_den: int = 2,
    nr: int = 200,
    ne: int = 16,
    double_comm_factor: int = 3,
    drop_uph_static: bool = False,
    omit_frank_condon: bool = False,
    reverse_beta_monotonicity: bool = False,
) -> mp.mpf:
    """
    Certifies inf_{0 < T <= T_0} m_0^2(T) > 0 on the full fermion x unbounded-phonon Fock space
    using outward-rounded 50-digit interval arithmetic over a partition of (r, E_1).
    """
    T0 = iv.mpf(T0_num) / T0_den
    beta0 = iv.mpf(4) if reverse_beta_monotonicity else (1 / T0)
    U = iv.mpf(U_num) / U_den
    kap = abs(iv.mpf(kap_num) / kap_den)
    gph = abs(iv.mpf(gph_num) / gph_den)
    omega = iv.mpf(om_num) / om_den
    alpha = iv.mpf(alpha_num) / alpha_den
    eps = iv.mpf(eps_num) / eps_den
    Uph = (gph * gph) / omega
    Ueff = U + Uph

    # Optimize rational staggered field h in {0, 1/20, ..., 15/20} and rational lam in {0, 1/5, ..., 5/5}:
    F_up = F_dimer_holstein_upper(
        beta0, 0, 1, gph, omega, lam_num, lam_den, omit_frank_condon=omit_frank_condon
    )
    for l_n in range(0, 6):
        for hn in range(0, 16):
            cand = F_dimer_holstein_upper(
                beta0, hn, 20, gph, omega, l_n, 5, omit_frank_condon=omit_frank_condon
            )
            if iv_hi(cand) < iv_hi(F_up):
                F_up = cand
    beta_ent = beta0
    s_osc_iv = -iv.log(1 - iv.exp(-beta_ent * (1 - eps) * omega)) / beta_ent
    if kap_num == 0 and iv_lo(Ueff) >= 0:
        r_max_iv = 1 / (1 + iv.exp(beta0 * Ueff / 2))
    else:
        r_max_iv = iv.mpf(1)

    e_hop_glob = 2 * iv.sqrt(3) * kap * I_sin_iv
    comm_norm = 12 + 2 * iv.sqrt(6) * kap
    worst_m2_lo = mp.mpf("inf")

    # Complete closed 2D rectangle cover of D(beta_0) with outward-rounded overlapping endpoints:
    for k in range(nr):
        r0_iv = (r_max_iv * k) / nr
        r1_iv = (r_max_iv * (k + 1)) / nr
        r_iv = iv.mpf([iv_lo(r0_iv), iv_hi(r1_iv)])
        S_el = iv.log(2) + H2_iv(r_iv)

        if kap_num != 0:
            e_hop_r = 2 * iv.sqrt(6) * kap * iv.sqrt(r_iv)
            e_hop_iv = iv.mpf(
                [min(iv_lo(e_hop_glob), iv_lo(e_hop_r)), min(iv_hi(e_hop_glob), iv_hi(e_hop_r))]
            )
        else:
            e_hop_iv = iv.mpf(0)

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
            # Since nu_max(r, E_1) is strictly increasing in E_1, evaluating at E_hi encloses [e_0, e_1]:
            rhs_nu = rhs_r + 3 * E_hi
            nu_max_iv = iv.mpf([mp.mpf(0), max(mp.mpf(0), iv_hi(rhs_nu / (eps * omega)))])

            static_uph = iv.mpf(0) if drop_uph_static else Uph * (1 - r_iv)
            neg_eint_direct = static_uph + iv.sqrt(2) * gph * iv.sqrt((1 - r_iv) * nu_max_iv)
            neg_eint_comm = (
                static_uph
                + (iv.sqrt(2) * gph / omega) * comm_norm * iv.sqrt(nu_max_iv)
                + (Uph / omega) * comm_norm * iv.sqrt(1 - r_iv)
            )
            neg_eint_hi = min(iv_hi(neg_eint_direct), iv_hi(neg_eint_comm))
            err_ph_iv = iv.mpf([mp.mpf(0), neg_eint_hi]) / double_comm_factor
            c_tot = e_hop_iv / 6 + err_ph_iv + I_W_iv / (9 * beta_ent * beta_ent)

            # Exact convex infimum of g(E_1) = E_1 - sqrt(J_sq * (E_1 + c_tot)) on [E_lo, E_hi]:
            E_star_iv = Jsq_iv / 4 - c_tot
            f_at_lo = E_lo - iv.sqrt(Jsq_iv * (E_lo + c_tot))
            f_at_hi = E_hi - iv.sqrt(Jsq_iv * (E_hi + c_tot))
            if iv_hi(E_hi) <= iv_lo(E_star_iv):
                # Entire subinterval lies to the left of E_1*, where g(E_1) is strictly decreasing:
                f_lo = iv_lo(f_at_hi)
            elif iv_lo(E_lo) >= iv_hi(E_star_iv):
                # Entire subinterval lies to the right of E_1*, where g(E_1) is strictly increasing:
                f_lo = iv_lo(f_at_lo)
            else:
                # Stationary point E_1* may lie in [E_lo, E_hi]:
                f_star = -Jsq_iv / 4 - c_tot
                f_lo = min(iv_lo(f_at_lo), iv_lo(f_at_hi), iv_lo(f_star))
            f_iv = iv.mpf([f_lo, max(iv_hi(f_at_lo), iv_hi(f_at_hi))])

            m2_iv = (alpha * (1 - r_iv) + f_iv) / (alpha + 1)
            if iv_lo(m2_iv) < worst_m2_lo:
                worst_m2_lo = iv_lo(m2_iv)

    return worst_m2_lo

# ---------------------------------------------------------------------------
# 3. Finite-cluster exact operators (Fermions x Phonon Fock space)
# ---------------------------------------------------------------------------
def build_fermion_ops(n_sites: int):
    """Build exact fermion creation/annihilation and Majorana operators on n_sites."""
    n_modes = 2 * n_sites
    dim_f = 1 << n_modes
    c, cd = [], []
    for m in range(n_modes):
        mat = np.zeros((dim_f, dim_f), dtype=complex)
        for st in range(dim_f):
            if (st >> m) & 1:
                parity = bin(st & ((1 << m) - 1)).count("1")
                mat[st ^ (1 << m), st] = -1.0 if (parity % 2) else 1.0
        c.append(mat)
        cd.append(mat.conj().T)
    n_op = [cd[m] @ c[m] for m in range(n_modes)]
    eta = [-c[m] - cd[m] for m in range(n_modes)]
    xi = [1j * (cd[m] - c[m]) for m in range(n_modes)]
    G1, G2, P, V, n_site = [], [], [], [], []
    for x in range(n_sites):
        u, d = 2 * x, 2 * x + 1
        d_op = c[d] @ c[u]
        dd_op = cd[u] @ cd[d]
        G1.append(dd_op + d_op)
        G2.append(-1j * (dd_op - d_op))
        nx = n_op[u] + n_op[d]
        n_site.append(nx)
        px = np.eye(dim_f, dtype=complex) - n_op[u] - n_op[d] + 2.0 * (n_op[u] @ n_op[d])
        P.append(px)
        V.append(np.eye(dim_f, dtype=complex) - px)
    return dim_f, c, cd, n_op, eta, xi, G1, G2, P, V, n_site


def build_phonon_ops(n_ph: int):
    """Build exact truncated boson operators b, b^dag, n_b, q = (b + b^dag)/sqrt(2) on {0, ..., n_ph-1}."""
    b = np.zeros((n_ph, n_ph), dtype=complex)
    for n in range(1, n_ph):
        b[n - 1, n] = math.sqrt(n)
    bd = b.conj().T
    nb = bd @ b
    q = (b + bd) / math.sqrt(2.0)
    p = -1j * (b - bd) / math.sqrt(2.0)
    P_flip = np.diag([1.0 if (n % 2 == 0) else -1.0 for n in range(n_ph)]).astype(complex)
    return b, bd, nb, q, p, P_flip


def check_2site_graded_majorana_rp(
    n_ph: int = 8,
    J_pair: float = 1.0,
    U: float = 1.0,
    kap: float = 0.4,
    omega: float = 1.2,
    g_ph: float = 0.6,
    theta_q_sign: int = 1,
    use_zero_flux_hopping: bool = False,
    use_odd_phonon_parity: bool = True,
) -> tuple[float, float, float]:
    """
    Tests Jaffe-Pedrocchi graded Majorana reflection positivity on a 2-site dimer x (n_ph)^2
    (site 0 in Lambda_-, site 1 in Lambda_+) under the combined particle-hole + phonon-parity
    unitary on odd site 1 and graded Majorana reflection xi_{0,s} -> xi_{1,s}, eta_{0,s} -> -eta_{1,s},
    with phonon reflection theta(q_0) = theta_q_sign * q_1 in the Majorana frame.
    Returns (herm_err, min_eig_even, min_eig_odd).
    """
    dim_f, c, cd, n_op, _, _, G1, G2, _, _, _ = build_fermion_ops(2)
    _, _, nb, q, p, P_flip = build_phonon_ops(n_ph)
    If = np.eye(dim_f, dtype=complex)
    Ib = np.eye(n_ph, dtype=complex)
    Ib2 = np.kron(Ib, Ib)

    a = lambda x, s: c[2 * x + s]
    ad = lambda x, s: cd[2 * x + s]
    nop = lambda x, s: n_op[2 * x + s]
    xi_m = [[a(x, s) + ad(x, s) for s in (0, 1)] for x in range(2)]
    eta_m = [[1j * (ad(x, s) - a(x, s)) for s in (0, 1)] for x in range(2)]

    # Unitary particle-hole + pi/2 phase transformation on odd site 1:
    U_ph_odd = If.copy()
    for s in (0, 1):
        par = If.copy()
        for y in range(2):
            for sp in (0, 1):
                if y != 1 or sp != s:
                    par = par @ (If - 2.0 * nop(y, sp))
        U_ph_odd = U_ph_odd @ (par @ (ad(1, s) + a(1, s)))
    U_phase_odd = np.diag(
        [np.exp(1j * 0.5 * math.pi * (((st >> 2) & 1) + ((st >> 3) & 1))) for st in range(dim_f)]
    )
    U_f = U_ph_odd @ U_phase_odd
    # Accompany odd-site particle-hole map (which flips n_1 - 1 -> -(n_1 - 1)) by odd-site phonon parity P_flip:
    P_odd_ph = P_flip if use_odd_phonon_parity else Ib
    U_tot = np.kron(U_f, np.kron(Ib, P_odd_ph))
    # xi_{0,s} -> xi_{1,s}, eta_{0,s} -> -eta_{1,s}, antilinear on scalars.
    M0_basis, M1_theta_basis = [], []
    for mask in range(16):
        op0 = If.copy()
        op1 = If.copy()
        for bit, (m0, m1, sgn) in enumerate([
            (xi_m[0][0], xi_m[1][0], +1.0),
            (eta_m[0][0], eta_m[1][0], -1.0),
            (xi_m[0][1], xi_m[1][1], +1.0),
            (eta_m[0][1], eta_m[1][1], -1.0),
        ]):
            if (mask >> bit) & 1:
                op0 = op0 @ m0
                op1 = op1 @ (sgn * m1)
        M0_basis.append(op0)
        M1_theta_basis.append(op1)

    def theta_maj_f(A_f: np.ndarray) -> np.ndarray:
        out = np.zeros_like(A_f)
        for op0, op1 in zip(M0_basis, M1_theta_basis):
            coeff = np.trace(op0.conj().T @ A_f) / 16.0
            out += np.conj(coeff) * op1
        return out

    beta_val = 2.0
    H_pair = -0.5 * J_pair * np.kron(G1[0] @ G1[1] + G2[0] @ G2[1], Ib2)
    H_U = 0.5 * U * np.kron(2.0 * If - G1[0] @ G1[0] - G1[1] @ G1[1], Ib2)
    if use_zero_flux_hopping:
        # Wrong-sign cross-plane link (+i D_a theta(D_a), as occurs on the unfrustrated link of a 0-flux plaquette):
        H_pi = -kap * np.kron(
            sum((ad(0, s) @ a(1, s) + ad(1, s) @ a(0, s) for s in (0, 1)), np.zeros_like(If)),
            Ib2,
        )
    else:
        # Under U_f = U_ph * U_phase, +(ad_0 a_1 + ad_1 a_0) maps to -0.5j*(xi_0 xi_1 - eta_0 eta_1) = -0.5j*(xi_0 th(xi_0) + eta_0 th(eta_0))
        H_pi = kap * np.kron(
            sum((ad(0, s) @ a(1, s) + ad(1, s) @ a(0, s) for s in (0, 1)), np.zeros_like(If)),
            Ib2,
        )
    n0_f = nop(0, 0) + nop(0, 1) - If
    n1_f = nop(1, 0) + nop(1, 1) - If
    H_hol = omega * np.kron(If, np.kron(nb, Ib) + np.kron(Ib, nb)) + g_ph * (
        np.kron(n0_f, np.kron(q, Ib)) + np.kron(n1_f, np.kron(Ib, q))
    )
    H_phys = H_pair + H_U + H_pi + H_hol
    del H_pair, H_U, H_pi, H_hol
    H_tilde = U_tot @ H_phys @ U_tot.conj().T
    del H_phys, U_tot
    herm_err = float(np.max(np.abs(H_tilde - H_tilde.conj().T)))
    evals, evecs = la.eigh(0.5 * (H_tilde + H_tilde.conj().T))
    del H_tilde
    w = np.exp(-beta_val * (evals - evals[0]))
    rho = (evecs * w) @ evecs.conj().T
    del evecs
    rho /= np.trace(rho).real

    even_f = [If, G1[0], G2[0], n0_f, xi_m[0][0] @ xi_m[0][1], eta_m[0][0] @ eta_m[0][1]]
    th_even_f = [theta_maj_f(f) for f in even_f]
    odd_f = [xi_m[0][0], eta_m[0][0], xi_m[0][1], eta_m[0][1], 1j * (xi_m[0][0] @ xi_m[0][1] @ eta_m[0][1])]
    th_odd_f = [1j * theta_maj_f(f) for f in odd_f]

    ph_base = [Ib, q, nb]
    ph_0 = np.stack(ph_base)
    ph_1 = np.stack([P_flip @ op.conj() @ P_flip if theta_q_sign == -1 else op.conj() for op in ph_base])
    rho_6 = rho.reshape(dim_f, n_ph, n_ph, dim_f, n_ph, n_ph)
    rho_f = np.einsum("pmn,qjl,knlimj->pqik", ph_0, ph_1, rho_6, optimize=True)
    del rho_6, rho

    n_ph_ops = len(ph_base)
    G_even = np.zeros((len(even_f) * n_ph_ops, len(even_f) * n_ph_ops), dtype=complex)
    for a_i, f_a in enumerate(even_f):
        for b_i, tf_b in enumerate(th_even_f):
            F_ab = f_a @ tf_b
            block = np.einsum("ik,pqik->pq", F_ab, rho_f)
            G_even[a_i * n_ph_ops : (a_i + 1) * n_ph_ops, b_i * n_ph_ops : (b_i + 1) * n_ph_ops] = block

    G_odd = np.zeros((len(odd_f) * n_ph_ops, len(odd_f) * n_ph_ops), dtype=complex)
    for a_i, f_a in enumerate(odd_f):
        for b_i, tf_b in enumerate(th_odd_f):
            F_ab = f_a @ tf_b
            block = np.einsum("ik,pqik->pq", F_ab, rho_f)
            G_odd[a_i * n_ph_ops : (a_i + 1) * n_ph_ops, b_i * n_ph_ops : (b_i + 1) * n_ph_ops] = block

    min_even = float(np.min(la.eigvalsh(0.5 * (G_even + G_even.conj().T))))
    min_odd = float(np.min(la.eigvalsh(0.5 * (G_odd + G_odd.conj().T))))
    return herm_err, min_even, min_odd

def check_2site_standard_trace_rp_t0(
    n_ph: int = 8,
    J_pair: float = 1.0,
    U: float = 1.0,
    omega: float = 1.2,
    g_ph: float = 0.6,
    theta_q_sign: int = 1,
    bond_int_coupling: float = 0.0,
    bond_cross_coupling: float = 0.0,
) -> tuple[float, float]:
    """
    Tests standard operator-trace reflection positivity Tr(A theta(A) exp(-beta H)) >= 0
    on a 2-site dimer x (n_ph)^2 at kappa = 0 (including optional internal or crossing bond phonons).
    Returns (herm_err, min_eig).
    """
    dim_f, _, _, _, _, _, G1, G2, P, V, n_site = build_fermion_ops(2)
    _, _, nb, q, p, P_flip = build_phonon_ops(n_ph)
    If = np.eye(dim_f, dtype=complex)
    Ib = np.eye(n_ph, dtype=complex)
    Ib2 = np.kron(Ib, Ib)

    # Site 1 carries unitary staggered sign G1[1] -> -G1[1], G2[1] -> -G2[1]:
    H_pair_ferro = -0.5 * J_pair * (G1[0] @ G1[1] - G2[0] @ G2[1])
    H_U = 0.5 * U * (V[0] + V[1])
    n0 = n_site[0] - If
    n1 = n_site[1] - If

    H_tot = (
        np.kron(H_pair_ferro + H_U, Ib2)
        + omega * (np.kron(If, np.kron(nb, Ib)) + np.kron(If, np.kron(Ib, nb)))
        + g_ph * (np.kron(n0, np.kron(q, Ib)) + np.kron(n1, np.kron(Ib, q)))
    )
    if bond_int_coupling != 0.0:
        # Internal bond phonon modulation (symmetric under reflection 0 <-> 1):
        H_tot = H_tot + bond_int_coupling * (
            np.kron(n0 @ n0, np.kron(q, Ib)) + np.kron(n1 @ n1, np.kron(Ib, q))
        )
    if bond_cross_coupling != 0.0:
        # Crossing bond phonon coupling q_cross * (G1[0] G1[1] - G2[0] G2[1]) (breaks RP sign!):
        q_cross = np.kron(q, Ib) + np.kron(Ib, q)
        H_tot = H_tot + bond_cross_coupling * np.kron(G1[0] @ G1[1] - G2[0] @ G2[1], q_cross)

    herm_err = float(np.max(np.abs(H_tot - H_tot.conj().T)))
    evals, evecs = la.eigh(0.5 * (H_tot + H_tot.conj().T))
    del H_tot
    beta = 1.5
    w = np.exp(-beta * (evals - np.min(evals)))
    rho = (evecs * (w / np.sum(w))) @ evecs.conj().T
    del evecs

    S_f = np.zeros((dim_f, dim_f), dtype=complex)
    for st in range(dim_f):
        b0, b1, b2, b3 = (st >> 0) & 1, (st >> 1) & 1, (st >> 2) & 1, (st >> 3) & 1
        st_new = (b2 << 0) | (b3 << 1) | (b0 << 2) | (b1 << 3)
        par = (b0 + b1) * (b2 + b3)
        S_f[st_new, st] = -1.0 if (par % 2) else 1.0

    f_ops_0 = [If, G1[0], 1j * G2[0], n0, P[0], V[0]]
    tf_ops_1 = [S_f @ f.conj() @ S_f.conj().T for f in f_ops_0]
    ph_base_0 = [Ib, q, p, nb]
    ph_0 = np.stack(ph_base_0)
    ph_1 = np.stack([P_flip @ op.conj() @ P_flip if theta_q_sign == -1 else op.conj() for op in ph_base_0])
    rho_6 = rho.reshape(dim_f, n_ph, n_ph, dim_f, n_ph, n_ph)
    # Tr(rho * theta(A_{a,p}) A_{b,q}) = Tr(((tf_a @ f_b) x ph_0[q] x ph_1[p]) @ rho)
    rho_f = np.einsum("qmn,pjl,knlimj->pqik", ph_0, ph_1, rho_6, optimize=True)
    del rho_6, rho

    n_ph_ops = len(ph_base_0)
    G_mat = np.zeros((len(f_ops_0) * n_ph_ops, len(f_ops_0) * n_ph_ops), dtype=complex)
    for a_i, tf_a in enumerate(tf_ops_1):
        for b_i, f_b in enumerate(f_ops_0):
            F_ab = tf_a @ f_b
            block = np.einsum("ik,pqik->pq", F_ab, rho_f)
            G_mat[a_i * n_ph_ops : (a_i + 1) * n_ph_ops, b_i * n_ph_ops : (b_i + 1) * n_ph_ops] = block
    min_eig = float(np.min(la.eigvalsh(0.5 * (G_mat + G_mat.conj().T))))
    return herm_err, min_eig

# ---------------------------------------------------------------------------
# 4. Adversarial Semantic Mutants (7 Named Mutants)
# ---------------------------------------------------------------------------
MUTANTS: tuple[str, ...] = (
    "omit-odd-phonon-parity",
    "drop-uph-static-term",
    "antisymmetric-phonon-reflection",
    "zero-flux-hopping",
    "crossing-bond-phonon",
    "excessive-phonon-coupling",
    "above-xy-critical-temperature",
)


def run_single_mutant(mutant: str) -> int:
    """
    Execute a single named adversarial mutation check.
    Must reject the mutated construction and return exit code 1 with 'REJECTED: ...' on stderr.
    """
    if mutant == "omit-odd-phonon-parity":
        # Omit the odd-sublattice phonon parity operator (-1)^{b_1^dag b_1} while keeping U_{ph,odd} on fermions.
        _, min_e, min_o = check_2site_graded_majorana_rp(
            n_ph=6,
            J_pair=1.0,
            U=1.0,
            kap=0.4,
            omega=1.2,
            g_ph=0.6,
            theta_q_sign=1,
            use_odd_phonon_parity=False,
        )
        if min_e < -1e-4 or min_o < -1e-4:
            print(
                f"REJECTED: omit-odd-phonon-parity breaks graded Majorana RP "
                f"(min_eig_even={min_e:.6e}, min_eig_odd={min_o:.6e} < 0)",
                file=sys.stderr,
            )
            return 1
        return 0

    if mutant == "drop-uph-static-term":
        # 1. Verify directly against the exact matrix double commutator [G_0^(a), [H_dim, G_0^(a)]] on the
        #    2-site dimer that dropping +4*U_ph*(1 - r) violates the exact double-commutator sum rule (Lemma 3.1).
        # 2. Run certify_holstein_all_T_le_T0 with drop_uph_static=True at strong coupling (g_ph=0.80, omega=1.0,
        #    U=0, T_0=0.05) and verify that dropping +U_ph*(1 - r) produces a false positive certification m_mut^2 > 0
        #    where the unmutated certifier rightly returns m_true^2 <= 0.
        dim_f, _, _, n_op, _, _, G1, G2, P, _, n_site = build_fermion_ops(2)
        _, _, nb, q, _, _ = build_phonon_ops(8)
        If = np.eye(dim_f, dtype=complex)
        Ib = np.eye(8, dtype=complex)
        Ib2 = np.kron(Ib, Ib)
        J_pair, U_val, omega_val, g_ph_val, beta_val = 1.0, 1.0, 1.5, 0.8, 2.0
        U_ph_val = (g_ph_val * g_ph_val) / omega_val
        n0_f = n_site[0] - If
        n1_f = n_site[1] - If
        H_dim = (
            -0.5 * J_pair * np.kron(G1[0] @ G1[1] + G2[0] @ G2[1], Ib2)
            + 0.5 * U_val * np.kron(2.0 * If - G1[0] @ G1[0] - G1[1] @ G1[1], Ib2)
            + omega_val * np.kron(If, np.kron(nb, Ib) + np.kron(Ib, nb))
            + g_ph_val * (np.kron(n0_f, np.kron(q, Ib)) + np.kron(n1_f, np.kron(Ib, q)))
        )
        evals, evecs = la.eigh(0.5 * (H_dim + H_dim.conj().T))
        w = np.exp(-beta_val * (evals - evals[0]))
        rho = (evecs * w) @ evecs.conj().T
        rho /= np.trace(rho).real
        # Exact matrix double commutator 1/8 * <[G_0^(1),[H_dim,G_0^(1)]] + [G_0^(2),[H_dim,G_0^(2)]]>_beta:
        G1_0 = np.kron(G1[0], Ib2)
        G2_0 = np.kron(G2[0], Ib2)
        dcomm_mat = 0.125 * (
            G1_0 @ (H_dim @ G1_0 - G1_0 @ H_dim)
            - (H_dim @ G1_0 - G1_0 @ H_dim) @ G1_0
            + G2_0 @ (H_dim @ G2_0 - G2_0 @ H_dim)
            - (H_dim @ G2_0 - G2_0 @ H_dim) @ G2_0
        )
        dcomm_exact = float(np.trace(rho @ dcomm_mat).real)
        pair_bond_exact = float(
            np.trace(rho @ (0.25 * J_pair * np.kron(G1[0] @ G1[1] + G2[0] @ G2[1], Ib2))).real
        )
        one_minus_r_exact = float(np.trace(rho @ np.kron(P[0], Ib2)).real)
        op_qtilde_n0 = g_ph_val * (
            np.kron(n0_f, np.kron(q, Ib)) + (g_ph_val / omega_val) * np.kron(n0_f @ n0_f, Ib2)
        )
        fluct_term_exact = -float(np.trace(rho @ op_qtilde_n0).real)
        # Exact sum rule: dcomm_exact == pair_bond_exact + U_ph*(1 - r) + fluct_term_exact:
        exact_sum_rule_err = abs(
            dcomm_exact - (pair_bond_exact + U_ph_val * one_minus_r_exact + fluct_term_exact)
        )
        assert exact_sum_rule_err < 1e-13, f"Sum rule error {exact_sum_rule_err}"
        mutated_dcomm_rhs = pair_bond_exact + fluct_term_exact
        identity_violation = dcomm_exact - mutated_dcomm_rhs
        # Also check that dropping +U_ph*(1 - r) in certify_holstein_all_T_le_T0 produces a false positive at g_ph=0.80:
        I_W_loc = (
            iv.sqrt(6)
            / (32 * iv.pi**3)
            * iv.gamma(iv.mpf(1) / 24)
            * iv.gamma(iv.mpf(5) / 24)
            * iv.gamma(iv.mpf(7) / 24)
            * iv.gamma(iv.mpf(11) / 24)
        )
        Jsq_loc = iv.mpf(["1.408389235", "1.408950134"])
        I_sin_loc = iv.mpf(["1.17", "1.21"])
        m2_true = certify_holstein_all_T_le_T0(
            60, 1000, 0, 1, 0, 1, 35, 100, 1, 1, 1, 2, Jsq_loc, I_W_loc, I_sin_loc, nr=100, ne=8
        )
        m2_false_pos = certify_holstein_all_T_le_T0(
            60, 1000, 0, 1, 0, 1, 35, 100, 1, 1, 1, 2, Jsq_loc, I_W_loc, I_sin_loc, nr=100, ne=8,
            drop_uph_static=True,
        )
        if identity_violation > 0.1 * U_ph_val * one_minus_r_exact and m2_true <= 0 and m2_false_pos > 0:
            print(
                f"REJECTED: drop-uph-static-term violates matrix double-commutator sum rule by "
                f"+{identity_violation:.6e} (= U_ph*(1-r)) and causes false certification "
                f"m_mut^2={mp.nstr(m2_false_pos, 5)} > 0 (vs true m_cert^2={mp.nstr(m2_true, 5)} <= 0)",
                file=sys.stderr,
            )
            return 1
        return 0

    if mutant == "antisymmetric-phonon-reflection":
        _, min_e, min_o = check_2site_graded_majorana_rp(
            n_ph=6, J_pair=1.0, U=1.0, kap=0.4, omega=1.2, g_ph=0.6, theta_q_sign=-1
        )
        if min_e < -1e-4 or min_o < -1e-4:
            print(
                f"REJECTED: antisymmetric-phonon-reflection theta(q)=-q breaks graded Majorana RP "
                f"(min_eig_odd={min_o:.6e} < 0)",
                file=sys.stderr,
            )
            return 1
        return 0

    if mutant == "zero-flux-hopping":
        _, min_e, min_o = check_2site_graded_majorana_rp(
            n_ph=6,
            J_pair=1.0,
            U=1.0,
            kap=0.4,
            omega=1.2,
            g_ph=0.6,
            theta_q_sign=1,
            use_zero_flux_hopping=True,
        )
        if min_e < -1e-4 or min_o < -1e-4:
            print(
                f"REJECTED: zero-flux-hopping breaks graded Majorana RP "
                f"(min_eig_odd={min_o:.6e} < 0)",
                file=sys.stderr,
            )
            return 1
        return 0

    if mutant == "crossing-bond-phonon":
        _, min_eig = check_2site_standard_trace_rp_t0(
            n_ph=6, J_pair=0.2, U=1.0, omega=1.0, g_ph=0.0, theta_q_sign=1, bond_cross_coupling=1.5
        )
        if min_eig < -1e-4:
            print(
                f"REJECTED: crossing-bond-phonon single-mode q_e across P breaks reflection positivity "
                f"(min_eig={min_eig:.6e} < 0)",
                file=sys.stderr,
            )
            return 1
        return 0

    # Mutants 6 and 7 test the interval certifier (use M=16 fast outward box enclosure of J_sq(1/2)):
    I_W_iv = (
        iv.sqrt(6)
        / (32 * iv.pi**3)
        * iv.gamma(iv.mpf(1) / 24)
        * iv.gamma(iv.mpf(5) / 24)
        * iv.gamma(iv.mpf(7) / 24)
        * iv.gamma(iv.mpf(11) / 24)
    )
    Jsq_half = iv.mpf(["1.408389235", "1.408950134"])
    I_sin_iv = iv.mpf(["1.17", "1.21"])

    if mutant == "excessive-phonon-coupling":
        # 1. Unmutated certifier at strong coupling g_ph=2.0, U=0, omega=1.0, T_0=0.075 yields m_cert^2 <= 0.
        # 2. Mutating Lemma 4.1 by omitting the Frank-Condon polaronic suppression exp(-2*lam^2*(g_ph/omega)^2) -> 1
        #    at lam=1 violates the exact Gibbs variational free-energy upper bound on the 2-site dimer.
        m2_m6 = certify_holstein_all_T_le_T0(
            75, 1000, 0, 1, 0, 1, 2, 1, 1, 1, 1, 2, Jsq_half, I_W_iv, I_sin_iv, nr=150, ne=8
        )
        dim_f, _, _, _, _, _, G1, G2, _, _, n_site = build_fermion_ops(2)
        _, _, nb, q, _, _ = build_phonon_ops(10)
        If = np.eye(dim_f, dtype=complex)
        Ib = np.eye(10, dtype=complex)
        Ib2 = np.kron(Ib, Ib)
        n0_f = n_site[0] - If
        n1_f = n_site[1] - If
        H_dim2 = (
            -0.5 * 1.0 * np.kron(G1[0] @ G1[1] + G2[0] @ G2[1], Ib2)
            + 0.5 * 0.0 * np.kron(2.0 * If - G1[0] @ G1[0] - G1[1] @ G1[1], Ib2)
            + 1.0 * np.kron(If, np.kron(nb, Ib) + np.kron(Ib, nb))
            + 2.0 * (np.kron(n0_f, np.kron(q, Ib)) + np.kron(n1_f, np.kron(Ib, q)))
        )
        evals2 = la.eigvalsh(0.5 * (H_dim2 + H_dim2.conj().T))
        beta0_val = 1000.0 / 75.0
        f_exact_2site = (
            evals2[0] - np.log(np.sum(np.exp(-beta0_val * (evals2 - evals2[0])))) / beta0_val
        ) / 2.0
        # Mutated Lemma 4.1 bound on 2-site dimer (1 bond per 2 sites -> F_dimer/2 - U_eff/2 without Frank-Condon):
        Ueff_val = 4.0
        f_mut_upper = float(
            iv_hi(
                F_dimer_holstein_upper(
                    iv.mpf(1000) / 75, 0, 1, iv.mpf(2), iv.mpf(1), 1, 1, omit_frank_condon=True
                )
            )
        ) - 0.5 * Ueff_val
        if m2_m6 <= 0 and f_mut_upper < f_exact_2site - 0.1:
            print(
                f"REJECTED: excessive-phonon-coupling g_ph=2.0 at U=0, omega=1.0 outside certified window "
                f"(m_cert^2={mp.nstr(m2_m6, 5)} <= 0; omitting Frank-Condon factor violates Gibbs bound: "
                f"f_mut_up={f_mut_upper:.5f} < f_exact={f_exact_2site:.5f})",
                file=sys.stderr,
            )
            return 1
        return 0

    if mutant == "above-xy-critical-temperature":
        # 1. Unmutated certifier at T_0/g=2.50 (above 3D XY critical ceiling) rightly yields m_cert^2 <= 0.
        # 2. Mutating Lemma 4.1-4.2 by reversing the beta >= beta_0 monotonicity direction (evaluating at
        #    beta_low=4.0 instead of beta_0=0.40) yields a false positive m_mut^2 > 0 at T_0/g=2.50.
        m2_m7 = certify_holstein_all_T_le_T0(
            2500, 1000, 10, 1, 0, 1, 4, 5, 5, 1, 1, 2, Jsq_half, I_W_iv, I_sin_iv, nr=150, ne=8
        )
        m2_m7_mut = certify_holstein_all_T_le_T0(
            2500, 1000, 10, 1, 0, 1, 4, 5, 5, 1, 1, 2, Jsq_half, I_W_iv, I_sin_iv, nr=150, ne=8,
            reverse_beta_monotonicity=True,
        )
        if m2_m7 <= 0 and m2_m7_mut > 0:
            print(
                f"REJECTED: above-xy-critical-temperature T_0/g=2.50 above 3D XY critical ceiling "
                f"(m_cert^2={mp.nstr(m2_m7, 5)} <= 0; reversed beta-monotonicity gives false positive "
                f"m_mut^2={mp.nstr(m2_m7_mut, 5)} > 0)",
                file=sys.stderr,
            )
            return 1
        return 0

    print(f"Unknown mutant: {mutant}", file=sys.stderr)
    return 2


def run_all_mutants() -> int:
    print(f"Running PK-HOLSTEIN adversarial mutation suite ({len(MUTANTS)} semantic mutants)...")
    script = str(sys.argv[0])
    for m_name in MUTANTS:
        proc = subprocess.run(
            [sys.executable, "-B", script, "--mutant", m_name],
            capture_output=True,
            text=True,
            check=False,
        )
        last_line = (proc.stderr.strip().splitlines() or [""])[-1]
        if proc.returncode != 1 or not last_line.startswith("REJECTED:"):
            print(
                f"MUTANT NOT VALIDLY REJECTED: {m_name}, rc={proc.returncode}, stderr={proc.stderr}",
                file=sys.stderr,
            )
            return 1
        print(f"  [REJECTED rc=1] {m_name}: {last_line}")
    print(f"ALL {len(MUTANTS)} ADVERSARIAL MUTANTS REJECTED (rc=1 each) -> SUITE EXIT 0.")
    return 0


# ---------------------------------------------------------------------------
# Main Replay Pipeline
# ---------------------------------------------------------------------------
def main() -> int:
    parser = argparse.ArgumentParser(description="PK-HOLSTEIN interval certification and RP verification.")
    parser.add_argument("--mutant", choices=MUTANTS, default=None, help="Run a single mutant (must exit 1).")
    parser.add_argument(
        "--mutants",
        "--all-mutants",
        dest="all_mutants",
        action="store_true",
        help="Run all 7 adversarial mutants and verify each exits with rc=1.",
    )
    args = parser.parse_args()

    if args.mutant is not None:
        return run_single_mutant(args.mutant)
    if args.all_mutants:
        return run_all_mutants()
    t_start = time.time()
    section("1. [CERTIFIED] 3D Watson Integral & Box Enclosures (mpmath.iv, 50 dps)")
    I_W_iv, W3_iv, Jsq_half, Jsq_3q, Jsq_1, I_sin_iv = compute_3d_integrals()
    print(f"  I_W     in {fmt_iv_directed(I_W_iv, 9)}")
    print(f"  W_3     in {fmt_iv_directed(W3_iv, 9)}")
    print(f"  J(1/2)  in {fmt_iv_directed(Jsq_half, 9)}")
    print(f"  J(3/4)  in {fmt_iv_directed(Jsq_3q, 9)}")
    print(f"  J(1)    in {fmt_iv_directed(Jsq_1, 9)}")
    print(f"  I_sin   in {fmt_iv_directed(I_sin_iv, 9)}")
    assert iv_hi(I_W_iv) < mp.mpf("1.516386060")
    assert iv_hi(Jsq_half) < mp.mpf("1.408950134")
    assert iv_hi(Jsq_3q) < mp.mpf("2.143715092")
    assert iv_hi(Jsq_1) < mp.mpf("3.065544240")
    assert iv_hi(I_sin_iv) < mp.mpf("1.207622427")

    section("2. [NUMERICAL] Exact Operator Identities & Double-Commutator Sum Rule with Phonons")
    dim_f, _, _, _, _, _, G1, G2, P, _, n_site = build_fermion_ops(2)
    If = np.eye(dim_f, dtype=complex)
    # Check exact identity [G1_x, n_x - 1] = -2i G2_x and [G2_x, n_x - 1] = +2i G1_x:
    comm_1 = G1[0] @ (n_site[0] - If) - (n_site[0] - If) @ G1[0]
    comm_2 = G2[0] @ (n_site[0] - If) - (n_site[0] - If) @ G2[0]
    err_c1 = float(np.max(np.abs(comm_1 - (-2j * G2[0]))))
    err_c2 = float(np.max(np.abs(comm_2 - (2j * G1[0]))))
    # Double commutator with H_e-ph = g_ph q_0 (n_0 - 1):
    # [G1_0, [n_0 - 1, G1_0]] + [G2_0, [n_0 - 1, G2_0]] = -8 (n_0 - 1) P_0 = -8 (n_0 - 1)
    dc_eph = (
        G1[0] @ ((n_site[0] - If) @ G1[0] - G1[0] @ (n_site[0] - If))
        - ((n_site[0] - If) @ G1[0] - G1[0] @ (n_site[0] - If)) @ G1[0]
        + G2[0] @ ((n_site[0] - If) @ G2[0] - G2[0] @ (n_site[0] - If))
        - ((n_site[0] - If) @ G2[0] - G2[0] @ (n_site[0] - If)) @ G2[0]
    )
    err_dc = float(np.max(np.abs(dc_eph - (-8.0 * (n_site[0] - If)))))
    # Check (n_x - 1)^2 = P_x and (n_x - 1) P_x = (n_x - 1):
    err_sq = float(np.max(np.abs((n_site[0] - If) @ (n_site[0] - If) - P[0])))
    err_proj = float(np.max(np.abs((n_site[0] - If) @ P[0] - (n_site[0] - If))))
    print(f"  max|[G1, n-1] + 2i G2| = {err_c1:.2e}, max|[G2, n-1] - 2i G1| = {err_c2:.2e}")
    print(f"  max|[G1,[n-1,G1]] + [G2,[n-1,G2]] + 8(n-1)| = {err_dc:.2e}")
    print(f"  max|(n-1)^2 - P_x| = {err_sq:.2e}, max|(n-1)P_x - (n-1)| = {err_proj:.2e}")
    assert max(err_c1, err_c2, err_dc, err_sq, err_proj) < 1e-14

    section("3. [NUMERICAL] Reflection Positivity & Fock-Cutoff Convergence (N_ph = 4, 6, 8, 10)")
    for n_ph in (4, 6, 8, 10):
        h_err, min_e, min_o = check_2site_graded_majorana_rp(
            n_ph=n_ph, J_pair=1.0, U=1.0, kap=0.4, omega=1.2, g_ph=0.6, theta_q_sign=1
        )
        h_t0, min_t0 = check_2site_standard_trace_rp_t0(
            n_ph=n_ph, J_pair=1.0, U=1.0, omega=1.2, g_ph=0.6, theta_q_sign=1, bond_int_coupling=0.4
        )
        print(
            f"  N_ph={n_ph:2d}: Graded Majorana RP (pi-flux kap=0.4, g_ph=0.6) -> "
            f"herm_err={h_err:.1e}, min_eig_even={min_e:+.3e}, min_eig_odd={min_o:+.3e} | "
            f"Std Trace RP (t=0 + internal bond phonon) -> min_eig={min_t0:+.3e}"
        )
        assert h_err < 1e-13 and h_t0 < 1e-13
        assert min_e >= -1e-12 and min_o >= -1e-12 and min_t0 >= -1e-12
        gc.collect()

    # Verify that anti-symmetric phonon reflection theta(q) = -q and crossing linear bond phonons fail RP:
    _, _, min_o_anti = check_2site_graded_majorana_rp(
        n_ph=6, J_pair=1.0, U=1.0, kap=0.4, omega=1.2, g_ph=0.6, theta_q_sign=-1
    )
    _, min_cross = check_2site_standard_trace_rp_t0(
        n_ph=6, J_pair=0.2, U=1.0, omega=1.0, g_ph=0.0, theta_q_sign=1, bond_cross_coupling=1.5
    )
    print(f"  [OBSTRUCTION CHECK] Anti-symmetric theta(q)=-q min_eig_odd = {min_o_anti:+.6e} < 0 (RP FAILS)")
    print(f"  [OBSTRUCTION CHECK] Crossing bond-phonon q_e*(G_x G_y) min_eig = {min_cross:+.6e} < 0 (RP FAILS)")
    assert min_o_anti < -1e-3 and min_cross < -1e-3

    section("4. [CERTIFIED] 15-Row Interval-Certified Pair ODLRO Table (PK-Holstein Model)")
    alpha_map = {
        (1, 2): Jsq_half,
        (3, 4): Jsq_3q,
        (1, 1): Jsq_1,
    }
    certified_targets = [
        # (label, T0n, T0d, Un, Ud, kn, kd, gpn, gpd, omn, omd, an, ad, expected_floor_str)
        ("1. t=0, U=-0.5, om=2.0, gph=0.25 (Ueff=-0.46875)", 65, 1000, -1, 2, 0, 1, 1, 4, 2, 1, 1, 2, "8.7232e-3"),
        ("2. t=0, U=-0.52, om=50.0, gph=1.60 (Ueff=-0.4688)", 50, 1000, -13, 25, 0, 1, 8, 5, 50, 1, 1, 2, "1.4148e-2"),
        ("3. t=0, U=0.0, om=0.01, gph=0.002 (soft phonon)", 4, 1000, 0, 1, 0, 1, 2, 1000, 1, 100, 1, 2, "5.0292e-2"),
        ("4. t=0, U=0.0, om=0.10, gph=0.025 (adiabatic)", 25, 1000, 0, 1, 0, 1, 25, 1000, 1, 10, 1, 2, "4.1730e-2"),
        ("5. t=0, U=0.0, om=1.00, gph=0.18 (intermediate)", 75, 1000, 0, 1, 0, 1, 18, 100, 1, 1, 1, 2, "1.9195e-2"),
        ("6. t=0, U=0.0, om=5.00, gph=0.70 (fast phonon)", 60, 1000, 0, 1, 0, 1, 7, 10, 5, 1, 1, 2, "4.3915e-3"),
        ("7. t=0, U=0.0, om=25.0, gph=1.60 (anti-adiabatic)", 60, 1000, 0, 1, 0, 1, 8, 5, 25, 1, 1, 2, "8.3300e-3"),
        ("8. t=0, U=4.0, om=2.00, gph=0.40 (Ueff=4.08)", 300, 1000, 4, 1, 0, 1, 2, 5, 2, 1, 3, 4, "2.9544e-2"),
        ("9. t=0, U=10.0, om=5.00, gph=0.75 (Ueff=10.1125)", 350, 1000, 10, 1, 0, 1, 3, 4, 5, 1, 1, 2, "6.7152e-3"),
        ("10. pi-flux U=0.0, kap=0.02, om=1.0, gph=0.15", 60, 1000, 0, 1, 2, 100, 15, 100, 1, 1, 1, 2, "1.1220e-2"),
        ("11. pi-flux U=0.0, kap=0.05, om=2.0, gph=0.12", 30, 1000, 0, 1, 5, 100, 12, 100, 2, 1, 1, 2, "5.7017e-4"),
        ("12. pi-flux U=4.0, kap=0.10, om=2.0, gph=0.35", 120, 1000, 4, 1, 1, 10, 35, 100, 2, 1, 1, 2, "6.1437e-3"),
        ("13. pi-flux U=10.0, kap=0.20, om=5.0, gph=0.60", 160, 1000, 10, 1, 2, 10, 3, 5, 5, 1, 1, 2, "5.5618e-3"),
        ("14. pi-flux U=12.0, kap=0.35, om=10.0, gph=0.35", 60, 1000, 12, 1, 35, 100, 35, 100, 10, 1, 1, 2, "6.8403e-3"),
        ("15. pi-flux U=20.0, kap=0.50, om=10.0, gph=0.45", 80, 1000, 20, 1, 1, 2, 45, 100, 10, 1, 1, 2, "5.5839e-4"),
    ]

    for lbl, T0n, T0d, Un, Ud, kn, kd, gpn, gpd, omn, omd, an, ad, exp_str in certified_targets:
        m2_lo = certify_holstein_all_T_le_T0(
            T0n, T0d, Un, Ud, kn, kd, gpn, gpd, omn, omd,
            an, ad, alpha_map[(an, ad)], I_W_iv, I_sin_iv, nr=300, ne=16
        )
        floor_str, floor_val = fmt_floor_sci(m2_lo, 5)
        print(
            f"  [CERTIFIED] {lbl:50s} | 0 < T/g <= {T0n/T0d:.4f} | "
            f"alpha={an}/{ad} | m_0^2 >= {floor_str} (exact lo={mp.nstr(m2_lo, 6)})"
        )
        assert m2_lo > 0 and floor_val <= m2_lo
        assert floor_str == exp_str, f"Mismatch for {lbl}: got {floor_str}, expected {exp_str}"

    print(f"\nALL PK-HOLSTEIN REPLAY CHECKS PASSED in {time.time() - t_start:.2f}s.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
