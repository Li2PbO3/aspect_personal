#!/usr/bin/env python3
"""Probe of the mixing-family anchors in the proposed ASPECT-internal
equilibrium initial condition module.

See INTERNAL_EQUILIBRIUM_IC_DESIGN.md sections 5.4.2 and 7.

Question answered here: in a region where the background is **all solid**
(phi_bg = 0), is the "tie line" family
        c_bar(s) = (1-s) A + s B,   A = c_s,bg,  B = c_l,bg
a usable one-dimensional family?

Answer: no.  When phi_bg = 0 the equilibrium solver returns
c_s,bg = c_bar_bg (which sums to 1) but c_l,bg = c_bar_bg / K, which is a
*hypothetical* liquid that does not even sum to one.  All members with s > 0
therefore have sum(c_bar) < 1, i.e. they are not valid bulk compositions, and
ASPECT's equilibrium solver does not renormalise.  The family is only a real
tie line where phi_bg > 0, which is exactly the region the python generator
keeps via its mask = (shape > 0) & (f_bg > 0).

The script also shows which family *can* create melt in an all-solid region:
enriching the rock, e.g. background -> pure cMORB.

NOTE: reuses the vectorised re-implementation of the C++ equilibrium kernel
from region5_carbon/case6_practical_0724/verify_teq_pressure/analyze.py
(solve_melt_fraction, K_PY).
"""

import sys

import numpy as np

sys.path.insert(0, "/home/zhiqianli/workspace/region5_carbon/case6_practical_0724/verify_teq_pressure")
from analyze import solve_melt_fraction, K_PY, RHO_CRUST, RHO_MANTLE, H_MOHO, G  # noqa: E402


def P_of_depth(d):                      # the python IC profile
    return RHO_CRUST * G * H_MOHO + RHO_MANTLE * G * (d - H_MOHO)


def Tm(P, k):
    if P <= k["Pthr"]:
        return k["T0"] + k["A"] * P + k["B"] * P * P
    Tthr = k["T0"] + k["A"] * k["Pthr"] + k["B"] * k["Pthr"] ** 2
    return Tthr + (k["A"] + 2 * k["B"] * k["Pthr"]) * (P - k["Pthr"])


C_BG = np.array([0.7, 0.3 - 200e-6 / 0.30, 200e-6 / 0.30])      # 200 ppm CO2


def eq(c, P, T_C):
    c = np.asarray(c, float)
    f = float(solve_melt_fraction(np.array([P]), np.array([T_C]), c[:, None], K_PY)[0])
    K = np.array([np.exp(np.clip((k["L"] / k["r"]) * (1.0 / T_C - 1.0 / Tm(P, k)), -700, 700))
                  for k in K_PY])
    cl = c / (f + (1 - f) * K)
    return f, K * cl, cl, K


def T_of_depth(d):                      # case6 PracticalTemperatureProfile
    k, q0, qm, hr, T0 = 3.35, 65e-3, 28e-3, 10e3, 273.15
    alpha, Cp, g = 2e-5, 1000.0, 9.81
    Tp = 1553.70 / np.exp(alpha * g * 139994.4 / Cp)
    return min(T0 + qm * d / k + ((q0 - qm) * hr / k) * (1 - np.exp(-d / hr)),
               Tp * np.exp(alpha * g * d / Cp))


def main():
    for depth in (110e3, 130e3, 225e3):
        P = P_of_depth(depth)
        T_C = T_of_depth(depth) - 273.15
        f_bg, cs, cl, K = eq(C_BG, P, T_C)
        print("=" * 78)
        print(f"depth {depth/1e3:.0f} km :  P = {P:.4e} Pa,  T = {T_C:.2f} C")
        print(f"  K = {np.round(K, 6)}      (K < 1 <=> highly incompatible)")
        print(f"  background: phi_bg = {f_bg:.6e}")
        print(f"    c_s,bg = {np.round(cs, 9)}   sum = {cs.sum():.6f}")
        print(f"    c_l,bg = {np.round(cl, 9)}   sum = {cl.sum():.6f}   <-- != 1 when f = 0")

        ss = (0.0, 0.25, 0.5, 0.75, 1.0)
        sums = [((1 - s) * cs + s * cl).sum() for s in ss]
        print("  sum(c_bar(s)) along the tie-line family:")
        print("    " + "  ".join(f"s={s:.2f}: {v:.4f}" for s, v in zip(ss, sums)))
        if min(sums) < 1.0 - 1e-12:
            print("    *** these members are NOT valid bulk compositions (sum != 1);")
            print("    *** ASPECT does not renormalise, so the tie-line family is")
            print("    *** undefined here and the Phi(s) row below is meaningless. ***")
        for name, A, B in (("tie line  (c_s,bg -> c_l,bg)", cs, cl),
                           ("bg -> pure cmorb", C_BG, np.array([0.0, 0.0, 1.0])),
                           ("bg -> pure morb", C_BG, np.array([0.0, 1.0, 0.0]))):
            phi = [eq((1 - s) * A + s * B, P, T_C)[0] for s in ss]
            ok = "phi*=1e-3 reachable" if max(phi) >= 1e-3 else "**unreachable**"
            print(f"    {name:32s} Phi(s) = {['%.3e' % v for v in phi]}  {ok}")


if __name__ == "__main__":
    main()
