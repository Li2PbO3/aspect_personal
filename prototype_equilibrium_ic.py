#!/usr/bin/env python3
"""Prototype/validation of the proposed ASPECT-internal equilibrium initial
condition algorithm (see INTERNAL_EQUILIBRIUM_IC_DESIGN.md, section 5.4).

NOTE: reuses the vectorised re-implementation of the C++ equilibrium kernel
from region5_carbon/case6_practical_0724/verify_teq_pressure/analyze.py
(function solve_melt_fraction, K_PRM/K_PY), so run it with that path
available (the script inserts it itself).


Validates, against the authoritative external IC file
(initial_composition_for_carbon_case06_uniphase.txt), that

  1. the analytic tie-line path  c_bar = (1-s) c_s,bg + s c_l,bg ,  s = phi*
     reproduces the file's bulk composition;
  2. the generic bisection path (used for all other mixing families)
     finds the same s;
  3. the round trip  Phi(P,T,c_bar) == phi*  holds.

Uses the same kernel parameterisation as the IC generator (June2026 cMORB).
"""
import sys, os
sys.path.insert(0, "/home/zhiqianli/workspace/region5_carbon/case6_practical_0724/verify_teq_pressure")
import numpy as np
from analyze import solve_melt_fraction, K_PY, RHO_CRUST, RHO_MANTLE, H_MOHO, G

CASE = "/home/zhiqianli/workspace/region5_carbon/case6_practical_0724"
ic = np.loadtxt(os.path.join(CASE, "initial_composition_for_carbon_case06_uniphase.txt"))
tp = np.loadtxt(os.path.join(CASE, "initial_temperature_for_carbon_case06_uniphase.txt"))
ic = ic[ic[:, 0] == 0]; tp = tp[tp[:, 0] == 0]
y = ic[:, 1]
phi_file = ic[:, 2]
cbar_file = np.stack([ic[:, 3], ic[:, 4], ic[:, 5]])          # (3, n)
T_C = tp[:, 2] - 273.15
P = RHO_CRUST * G * H_MOHO + RHO_MANTLE * G * (250000.0 - y - H_MOHO)

C_BG = np.array([0.7, 0.3 - 200e-6/0.30, 200e-6/0.30])         # 200 ppm CO2
c_bg = np.repeat(C_BG[:, None], y.size, axis=1)

def state(c, P, T):
    f = solve_melt_fraction(P, T, c, K_PY)
    K = np.stack([np.exp(np.clip((k["L"]/k["r"])*(1.0/T - 1.0/np.where(P > k["Pthr"],
            (k["T0"]+k["A"]*k["Pthr"]+k["B"]*k["Pthr"]**2)+(k["A"]+2*k["B"]*k["Pthr"])*(P-k["Pthr"]),
             k["T0"]+k["A"]*P+k["B"]*P*P)), -700, 700)) for k in K_PY], axis=0)
    c_l = c/(f + (1-f)*K)
    c_s = K*c_l
    return f, c_s, c_l

# background state
f_bg, cs_bg, cl_bg = state(c_bg, P, T_C)

# --- (1) analytic tie-line path:  s = phi* ---------------------------------
s_ana = phi_file.copy()
cbar_ana = (1.0 - s_ana)*cs_bg + s_ana*cl_bg
f_ana, cs_ana, cl_ana = state(cbar_ana, P, T_C)

# --- (2) generic bisection on s -------------------------------------------
def phi_of_s(s):
    c = (1.0 - s)*cs_bg + s*cl_bg
    return solve_melt_fraction(P, T_C, c, K_PY)

s_bis = np.zeros_like(phi_file)
interior = (y > 1000) & (y < 149000)
idx = np.where(interior)[0]

# vectorised bisection: all nodes at once (Phi is monotone in s)
lo = np.zeros(y.size); hi = np.ones(y.size)
for _ in range(70):
    mid = 0.5*(lo+hi)
    c = (1.0-mid)*cs_bg + mid*cl_bg
    f = solve_melt_fraction(P, T_C, c, K_PY)
    below = f < phi_file
    lo = np.where(below, mid, lo)
    hi = np.where(below, hi, mid)
s_bis = np.where(phi_file > 0.0, 0.5*(lo+hi), 0.0)

print(f"nodes (interior): {idx.size}")
print(f"[1] c_bar(analytic) - c_bar(file) : max={np.abs(cbar_ana-cbar_file)[:,idx].max():.3e}"
      f"  (per-component max={np.abs(cbar_ana-cbar_file)[:,idx].max(axis=1)})")
print(f"[1] phi(round trip) - phi*(file)  : max={np.abs(f_ana-phi_file)[idx].max():.3e}")
print(f"[2] s(bisection)   - s(analytic)  : max={np.abs(s_bis-s_ana)[idx].max():.3e}")
print(f"[2] c_bar(bisection) - c_bar(file): max={np.abs(((1-s_bis)*cs_bg+s_bis*cl_bg)-cbar_file)[:,idx].max():.3e}")
ipk = int(np.argmin(np.abs(y - 25000.0)))
print(f"\nwave peak (depth {250-y[ipk]/1e3:.1f} km):")
print(f"  phi* = {phi_file[ipk]:.6e}  s_bis = {s_bis[ipk]:.12f}")
print(f"  c_bar file      = {np.round(cbar_file[:,ipk], 12)}")
print(f"  c_bar analytic  = {np.round(cbar_ana[:,ipk], 12)}")
print(f"  c_bar bisection = {np.round(((1-s_bis)*cs_bg+s_bis*cl_bg)[:,ipk], 12)}")
print(f"  c_s,bg = {np.round(cs_bg[:,ipk],9)}   c_l,bg = {np.round(cl_bg[:,ipk],9)}")
print(f"  phi_bg = {f_bg[ipk]:.6e}")
nad = int(np.sum(phi_file[idx] <= 0.0))
print(f"\nall-solid nodes in interior: {nad} ; max phi* = {phi_file.max():.4e}")
