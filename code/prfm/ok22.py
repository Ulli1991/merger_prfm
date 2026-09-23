"""Fits from Ostriker & Kim 2022 (ApJ 936, 137), TIGRESS calibration.
P in K cm^-3 (P/k_B), Upsilon in km/s, Sigma_SFR in Msun kpc^-2 yr^-1."""
import numpy as np
def ups_th(PDE):   return 10 ** (-0.506 * np.log10(PDE) + 4.45)      # eq. 26a
def ups_turb(PDE): return 10 ** (-0.060 * np.log10(PDE) + 2.81)      # eq. 26b
def ups_tot(PDE):  return 10 ** (-0.212 * np.log10(PDE) + 3.86)      # eq. 26c
def sfr_of_PDE(PDE):  return 10 ** (1.21 * np.log10(PDE) - 7.66)     # eq. 28c
def sfr_of_W(W):      return 10 ** (1.17 * np.log10(W) - 7.32)       # eq. 28b (W_2p)
def sfr_of_Ptot(P):   return 10 ** (1.18 * np.log10(P) - 7.43)       # eq. 28a (P_tot,2p)
def sigma_eff_PDE(PDE): return 12.0 * (PDE / 1e4) ** 0.22            # sec 5, for PDE > 1e4
def sigma_eff_Ptot(P):  return 9.8 * (P / 1e4) ** 0.15
# yields vs Sigma_SFR (eq. 25): Ups_th = 110 (S/0.01)^-0.4, Ups_turb = 330 (S/0.01)^-0.05, Ups_tot = 740 (S/0.01)^-0.18
def ups_th_sfr(S):   return 110.0 * (S / 0.01) ** -0.4
def ups_turb_sfr(S): return 330.0 * (S / 0.01) ** -0.05
def ups_tot_sfr(S):  return 740.0 * (S / 0.01) ** -0.18
# OK22 conventions: "2p" = warm+cold gas with T < 2e4 K; SFR from stars younger than 40 Myr;
# midplane values averaged over |z| < ~ H/2 and over time; P_DE = pi G Sigma^2/2 + Sigma (2 G rho_sd)^1/2 sigma_eff

# ---- TIGRESS-NCR, metallicity-dependent (Kim et al. 2024, eqs 14-18).  W4 = W / (1e4 k_B cm^-3 K), Zp = Z/Z_sun.
# Calibrated over Zp 0.1-3, Sigma_gas 5-150, Sigma_star 1-50 Msun/pc^2, with ionising radiation ("early feedback")
# that TIGRESS-classic (and this run) does not have.
def ncr_ups_th(W4, Zp=0.1):    return 390.0 * W4 ** -0.46 * Zp ** -0.53
def ncr_ups_turb(W4, Zp=0.1):  return 561.0 * W4 ** -0.21 * Zp ** -0.04
def ncr_ups_mag(W4, Zp=0.1):   return 578.0 * W4 ** -0.40 * Zp ** -0.44
def ncr_ups_nonth(W4, Zp=0.1): return 1.17e3 * W4 ** -0.22 * Zp ** -0.18
def ncr_ups_tot(W4, Zp=0.1):   return 1.65e3 * W4 ** -0.29 * Zp ** -0.27
def ncr_eps_dyn(W4, Zp=0.1):   return 0.0071 * W4 ** 0.41 * Zp ** 0.30
def jeff_eps_dyn(P4):          return 0.012 * P4 ** 0.43      # Jeffreson et al. 2026, from TIGRESS-classic
