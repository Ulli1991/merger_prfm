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
