"""Excursion-set mass function of self-gravitating structures in a turbulent,
rotating, vertically stratified gas layer (Hopkins 2012, MNRAS 423, 2037),
driven by the PRFM state of a patch:

    Sigma  [Msun/pc^2]   gas surface density
    h      [pc]          gas scale height (measured or from equilibrium)
    vt     [km/s]        turbulent velocity dispersion on scale h
    cs     [km/s]        sound speed of the star-forming (cold) phase
    kappa, Omega [1/Myr] epicyclic and orbital frequency

Log-density field on scale R is Gaussian with variance
    S(R) = int_{k<1/R} ln[1 + b^2 v_t^2(k) / (c_s^2 + kappa^2/k^2)] dln k
    v_t^2(k) = vt^2 min[1, (k h)^-p]                (p = 1: Burgers-like)
Collapse barrier (finite-thickness, rotating dispersion relation):
    rho_crit/rho_0 = (Q / 2 kt) (1 + h/R) [ sig^2(R)/sig^2(h) h/R + kt^2 R/h ]
    Q = sig(h) kappa / (pi G Sigma),  kt = kappa/Omega,  sig^2 = c_s^2 + v_t^2
    B(S) = ln(rho_crit/rho_0) + S/2
First-crossing distribution f(S) from the Zhang & Hui (2006) integral equation.
Mass of a structure of scale R:
    M(R) = 4 pi rho_crit h^3 [ R^2/(2h^2) + (1 + R/h) e^{-R/h} - 1 ]
Mass function per unit area:  dn/dM = (Sigma / M) f(S) |dS/dM|.
"""
import numpy as np

G = 4.498e-3          # pc^3 Msun^-1 Myr^-2
KMS = 1.0227          # km/s -> pc/Myr

def variance(R, h, vt, cs, kappa, b2=0.75, p=1.0, nk=400):
    """S(R) for an array of R [pc]; velocities in pc/Myr, kappa 1/Myr."""
    kmin = 1e-3 / max(h, 1e-3); kmax = 1.0 / R.min()
    lk = np.linspace(np.log(kmin), np.log(kmax), nk); k = np.exp(lk)
    vt2 = vt ** 2 * np.minimum(1.0, (k * h) ** (-p))
    integrand = np.log(1 + b2 * vt2 / (cs ** 2 + kappa ** 2 / k ** 2))
    cum = np.concatenate([[0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(lk))])
    return np.interp(np.log(1.0 / R), lk, cum)

def barrier(R, h, vt, cs, kappa, Omega, Sigma, p=1.0):
    """ln(rho_crit/rho_0) for R [pc]; returns also rho_crit/rho_0."""
    sig2R = cs ** 2 + vt ** 2 * np.minimum(1.0, (R / h) ** p)
    sig2h = cs ** 2 + vt ** 2
    Q = np.sqrt(sig2h) * kappa / (np.pi * G * Sigma)
    kt = kappa / max(Omega, 1e-9)
    x = R / h
    rc = (Q / (2 * kt)) * (1 + 1 / x) * (sig2R / sig2h / x + kt ** 2 * x)
    return np.log(rc), rc, Q

def first_crossing(S, B, n=1500):
    """First-crossing distribution f(S) (per unit S) for a moving barrier B(S),
    Zhang & Hui (2006) discretisation on a uniform S grid. Returns f on the input S grid."""
    Su = np.linspace(S.min(), S.max(), n); dS = Su[1] - Su[0]
    Bu = np.interp(Su, S, B); dB = np.gradient(Bu, Su)
    def P0(d, s):
        return np.exp(-d * d / (2 * s)) / np.sqrt(2 * np.pi * s)
    g1 = (Bu / Su - 2 * dB) * P0(Bu, Su)
    f = np.zeros(n); f[0] = max(g1[0], 0.0)
    for i in range(1, n):
        j = np.arange(i)
        ds = Su[i] - Su[j]; db = Bu[i] - Bu[j]
        g2 = (2 * dB[i] - db / ds) * P0(db, ds)
        g2_ii = dB[i] / np.sqrt(2 * np.pi * dS)            # diagonal limit S'->S
        w = np.full(i, dS); w[0] = 0.5 * dS
        f[i] = (g1[i] + np.sum(w * f[j] * g2)) / (1.0 - 0.5 * dS * g2_ii)
    f = np.maximum(f, 0.0)
    return np.interp(S, Su, f)

def mass_function(Sigma, h, vt_kms, cs_kms, kappa, Omega, b2=0.75, p=1.0, nR=350, Rmax_h=300.0, rho0=None, rho_fb=None):
    """Returns dict(M [Msun], dndM [1/(Msun pc^2)], R, S, B, f, Q). Inputs: Sigma Msun/pc^2, h pc,
    velocities km/s, kappa/Omega 1/Myr."""
    vt = vt_kms * KMS; cs = cs_kms * KMS
    mach = vt / cs
    Rs = h * mach ** (-2.0 / p) * 0.3                # start below the sonic scale
    R = np.logspace(np.log10(max(Rs, 1e-3 * h)), np.log10(Rmax_h * h), nR)[::-1]   # large -> small (S increasing)
    S = variance(R, h, vt, cs, kappa, b2, p)
    lnrc, rc, Q = barrier(R, h, vt, cs, kappa, Omega, Sigma, p)
    rho0 = Sigma / (2 * h) if rho0 is None else rho0
    if rho_fb is not None:      # feedback-limited: must also reach rho_fb (t_ff = t_fb) -> combined barrier
        lnrc = np.maximum(lnrc, np.log(rho_fb / rho0)); rc = np.exp(lnrc)
    B = lnrc + S / 2
    ok = np.concatenate([[True], np.diff(S) > 1e-9]); R, S, B, rc = R[ok], S[ok], B[ok], rc[ok]
    S = np.maximum(S, 1e-6)
    f = first_crossing(S, B)
    x = R / h
    M = 4 * np.pi * rc * rho0 * h ** 3 * (x ** 2 / 2 + (1 + x) * np.exp(-x) - 1)
    dSdM = np.abs(np.gradient(S, M))
    dndM = Sigma / M * f * dSdM
    return dict(M=M, dndM=dndM, R=R, S=S, B=B, f=f, Q=Q, rc=rc, mach=mach)

def mle_slope_from_mf(M, dndM, mmin, mmax=None):
    """alpha of a power law fitted (max-likelihood, continuous) to a tabulated dN/dM above mmin."""
    s = (M >= mmin) & (dndM > 0) & ((M <= mmax) if mmax else True)
    if s.sum() < 4: return np.nan
    w = dndM[s] * np.abs(np.gradient(M[s]))
    return 1 + w.sum() / np.sum(w * np.log(M[s] / mmin))

if __name__ == '__main__':
    # --- validation of the first-crossing solver against closed forms
    S = np.linspace(0.02, 6, 300)
    for B0, beta in [(1.5, 0.0), (1.0, 0.4)]:
        B = B0 + beta * S
        f = first_crossing(S, B)
        fa = B0 / (np.sqrt(2 * np.pi) * S ** 1.5) * np.exp(-(B0 + beta * S) ** 2 / (2 * S))
        print(f'barrier B={B0}+{beta}S: max |f_num/f_an - 1| in [0.3,5] = {np.max(np.abs(f[(S>0.3)&(S<5)]/fa[(S>0.3)&(S<5)]-1)):.3f}, '
              f'integral num {np.trapezoid(f,S):.3f} an {np.trapezoid(fa,S):.3f}')
    # --- fiducial patch states (kappa, Omega in km/s/kpc -> 1/Myr)
    K = 1.0227e-3
    for lab, kw in [('dwarf  ', dict(Sigma=5.0, h=150.0, vt_kms=6.0, cs_kms=0.2, kappa=60 * K, Omega=40 * K)),
                    ('merger ', dict(Sigma=30.0, h=100.0, vt_kms=8.0, cs_kms=0.2, kappa=100 * K, Omega=70 * K)),
                    ('MW-like', dict(Sigma=12.0, h=100.0, vt_kms=7.0, cs_kms=0.2, kappa=37 * K, Omega=26 * K))]:
        r = mass_function(**kw)
        M, d = r['M'], r['dndM']; ok = d > 0
        if not ok.any():
            print(lab, 'no crossings: S range', r['S'].min(), r['S'].max(), 'B range', r['B'].min(), r['B'].max()); continue
        Mc = M[ok][np.argmax(M[ok] ** 2 * d[ok])]
        al = mle_slope_from_mf(M, d, 1e3, Mc)
        print(f'{lab}: Q={r["Q"]:.2f} Mach={r["mach"]:.0f}  S max={r["S"].max():.2f}  B min={r["B"].min():.2f}  '
              f'M range {M[ok].min():.2e}-{M[ok].max():.2e}  cutoff {Mc:.2e} Msun  alpha(1e3..Mc)={al:.2f}  int f dS={np.trapezoid(r["f"], r["S"]):.2f}')


# ---------------------------------------------------------------- alternative object models
def mass_function_threshold(Sigma, h, vt_kms, cs_kms, kappa, Omega, Bth, b2=0.75, p=1.0, nR=350, Rmax_h=300.0, rho0=None):
    """Same log-density field, but objects = connected regions above a fixed density
    threshold rho_th = rho0 e^{Bth} (first crossing of a CONSTANT barrier; mass at fixed
    density M = 4 pi rho_th h^3 [x^2/2 + (1+x)e^-x - 1]).  Closed-form f(S)."""
    vt = vt_kms * KMS; cs = cs_kms * KMS
    mach = vt / cs
    Rs = h * mach ** (-2.0 / p) * 0.3
    R = np.logspace(np.log10(max(Rs, 1e-3 * h)), np.log10(Rmax_h * h), nR)[::-1]
    S = variance(R, h, vt, cs, kappa, b2, p)
    ok = np.concatenate([[True], np.diff(S) > 1e-9]); R, S = R[ok], S[ok]
    S = np.maximum(S, 1e-6)
    B = Bth + S / 2                      # constant density threshold in ln(rho/rho0); mean of ln rho is -S/2
    f = Bth / np.sqrt(2 * np.pi * S ** 3) * np.exp(-B ** 2 / (2 * S))      # exact: linear barrier Bth + S/2
    rho0 = Sigma / (2 * h) if rho0 is None else rho0
    x = R / h
    M = 4 * np.pi * rho0 * np.exp(Bth) * h ** 3 * (x ** 2 / 2 + (1 + x) * np.exp(-x) - 1)
    dSdM = np.abs(np.gradient(S, M))
    return dict(M=M, dndM=Sigma / M * f * dSdM, R=R, S=S, B=B, f=f, mach=mach)

def walks_mass_function(S, B, M, nwalk=200000, margin=np.log(2.0), seed=1):
    """Monte-Carlo random walks (sharp-k filter => Markovian) on the uniform-S grid.
    Returns (first-crossing counts per S bin, all-crossing counts) where an 'all' crossing
    is any upward crossing of the barrier after the walk has been at least `margin`
    below it since the previous crossing (hierarchical sub-objects).  Validates f(S)."""
    rng = np.random.default_rng(seed)
    n = len(S); dS = np.diff(S, prepend=0.0)
    first = np.zeros(n); allc = np.zeros(n)
    chunk = 20000
    for c0 in range(0, nwalk, chunk):
        m = min(chunk, nwalk - c0)
        d = np.zeros(m); above_before = np.zeros(m, bool); armed = np.ones(m, bool); done = np.zeros(m, bool)
        for i in range(n):
            d = d + rng.normal(0.0, np.sqrt(max(dS[i], 1e-12)), m)
            above = d >= B[i]
            cross = above & armed
            first[i] += np.sum(cross & ~done); done |= cross
            allc[i] += np.sum(cross)
            armed = np.where(cross, False, armed)
            armed |= (d < B[i] - margin)
    return first / nwalk, allc / nwalk


def rho_feedback(t_fb_myr):
    """density at which the free-fall time equals t_fb [Msun/pc^3]."""
    return 3 * np.pi / (32 * G * t_fb_myr ** 2)
