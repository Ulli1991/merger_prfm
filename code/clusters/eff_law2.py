"""Close the efficiency gap?  With the pre-onset clump state now including B, c_s, v_A and the measured density-PDF width sig_s:
(1) magnetic support: alpha_tot = alpha_kin (1 + 1/M_A^2) with M_A = sigma_3d / v_A; plasma beta = 2 c_s^2 / v_A^2; does eps_int vs
    alpha_tot collapse the phase offsets that eps_int vs alpha_kin leaves (0.8 dex)?
(2) Federrath & Klessen 2012 multi-freefall (PN) efficiency per free-fall time:
    eps_ff = eps/(2 phi_t) exp(3 sig_s^2/8) [1 + erf((sig_s^2 - s_crit)/sqrt(2 sig_s^2))], s_crit = ln(0.067 theta^-2 alpha_vir M^2 f(beta)),
    f(beta) = (1 + 1/beta)^-1, eps/(2 phi_t) = 0.245, theta = 0.35 (their fit).  sig_s: measured (std of ln n in the clump) and the
    formula ln(1 + b^2 M^2 beta/(beta+1)) with b = 0.4.  Observed per-free-fall efficiency: eps_int * t_ff(n_mean,on) / t_SF with
    t_SF = lifetime after onset (min 1 Myr).  Compare: Spearman, per-phase offsets log(obs/pred), pooled residual.
(3) phase offsets of every model side by side."""
import sys, numpy as np
sys.path.insert(0, '/raven/u/uli/dwarf_merger'); import common as C
from scipy.stats import spearmanr
from scipy.special import erf
EP = [*C.PHASES]
G_PC = 4.30091e-3
for cat in ('cl100', 'cl10'):
    d = np.load(f'{C.DATADIR}/clouds/lineages_{cat}.npz'); L = d['rows']; col = {n: i for i, n in enumerate(d['names'])}
    g = (L[:, col['Mstar']] > 0) & (L[:, col['M_max']] >= 300) & (L[:, col['t_start']] > 5) & (L[:, col['t_end']] < 221) & np.isfinite(L[:, col['alpha_on']]) & (L[:, col['alpha_on']] > 0) & np.isfinite(L[:, col['vA_on']]) & (L[:, col['vA_on']] > 0)
    Lg = L[g]; c = lambda k: Lg[:, col[k]]
    eps = c('eps_int'); al = c('alpha_on'); sig = c('sig_on'); vA = c('vA_on'); cs = c('cs_on'); sigs = c('sigs_on'); ton = c('t_onset'); nmean = c('n_mean_on')
    MA = sig / vA; beta = 2 * cs ** 2 / vA ** 2; Mach = sig / cs; al_tot = al * (1 + 1 / MA ** 2)
    rho_msun = nmean * 1.4 * 1.6726e-24 / (1.989e33 / 3.086e18 ** 3); tff = np.sqrt(3 * np.pi / (32 * G_PC * rho_msun)) * 3.086e13 / 3.156e13   # Myr
    tSF = np.maximum(c('t_end') - ton + 10., 1.0)          # star-forming duration incl. the 10 Myr tail
    eps_ff_obs = eps * tff / tSF
    print(f'\n===== {cat}: {g.sum()} SF lineages =====')
    print('pre-onset clump state per phase (median): alpha_kin, M_A, beta, Mach, sig_s(meas), sig_s(formula b=0.4), B [muG], t_ff [Myr], t_SF [Myr]')
    sigs_f = np.sqrt(np.log(1 + 0.16 * Mach ** 2 * beta / (beta + 1)))
    for t0, t1, lab in EP:
        m = (ton >= t0) & (ton < t1)
        print(f'  {lab:20s} n={m.sum():4d}  alpha {np.median(al[m]):5.1f}  M_A {np.median(MA[m]):4.2f}  beta {np.median(beta[m]):5.2f}  Mach {np.median(Mach[m]):4.1f}  sig_s {np.median(sigs[m]):.2f} / {np.median(sigs_f[m]):.2f}  B {np.median(c("B_on")[m]):5.1f}  t_ff {np.median(tff[m]):4.1f}  t_SF {np.median(tSF[m]):4.1f}')
    # (1) alpha_kin vs alpha_tot
    le = np.log10(eps)
    print('\n(1) pooled Spearman(log eps_int, x) and the phase offsets of log eps_int from a pooled power-law fit in x:')
    for name, x in (('alpha_kin', al), ('alpha_tot = alpha(1+1/M_A^2)', al_tot), ('M_A', MA), ('beta', beta), ('Mach', Mach), ('sig_s', sigs), ('B', c('B_on'))):
        lx = np.log10(x); cf = np.polyfit(lx, le, 1); res = le - np.polyval(cf, lx)
        offs = [np.mean(res[(ton >= t0) & (ton < t1)]) for t0, t1, _ in EP]
        print(f'  {name:30s} rho {spearmanr(le, lx)[0]:+.2f}  slope {cf[0]:+.2f}  residual {np.std(res):.2f} dex  phase offsets ' + ' '.join(f'{o:+.2f}' for o in offs) + f'  (range {max(offs)-min(offs):.2f})')
    # two-variable fits
    for name, X in (('alpha_kin + beta', np.column_stack([np.log10(al), np.log10(beta)])), ('alpha_kin + M_A', np.column_stack([np.log10(al), np.log10(MA)])), ('alpha_kin + Mach', np.column_stack([np.log10(al), np.log10(Mach)])), ('alpha_tot + Mach', np.column_stack([np.log10(al_tot), np.log10(Mach)]))):
        A = np.column_stack([np.ones(len(le)), X]); cf = np.linalg.lstsq(A, le, rcond=None)[0]; res = le - A @ cf
        offs = [np.mean(res[(ton >= t0) & (ton < t1)]) for t0, t1, _ in EP]
        print(f'  {name:30s} coefs {cf[1]:+.2f} {cf[2]:+.2f}  residual {np.std(res):.2f} dex  phase offsets ' + ' '.join(f'{o:+.2f}' for o in offs) + f'  (range {max(offs)-min(offs):.2f})')
    # (2) FK12
    print('\n(2) Federrath & Klessen 2012 multi-freefall PN efficiency vs observed eps_ff = eps_int t_ff / t_SF')
    def fk12(ss, alpha, mach, b_):
        scrit = np.log(0.067 / 0.35 ** 2 * alpha * mach ** 2 / (1 + 1 / b_)); s2 = ss ** 2
        return 0.245 * np.exp(3 * s2 / 8) * (1 + erf((s2 - scrit) / np.sqrt(2 * s2)))
    lo = np.log10(np.maximum(eps_ff_obs, 1e-6))
    for name, ss, alpha_use, f_b in (('measured sig_s, alpha_kin, with beta', sigs, al, beta), ('formula sig_s, alpha_kin, with beta', sigs_f, al, beta), ('measured sig_s, alpha_kin, no field (beta->inf)', sigs, al, np.full_like(beta, 1e6)), ('measured sig_s, alpha_tot, with beta', sigs, al_tot, beta)):
        pred = fk12(ss, alpha_use, Mach, f_b); lp = np.log10(np.maximum(pred, 1e-8)); ok = pred > 1e-8
        res = lo[ok] - lp[ok]; offs = [np.mean(res[((ton >= t0) & (ton < t1))[ok]]) for t0, t1, _ in EP]
        print(f'  {name:48s} rho {spearmanr(lo[ok], lp[ok])[0]:+.2f}  median log(obs/pred) {np.median(res):+.2f}  scatter {np.std(res):.2f} dex  phase offsets ' + ' '.join(f'{o:+.2f}' for o in offs) + f'  (range {max(offs)-min(offs):.2f})   pred med per phase ' + ' '.join(f'{np.median(pred[(ton >= t0) & (ton < t1)]):.3f}' for t0, t1, _ in EP))
    print('  observed eps_ff median per phase: ' + ' '.join(f'{np.median(eps_ff_obs[(ton >= t0) & (ton < t1)]):.4f}' for t0, t1, _ in EP) + '   (eps_int: ' + ' '.join(f'{np.median(eps[(ton >= t0) & (ton < t1)]):.4f}' for t0, t1, _ in EP) + ')')
