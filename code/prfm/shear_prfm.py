"""PRFM + mean-flow (shear) driving.  Turbulent energy balance per unit area of the layer:
    Sigma sigma^3 / (2H)  =  Edot_fb  +  c_s Sigma sigma H S^2        (S = mean-flow shear rate, c_s = O(1) efficiency)
feedback only:  rho sigma_fb^2 = Ups_turb Sigma_SFR   (OK22),  rho = Sigma / 2H
=> sigma^3 - 2 c_s H^2 S^2 sigma - sigma_fb^3 = 0 ;  limits  P_turb -> Ups_turb Sigma_SFR  (feedback)  or  c_s Sigma H S^2  (flow)
regulation:  W = Ups_th Sigma_SFR + rho sigma^2 + P_mag   =>  Sigma_SFR = (W - P_flow) / Ups   with  P_flow = c_s Sigma H S^2 + P_mag
Tests on the patch catalogue (pressure_sources.npz): (a) quiescent patches: P_turb vs Sigma H S^2; (b) c_s fit on all patches;
(c) does W_fb = W - P_flow predict the SFR and the quiescent fraction better than W?
usage: python shear_prfm.py [shear_key]   (default 'shear' = particle-scale gradient; later 'shear_patch' from turb_decomp)"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C, ok22
from scipy.stats import spearmanr
PU = 4.903245584337325e-06
R = dict(np.load(f'{C.DATADIR}/prfm/pressure_sources.npz'))
key = sys.argv[1] if len(sys.argv) > 1 else 'shear'
ok = (R['t'] > 10) & np.isfinite(R[key]) & (R[key] > 0) & (R['H'] > 0)
for k in R: R[k] = R[k][ok]
Sig, H, S, Pt, Pm, SFR, W, P = R['Sigma_gas'], R['H'], R[key], R['Pturb_2p'], R['Pmag_2p'], R['SigSFR_10'], R['W_2p'], R['Ptot_2p']
rho = R['rho_mid_2p'] * 1e9                                  # Msun/kpc^3
X = Sig * 1e6 * H * S**2 * PU                                # Sigma H S^2 in K cm^-3
sig_z = np.sqrt(Pt / PU / np.maximum(rho, 1e-30))            # km/s
EP = [*C.PHASES]
print(f'shear key = {key}\n(a) QUIESCENT patches (SigSFR_10 = 0): log P_turb = a + b log(Sigma H S^2); c_s = P_turb/(Sigma H S^2); sigma_z/(H S)')
print('phase                n_q   b      a     scatter  median c_s (16-84)        median sig_z/(HS)   rho_s(Pturb, S | Sigma H fixed)')
def partial(y, x, z):
    r = lambda v: np.argsort(np.argsort(v)).astype(float); ry, rx, rz = r(y), r(x), r(z); A = np.column_stack([np.ones(len(rz)), rz])
    ey = ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0]; ex = rx - A @ np.linalg.lstsq(A, rx, rcond=None)[0]; return np.corrcoef(ey, ex)[0, 1]
for t0, t1, lab in EP:
    q = (R['t'] >= t0) & (R['t'] < t1) & (SFR == 0) & (Pt > 0)
    if q.sum() < 20: print(lab, 'too few'); continue
    b, a = np.polyfit(np.log10(X[q]), np.log10(Pt[q]), 1); res = np.log10(Pt[q]) - (a + b * np.log10(X[q]))
    cs = Pt[q] / X[q]; r = sig_z[q] / (H[q] * S[q])
    print(f'{lab:20s} {q.sum():5d}  {b:5.2f}  {a:5.2f}   {res.std():.2f}    {np.median(cs):6.2f} ({np.percentile(cs,16):.2f}-{np.percentile(cs,84):.2f})      {np.median(r):5.2f}              {partial(Pt[q], S[q], Sig[q]*H[q]):+.2f}')
print('\n(b) ALL patches: P_turb = Ups_turb,OK22(P) Sigma_SFR + c_s Sigma H S^2 ; c_s by least squares in log space; variance explained')
for t0, t1, lab in EP:
    m = (R['t'] >= t0) & (R['t'] < t1) & (Pt > 0)
    Pfb = ok22.ups_turb(P[m]) * 1e5 * SFR[m] * (1.989e33 / 3.156e7 / (3.0857e21)**2) / 1.3807e-16
    grid = np.logspace(-2, 1.5, 71); best = min(grid, key=lambda c: np.std(np.log10(Pt[m]) - np.log10(Pfb + c * X[m])))
    s_fb = np.std(np.log10(Pt[m][Pfb > 0]) - np.log10(Pfb[Pfb > 0])); s_both = np.std(np.log10(Pt[m]) - np.log10(Pfb + best * X[m])); s_flow = np.std(np.log10(Pt[m]) - np.log10(best * X[m]))
    off_fb = np.median(np.log10(Pt[m][Pfb > 0]) - np.log10(Pfb[Pfb > 0])); off_both = np.median(np.log10(Pt[m]) - np.log10(Pfb + best * X[m]))
    print(f'{lab:20s} n={m.sum():5d}  c_s={best:5.2f}   scatter[dex]: feedback-only {s_fb:.2f} (offset {off_fb:+.2f}, SF patches only)  flow-only {s_flow:.2f}  both {s_both:.2f} (offset {off_both:+.2f})')
print('\n(c) SFR prediction: W_fb = max(W - c_s Sigma H S^2 - P_mag, 0).  Binned <SigSFR> vs OK22 sfr_of_W(bin) using W or W_fb as the abscissa; quiescent-patch prediction')
cs_glob = 1.0
for t0, t1, lab in EP:
    m = (R['t'] >= t0) & (R['t'] < t1)
    Wfb = np.maximum(W[m] - cs_glob * X[m] - Pm[m], 0)
    def binned_rms(x, s):
        edges = np.logspace(2.5, 5.75, 14); out = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            k = (x >= lo) & (x < hi)
            if k.sum() >= 30 and s[k].mean() > 0: out.append(np.log10(s[k].mean()) - np.log10(ok22.sfr_of_W(np.sqrt(lo * hi))))
        out = np.array(out); return (np.sqrt(np.mean(out**2)), np.mean(out), len(out)) if len(out) else (np.nan, np.nan, 0)
    rW = binned_rms(W[m], SFR[m]); rF = binned_rms(Wfb[m > -1] if False else Wfb, SFR[m])
    pred_q = Wfb <= 0; obs_q = SFR[m] == 0
    tp = np.mean(obs_q[pred_q]) if pred_q.any() else np.nan
    print(f'{lab:20s} rms log offset from OK22 (bins): using W {rW[0]:.2f} (mean {rW[1]:+.2f}, {rW[2]} bins) | using W_fb {rF[0]:.2f} (mean {rF[1]:+.2f}, {rF[2]} bins)'
          f' | predicted quiescent {pred_q.mean():.2f} of patches, observed quiescent {obs_q.mean():.2f}; of predicted-quiescent, {tp:.2f} are quiescent')
