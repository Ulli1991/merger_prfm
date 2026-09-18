import base64
FIG='/ptmp/uli/dwarf_merger/prfm/figs'
def img(name,alt): return f'<div class="frame"><img alt="{alt}" src="data:image/png;base64,{base64.b64encode(open(f"{FIG}/{name}.png","rb").read()).decode()}"></div>'
steps=[
 dict(n='0',name='story_0_sep',q='What happens when, in this merger?',h='The merger as a clock',alt='Nuclear separation against time',
  eq='four phases:&nbsp; before first passage (&lt; 50 Myr) · between passages (50–106) · coalescence (106–150) · remnant (&gt; 150)',
  see='Two gas-rich dwarfs on a bound orbit. Two close passages, first coalescence at 106 Myr, a re-separation to 1.6 kpc, and a relaxed remnant after 150 Myr.',
  means='Every plot below uses this time axis and these four phases. Read the phases as four different physical situations: an isolated disc, a disturbed disc, a collision, and a relaxed remnant.'),
 dict(n='1',name='story_1_equilibrium',q='Is the gas layer held up by its pressure?',h='Yes, in every phase but coalescence, once the weight is that of the layer',alt='Pressure over weight against time, full column versus layer',
  eq='balance:&nbsp; P<sub>tot</sub> = 𝒲 &nbsp;&nbsp;with&nbsp;&nbsp; P<sub>tot</sub> = P<sub>th</sub> + ρσ<sub>z</sub><sup>2</sup> + (B<sup>2</sup>−2B<sub>z</sub><sup>2</sup>)/8π &nbsp;&nbsp;and&nbsp;&nbsp; 𝒲 = ∫<sub>|z|&lt;0.5 kpc</sub> ρ g<sub>z</sub> dz',
  see='The light curve is the weight of the full 3 kpc column: after coalescence the pressure falls to a third of it. The dark curve is the weight of the gas layer alone, within 0.5 kpc of the midplane: 1.00, 0.96, 0.66, 0.96 in the four phases.',
  means='The full-column failure is geometry, not missing physics. Five to fifteen per cent of the gas mass is tidal debris high above the plane; because dark-matter gravity grows with height, that mass carries 40 to 50 per cent of the column weight. The debris is not a layer in any orientation, streams through the plane at 25 to 30 km/s, and forms no stars. For the layer itself, PRFM’s first claim holds to five per cent except during the 40 Myr of coalescence, when even the layer is out of balance by a third.'),
 dict(n='2',name='story_2_feedback',q='Is that pressure supplied by star formation?',h='Only in the disc phase and inside the bursts',alt='Share of the midplane pressure supplied by feedback against time',
  eq='regulation:&nbsp; P<sub>tot</sub> = Υ<sub>tot</sub>(P) Σ<sub>SFR</sub>,&nbsp; Υ<sub>tot</sub> = 10<sup>3.86</sup> P<sup>−0.21</sup> km s<sup>−1</sup> &nbsp;&nbsp;⇒&nbsp;&nbsp; Σ<sub>SFR</sub> = 10<sup>−7.32</sup> 𝒲<sup>1.17</sup> M<sub>☉</sub> yr<sup>−1</sup> kpc<sup>−2</sup>',
  see='The dark area is the fraction of the measured pressure that the calibrated yield times the actual star formation rate accounts for. It is 0.5 to 0.95 before coalescence. Afterwards it drops to 1 to 15 per cent for two long stretches, 105 to 125 and 145 to 195 Myr, and comes back to one only inside the two bursts.',
  means='PRFM’s second claim, that feedback supplies the pressure, holds in the disc and fails between the bursts. Between bursts the layer is supported by something else and needs no stars: the star formation rate sits one to two dex below the PRFM value for the measured pressure. Star formation becomes intermittent, and the next two steps identify the something else.'),
 dict(n='3',name='story_3_shear',q='What supports the layer when feedback does not?',h='Turbulence driven by the shear of the merger flow',alt='Turbulent pressure against surface density times thickness times shear rate squared for quiescent patches',
  eq='energy balance:&nbsp; Σσ<sup>3</sup>/2H = Ė<sub>fb</sub> + ε<sub>S</sub> Σ σ H S<sup>2</sup> &nbsp;&nbsp;⇒&nbsp;&nbsp; P<sub>turb</sub> → ε<sub>S</sub> Σ H S<sup>2</sup> when Ė<sub>fb</sub> = 0,&nbsp; ε<sub>S</sub> = 0.02',
  see='Each point is a patch with no star formation in the last 10 Myr. Its turbulent pressure, with the bulk part of the motion removed, scales with Σ H S², where S is the shear rate of the mean flow measured across the patch. The line is the closure with one constant.',
  means='Dissipation on a crossing time balances feedback injection plus the work the mean-flow shear does on the turbulence. With ε<sub>S</sub> = 0.02 the closure reproduces the turbulent pressure to 0.5 dex with no offset in all four phases, where feedback alone is off by up to two dex. The constant is fitted, not derived; an order-unity argument overpredicts by fifty. About half of the measured vertical velocity variance is bulk flow rather than turbulence proper.'),
 dict(n='4',name='story_4_dynamo',q='What else grows with the merger?',h='A dynamo field that reaches equipartition in the remnant',alt='Magnetic field strength and magnetic pressure share against time',
  eq='P<sub>mag</sub> = B<sup>2</sup>/8π;&nbsp; seed 10 nG,&nbsp; B ∝ e<sup>t/10 Myr</sup> → 1.5 μG by 60 Myr,&nbsp; 5 μG at coalescence,&nbsp; 12–20 μG mass-weighted in the remnant;&nbsp; P<sub>mag</sub>/(P<sub>th</sub>+P<sub>turb</sub>): 0.05 → 0.5–1.4',
  see='The solid curve is the mass-weighted field over the patches, dominated by the dense star-forming gas; the dashed curve is the magnetic share of the pressure in the median patch. The table gives both per density class.',
  means='The seed is not the issue: the dynamo saturates within 60 Myr, before the merger matters. What grows afterwards is the saturation level, with compression and the stronger turbulence. In the disc phase the field is a few per cent of the pressure and the PRFM test is unaffected. In the remnant the magnetic pressure equals or exceeds thermal plus turbulent, and the Maxwell stress carries 15 to 30 per cent of the vertical support. It enters the story as a measured term, not a modelled one.',
  extra='''<div class="tablewrap"><table>
<thead><tr><th>phase</th><th>diffuse 1–10 M<sub>☉</sub> pc<sup>−2</sup></th><th>10–30</th><th>&gt; 30</th><th>star-forming patches</th><th>mass-weighted mean</th></tr></thead>
<tbody>
<tr><td>before first passage</td><td>0.04 μG · 0.00</td><td>1.0 μG · 0.01</td><td>—</td><td>0.08 μG · 0.00</td><td>0.3 μG</td></tr>
<tr><td>between passages</td><td>0.6 μG · 0.03</td><td>2.3 μG · 0.06</td><td>—</td><td>1.4 μG · 0.04</td><td>1.4 μG</td></tr>
<tr><td>coalescence</td><td>1.5 μG · 0.24</td><td>5.1 μG · 0.40</td><td>12 μG · 0.40</td><td>3.9 μG · 0.33</td><td>5.5 μG</td></tr>
<tr><td>remnant</td><td>1.8 μG · 0.58</td><td>9.5 μG · 0.95</td><td>24 μG · 1.16</td><td>6.4 μG · 0.79</td><td>13 μG</td></tr>
</tbody></table></div>
<p class="note">Median field strength per density class and, after the dot, magnetic pressure over thermal plus turbulent pressure in the same patches.</p>'''),
 dict(n='5',name='story_5_mf',q='What do the star clusters record?',h='A mass function with a fixed shape and a rising top',alt='Bound-cluster mass function per phase with Schechter fits',
  eq='dN/dM ∝ M<sup>−α</sup> exp(−M/M<sub>c</sub>):&nbsp; α = 1.32, 1.20, 1.40, 1.34;&nbsp; M<sub>c</sub> = 3.8×10<sup>3</sup>, 5.9×10<sup>3</sup>, 1.7×10<sup>4</sup>, 1.0×10<sup>5</sup> M<sub>☉</sub>',
  see='Bound clusters above 300 M<sub>☉</sub> per phase, with Schechter fits. The slope is the same in every phase while the truncation rises by a factor 25.',
  means='The merger moves the top of the mass function, not its shape. The function is curved everywhere: the local slope is 1.0 to 1.4 below a thousand solar masses and reaches 2 only within a factor of a few of the truncation, which is where the canonical −2 of single power-law fits comes from. The slope has no theory in this picture. The rest of the story is about the top.'),
 dict(n='6',name='story_6_reservoir',q='What sets the mass of the largest cluster in a patch?',h='Half of what the patch formed',alt='Most massive bound cluster against young stellar mass per patch',
  eq='M<sub>max</sub> = f<sub>max</sub> M<sub>young</sub>,&nbsp;&nbsp; f<sub>max</sub> = 0.53 (16–84 per cent: 0.27–0.88),&nbsp; rank correlation 0.88 over 1786 patches',
  see='Each point is a patch: the most massive bound cluster against all stars the patch formed in the last 10 Myr. The relation is tight and the ratio does not depend on weight, Toomre Q or any fragmentation scale.',
  means='Star formation proceeds as discrete events, and each event puts about half its stars into one bound cluster. So the top of the mass function is the size of the largest event, a reservoir limit, not a fragmentation scale of the turbulent gas. Above events of 10<sup>5</sup> M<sub>☉</sub>, in the coalescence bursts, the capture fraction collapses to a few per cent and the event fragments into many clusters.'),
 dict(n='7',name='story_7_ceiling',q='Can PRFM predict the largest cluster?',h='Yes: the ceiling follows from steps 2 and 6 with nothing fitted',alt='Largest bound cluster against the layer weight of the heaviest star-forming patches with the PRFM prediction',
  eq='M<sub>max</sub> = f<sub>max</sub> · Σ<sub>SFR</sub>(𝒲) · A · τ &nbsp;&nbsp;⇒&nbsp;&nbsp; log M<sub>max</sub> = 1.17 log 𝒲 − 1.22 &nbsp;&nbsp;(f<sub>max</sub> = 0.5, A = 0.25 kpc<sup>2</sup>, τ = 10 Myr)',
  see='For each 25 Myr interval, the largest bound cluster formed against the layer weight of the heaviest star-forming patches, the 90th percentile. The line is the PRFM star formation rate at that weight, times the patch area and the 10 Myr formation window, times one half. Nothing is fitted.',
  means='Over the seven intervals with star formation the prediction matches to a median +0.06 dex with 0.22 dex scatter, and the measured slope is 1.07 against 1.17. Using the median weight instead underpredicts by 3.5, which is step 6 again: the largest cluster forms in the heaviest patch. The one open circle is the quiescent interval, 160 to 185 Myr: the highest weight of all and a largest cluster fifty times below the line, because star formation was off (step 2). The weight sets the ceiling of the cluster mass function; the duty cycle of the regulation decides whether it is reached.'),
]
body=''
for s in steps:
    body+=f'''<section>
  <div class="eyebrow">step {s['n']} · {s['q']}</div>
  <h2>{s['h']}</h2>
  {img(s['name'],s['alt'])}
  <div class="eq">{s['eq']}</div>
  <p><b>What the plot shows.</b> {s['see']}</p>
  <p><b>What it means.</b> {s['means']}</p>
  {s.get('extra','')}
</section>
'''
html='''<title>Dwarf Merger Cluster MF</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Serif:wght@500;600&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400&display=swap">
<style>
:root{--bg:#f3f5f7;--panel:#ffffff;--ink:#1b2430;--muted:#5b6b7a;--line:#d5dbe1;--accent:#1f6f78;--tint:#e8eef2;}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#141a21;--panel:#1c242d;--ink:#e6ebef;--muted:#9aa8b5;--line:#2c3742;--accent:#5fc0c8;--tint:#232d38;}}
:root[data-theme="dark"]{--bg:#141a21;--panel:#1c242d;--ink:#e6ebef;--muted:#9aa8b5;--line:#2c3742;--accent:#5fc0c8;--tint:#232d38;}
body{background:var(--bg);color:var(--ink);font-family:"IBM Plex Sans",system-ui,sans-serif;font-size:15.5px;line-height:1.55;padding-inline:24px;padding-block:36px 56px;}
main{max-width:880px;margin:0 auto;display:grid;gap:52px;}
h1{font-family:"IBM Plex Serif",Georgia,serif;font-weight:600;font-size:1.75rem;margin:0 0 10px;text-wrap:balance;}
h2{font-family:"IBM Plex Serif",Georgia,serif;font-weight:600;font-size:1.2rem;margin:0;text-wrap:balance;}
p{max-width:72ch;margin:0;}
p b{font-weight:600;}
.summary{display:grid;gap:10px;max-width:72ch;}
.summary p{color:var(--ink);}
.eyebrow{font-family:"IBM Plex Mono",monospace;font-size:.76rem;letter-spacing:.06em;text-transform:uppercase;color:var(--accent);}
section{display:grid;gap:14px;}
.frame{background:#ffffff;border:1px solid var(--line);padding:10px;overflow-x:auto;}
.frame img{display:block;max-width:100%;height:auto;margin:0 auto;}
.eq{font-family:"IBM Plex Mono",monospace;background:var(--tint);padding:10px 14px;font-size:.92rem;line-height:1.6;overflow-x:auto;}
.defs{background:var(--panel);border:1px solid var(--line);padding:14px 18px;display:grid;grid-template-columns:auto 1fr;gap:6px 16px;font-size:.92rem;max-width:72ch;}
.defs dt{font-family:"IBM Plex Mono",monospace;color:var(--accent);white-space:nowrap;}
.defs dd{margin:0;}
.chain{background:var(--tint);border-left:3px solid var(--accent);padding:10px 16px;display:grid;gap:6px;max-width:72ch;}
.chain .row{display:grid;grid-template-columns:2.2em 1fr;gap:8px;align-items:start;}
.chain .row>span{font-family:"IBM Plex Mono",monospace;color:var(--accent);}
.tablewrap{overflow-x:auto;}table{border-collapse:collapse;font-size:.88rem;font-variant-numeric:tabular-nums;}th,td{padding:5px 12px;text-align:right;border-bottom:1px solid var(--line);white-space:nowrap;}th:first-child,td:first-child{text-align:left;}thead th{font-weight:500;color:var(--muted);font-family:"IBM Plex Mono",monospace;font-size:.76rem;letter-spacing:.04em;text-transform:uppercase;}
.note{color:var(--muted);font-size:.9rem;}
.meta{display:flex;flex-wrap:wrap;gap:10px 28px;font-family:"IBM Plex Mono",monospace;font-size:.78rem;color:var(--muted);border-top:1px solid var(--line);padding-top:14px;}
.meta span b{color:var(--ink);font-weight:500;}
</style>
<main>
<header>
  <div class="eyebrow">dwarf merger · pressure-regulated star formation · cluster mass function</div>
  <h1>From the weight of the gas layer to the most massive star cluster</h1>
  <div class="summary">
    <p><b>The argument in five sentences.</b> Pressure-regulated star formation (PRFM) says a gas layer is held up by its own pressure, and that the pressure is made by the stars it forms, so the star formation rate follows from the weight of the layer. In an isolated dwarf disc both parts hold. In the merger the first part still holds for the gas layer once tidal debris is left out of the weight, but the second part fails between bursts: the pressure is then made by the shear of the merger flow and by a magnetic field, not by stars, and star formation switches off until that support gives way. When it switches on, the heaviest patches form stars at the PRFM rate, and half of what a patch forms in 10 Myr ends up in its largest cluster. Multiplying those two measured facts predicts the most massive cluster in the galaxy from the weight of its heaviest gas, with nothing fitted, and the one interval with no star formation falls fifty times below that ceiling.</p>
  </div>
</header>

<section>
  <div class="eyebrow">before the plots · what is measured</div>
  <dl class="defs">
    <dt>patch</dt><dd>a 0.5 × 0.5 kpc column of gas, in the frame of each galaxy (one frame after coalescence), 61 200 patch-snapshots in total; only patches with more than 1 M<sub>☉</sub> pc<sup>−2</sup> of gas and no contamination by the companion are used</dd>
    <dt>𝒲</dt><dd>weight: the column integral of gas density times vertical gravity, ∫ρ g<sub>z</sub> dz, over the gas within 0.5 kpc of the midplane; gravity is solved directly for gas, stars and dark matter, and is 90 to 96 per cent dark matter</dd>
    <dt>P<sub>tot</sub></dt><dd>midplane pressure of gas below 2×10<sup>4</sup> K: thermal, plus vertical turbulent ρσ<sub>z</sub><sup>2</sup>, plus the vertical Maxwell stress</dd>
    <dt>Σ<sub>SFR</sub>, Υ</dt><dd>star formation rate per area from stars younger than 10 Myr; Υ is the pressure per unit star formation rate calibrated by Ostriker &amp; Kim (2022)</dd>
    <dt>S, H, Σ, ε<sub>S</sub></dt><dd>shear rate of the mean flow across a patch, gas scale height, gas surface density; ε<sub>S</sub> is the dimensionless efficiency with which shear work feeds the turbulence (not a sound speed)</dd>
    <dt>M<sub>young</sub>, M<sub>max</sub></dt><dd>stars formed in a patch in the last 10 Myr; its most massive bound cluster (friends-of-friends, energy-bound)</dd>
    <dt>A, τ, f<sub>max</sub></dt><dd>patch area 0.25 kpc<sup>2</sup>; cluster formation window 10 Myr; fraction of an event in its largest cluster, 0.5</dd>
  </dl>
</section>
'''+body+'''
<section>
  <div class="eyebrow">the chain, and what kind of statement each link is</div>
  <div class="chain">
    <div class="row"><span>1</span><p>Layer in vertical equilibrium in every phase but coalescence — <b>measured</b>; the full-column failure is debris geometry.</p></div>
    <div class="row"><span>2</span><p>Pressure made by feedback only in the disc and in bursts — <b>measured</b>; PRFM regulation fails between bursts.</p></div>
    <div class="row"><span>3</span><p>Between bursts the support is shear-driven turbulence — <b>empirical closure</b>, one constant ε<sub>S</sub> = 0.02 in all phases.</p></div>
    <div class="row"><span>4</span><p>… plus a dynamo field at equipartition in the remnant — <b>measured</b>, not modelled.</p></div>
    <div class="row"><span>5</span><p>Cluster mass function: fixed slope, rising truncation — <b>measured</b>; the slope has no theory here.</p></div>
    <div class="row"><span>6</span><p>Largest cluster = half its patch’s event — <b>measured</b>, independent of environment.</p></div>
    <div class="row"><span>7</span><p>Largest cluster = ½ · Σ<sub>SFR</sub>(𝒲) · A · τ at the heaviest star-forming weight — <b>derived from 2 and 6, no free parameter</b>, matches to 0.06 dex.</p></div>
  </div>
  <p class="note">Caveats: one run, one resolution; the low-mass cluster slope in the merger phases needs the particle-number convergence redone per phase; the disc-phase field saturation is not established; the ceiling rests on seven intervals in which both quantities rise through the merger, though the two non-monotonic intervals follow it. That the largest cluster in solar masses roughly equals the weight in K cm<sup>−3</sup> is a coincidence of units.</p>
</section>
<div class="meta"><span>notes <b>paper/notes_2026-09-13_prfm_merger.md</b></span><span>narrative <b>paper/story_prfm_to_clusters.md</b></span><span>figures <b>prfm/story_figs.py</b></span><span>page <b>prfm/build_story_page.py</b></span><span>updated <b>2026-09-14</b></span></div>
</main>
'''
out='/tmp/claude-28632/-raven-u-uli-dwarf-merger/132f9852-368d-4946-8cb1-7274f158829e/scratchpad/prfm_story.html'; open(out,'w').write(html); print(out, len(html)//1024,'KB')
