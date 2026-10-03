# Module 3: verification against `dicty_reactions.jl` and the cited papers

## Julia (lines 1146-1518; 3.3b in L2041-2046, brake in L2185-2239)

| Item | tex | Julia | Result |
|---|---|---|---|
| RasG, GEF_R, GAP, brake pools | 38 500, 138 600, 46 200, 1540 /voxel | L1149-1203 | ok |
| k_E,on, k_E,off, k_R,on | 3.8e-3, 0.5, 0.333 | L1204-1206 | ok (tau_E = 2 s) |
| k_I,on (3.4e), k_I,off, k_I,b | 2.53e-4, 0.033 (tau_I 30 s), 5e-4 | L1207-1216, L1276 | ok |
| K_NF1 | 0.2 uM | L1436 | ok |
| k_cat | 15.8 /s, 11 x 1.4 | L1443: recomputed from k_rasGon/M x (E/I)_rest x RasG/theta with Gbg_rest = 6/voxel (from the 2.0/2.3/2.5 balance), E/I = 2.974: 15.84 | ok |
| k_PIP3 | 0.04 | 40/PIP2_SCALE (1000.3) = 0.03999 | ok |
| k_B,on, k_B,off, k_B,hyd | 24.1, 0.1, 39.8 | 0.1 x 4817.6/20 = 24.09; 0.6 x 10.6/((1540/2)/4817.6) = 39.79 | ok |
| D | 0.1 (GEF_R, RasG), 20 (GAP, B) | L1495-1514 | ok, incl. the caveat that B* is a membrane species with D = 20 |
| reaction 3.4 (RC-driven) | not in the table | defined (L1310) but used only for DICTY_M3_INHIB=rc | ok, default is gbg |
| the brake is read by PKBR1* (3.3c) | yes | L2236 (default since 2026-09-30) | ok |

## Papers

| Statement | Paper | Result |
|---|---|---|
| Ras rapid, transient activation by cAMP | Kae 2004, abstract; RasC peak at 5 s, back by 45 s; RasG within 5 s | ok |
| RasG upstream activator of PI3K | Sasaki 2004, abstract; Takeda 2012, p. 1 | ok |
| near-perfect adaptation, only IFFL fits, proportional activation, RasGAP as global inhibitor | Takeda 2012, abstract, p. 2-3 | ok |
| integral control gives oscillations and step-dependent kinetics | Takeda 2012, p. 3 | ok |
| RasGEF activation faster than RasGAP | Takeda 2012, p. 3 | ok |
| Ras activator/inactivator (RasGEF, RasGAP) both driven by receptor/G protein | Takeda 2012, Fig. 3A legend, PDF p. 4 ("leading to the activation of both RasGEF and RasGAP, which activates and inactivates Ras"); Meier-Schellersheim 2006, PDF p. 12 ("direct activation of Ras through Gβγ. In addition, Gβγ activates RasGAP"; modelling paper) | ok; Levchenko 2002 / Ma 2004 describe LEGI for PI3K/PTEN, not Ras, and are cited only for the LEGI concept (unmarked) |
| fold-change detection | Goentoro 2009, abstract | ok |
| DdNF1 major Ras regulator, unregulated Ras in nfaA- | Zhang 2008, abstract | ok |
| rectification with zero-order ultrasensitivity | Nakajima 2014, abstract | ok |
| NF1 KM 0.3 uM, kcat 1.4 /s | Wiesmuller 1992 Table I, row NF1-333 (kcat 1.4 /s, KM 0.3 uM); Coyle 2016 p. 11 | ok |
| enhanced basal Ras activity | Zhang 2008 (nfaA-: "elevated basal levels"); Xu 2021 (c2gapA-: "enhanced basal Ras activity") | ok |
| RasGEFR activates RasG; Gbeta needed for the initial response | Kae 2007, abstract; Kortholt 2013, abstract | ok |
| PIP3 -> Ras positive feedback | Fukushima 2019, abstract ("positive-feedback regulation of Ras-GTP by the downstream PIP3"); Sasaki 2004 shows only Ras -> PI3K | ok, Sasaki supports the Ras -> PI3K half |
| amplification between G protein and Ras | Kataria 2013, abstract | ok |
| pkbA-/pkbr1-: RasC increased, fails to adapt by 40 s; PKB/PKBR1 phosphorylate Sca1 | Charest 2010, p. 8 and abstract | ok |
| PKBs regulate other Ras and PI3K; delayed negative feedback | Miao 2017, p. 3 and p. 11 | ok |
| "mechanism for the RasG reporter is not identified" | Miao 2017, p. 11 only says PKBs "indirectly inhibit Ras GEF Aimless and directly activate PI(5)K" (previous reports) | supported only indirectly (mechanisms are proposed, not established) |
| excitable, fires without cues; stimulus biases threshold | Fukushima 2019, Miao 2017, Devreotes 2017 (p. 1, p. 12 "biased excitable network"), Iglesias 2012 abstract | ok |

## 2026-10-02: PIP3 feedback and brake re-labelled

* 3.3b (PIP3 -> RasG-GTP) replaced by 3.3g (PIP3 -> G-protein dissociation -> Gβγ) as default (`DICTY_M3_PIP3PATH=gbg`; `ras` restores 3.3b). The placement is a proposal; the sources support only the loop as a whole (Sasaki 2004, Fukushima 2019) and an amplification between G protein and Ras (Kataria 2013). `k_pip3Gbg` is untuned.
* 3.3c/d/e labelled *proposed*: Charest 2010 and Cai 2010 concern RasC (Module 8, 8.2s); Miao 2017 only names the PKBs as candidates (RBD is isoform-unspecific); Takeda 2012 rejects an integral-type feedback. No source documents a PKB/PKBR1 or PIP3 brake of RasG. Decision rule for keeping it: see the 3.3c/d/e block in `dicty_reactions.jl`.
* Marks M3-ch-adapt, M3-sca1, M3-miao-ext, M3-miao-mech are kept but the text now states what each supports (RasC, candidate hint, mechanism unknown).

## 2026-10-02: sources added for the PIP3 -> G-protein question

* New: `vanhemert2010` (Van Hemert et al. 2010, J Cell Sci 123:2922; library, `M3/`, INDEX, Bibliography). `elzie2009` is now also used in M3. The second copy `2597.pdf` was identical to `all/elzie2009.pdf` and moved to `sources/_duplicates/`.
* New claims: g_persist (Elzie 2009 p. 2), gbg_pi3k_indep (Van Hemert 2010 p. 6; mobility, not dissociation), sb_pip3_indep (Kortholt 2013 p. 1), ras_pip3_indep (Fukushima 2019 p. 2). Marked PDFs regenerated.
* Result: no source measures a PIP3/PI3K effect on the G protein; the four statements point against PIP3 -> Gβγ (3.3g) and do not place the measured PIP3 -> Ras effect.

## 2026-10-02 (evening): adaptation vs. brake, PIP3 feedback back on Ras, tau_I = 10 s

Changes in `dicty_reactions.jl` (no simulation was run; the module is not re-calibrated):

| Item | tex | Julia | Result |
|---|---|---|---|
| PIP3 feedback default | 3.3b, k_PIP3,R = 0.04 | `M3_PIP3PATH` default `ras`; `k_pip3Ras` = 40/PIP2_SCALE = 0.03999 | ok; 3.3g (`gbg`) kept as option |
| k_I,off | 0.1 /s (tau_I 10 s) | `k_gapOff` 0.1 | ok (Takeda 2012 Supp. Table S1, p. 22) |
| k_I,on (3.4e), k_I,b | 7.67e-4, 1.52e-3 | 0.000253 x 3.0303 = 7.667e-4; 5e-4 x 3.0303 = 1.515e-3 | ok; rc-path k_gapOn 0.08 -> 0.2424 likewise (not in the table, not default) |
| k_cat | 15.8 /s | `_fI_rest_G` depends on k_gapOnG/k_gapOff only -> unchanged | ok |
| brake reader | 3.3c_R, RasG-GTP | `M3_BRAKE_READ` default `ras`, `r3_3cR` (renamed from `r3_3cG` 2026-10-03; PKB readers 3.3c′ now `r3_3cP1`, `r3_3cPA`) | ok; `pkb`, `pip3` remain |
| k_B,on | 0.0753 | 0.1 x 4817.6 / 6400 = 0.07528 | ok (`BRAKE_HALF_RAS` = 6400, design anchor) |
| k_B,off, k_B,hyd, B_tot, D | 0.1, 39.8, 1540, 20 | unchanged | ok |
| Sca1 (M8) | 8.2p/8.2q/8.2s with S | `r8_2p`, `r8_2pA`, `r8_2q`, `r8_2s` on `Sca1_cyto`/`Sca1a_mem`; k_S,on = 24.09, k_S,off = 0.1, S_tot = 1540, D = 20 | ok; gated only on `DICTY_M8_SCA1FB` > 0 |
| Julia syntax | | `Meta.parseall`: 0 errors | the file cannot be loaded standalone here (engine and `cell_steady_state.jl` missing) |

New statements (claim, paper, page in `sources/M3/<key>.pdf`), all located by `mark_pdfs.py` without warning:

| Claim | Paper, page |
|---|---|
| g_persist (text corrected: Elzie did not measure PI3K, "would have subsided") | janetopoulos2001 p. 1; elzie2009 p. 2 |
| g_steps, reset | xu2022 p. 3 |
| nfblb_def | xu2022 p. 1 |
| c2gap_ras | xu2022 p. 2 |
| nfblb_fit | xu2022 p. 2, p. 4 |
| ras_no_down | kortholt2011 p. 1, p. 2 |
| excit_intr | fukushima2019 p. 2, p. 7; iwamoto2025 p. 2 |
| dnf | fukushima2019 p. 1 |
| pip3_gef | fukushima2019 p. 10 (Eqn 11 text) |
| refr | artemenko2016 p. 4 (text), p. 5 (Fig. 3 legend: LimE, paired 2-s mechanical stimuli) |
| nf1_rasg | zhang2008 p. 1 |
| undershoot | takeda2012 p. 2 |
| tau_fit | takeda2012_supp p. 22 (Table S1, rectangle), p. 6 (two passages) |
| tau_pair | takeda2012_supp p. 5 |
| int_def | takeda2012_supp p. 3 |
| slow_buf | takeda2012_supp p. 9 |
| gap_gbg | meierschellersheim2006 p. 5 |
| pip3_small | li2018 p. 4 |
| li_gap | li2018 p. 1 |
| ly_red | sasaki2004 p. 7 |
| gef_pfb | iwamoto2025 p. 10 |
| rbd | kae2004 p. 1; iwamoto2025 p. 2 |
| sca1_rasc | charest2010 p. 4 |

Checks: every `\mk/\mkt/\mkn` claim of sections/m3.tex and its tables exists in `marks_claims.tex`; every claim of m3.json is used in the TeX and has at least one mark; every `\cite` key exists in `Bibliography.bib` (new: xu2022, kortholt2011, artemenko2016, iwamoto2025, meierschellersheim2006). LaTeX (scratch copy, `alpha` style because `labelalpha.bst` is not installed here) compiles; the only error is the empty bibliography, which follows from the `\cite` redefinition in Main.tex.

Not in the paper (no PDF in the repo, so not markable): Huang 2013 (refractory period of the RBD signal), Sasaki 2007 (G-protein-independent Ras/PI3K loop), Edwards 2018, Lee 2010, Xu 2017, Xu 2026; see `sources/M3/SYNTHESE_Rueckkopplung.md`.

Internal check, not a source and not in the paper: Takeda's own IFFL (Supp. Text S2, Table S1) with two 2-s pulses gives a second/first peak ratio of 0.48 at 9-12 s and 0.85 at 45 s, i.e. relative but no absolute refractoriness; with tau_GAP = 30 s the step response falls below 5 % of the peak only after about 68 s (measured: < 35 s).

Open: re-measure the driven RasG-GTP and set `BRAKE_HALF_RAS`; check rest state (low branch), step return < 35 s, no post-peak oscillation, reset < 1 min; re-check Modules 0/4/8 after tau_I = 10 s.
