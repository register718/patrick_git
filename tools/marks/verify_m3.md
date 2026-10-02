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
| LEGI | Levchenko 2002 abstract ("coordinately controlled by the G-protein activation"); Ma 2004 abstract ("controlled by receptor occupancy") | ok; the two papers name different control variables (G protein / receptor occupancy), the tex follows Levchenko |
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
