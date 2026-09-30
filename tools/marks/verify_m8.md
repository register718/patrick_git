# Module 8: verification against `dicty_reactions.jl` and the cited papers

## Julia (lines 2243-2634 and 2900-3020)

| Item | tex | Julia | Result |
|---|---|---|---|
| [GEF_A]tot | 1600 | GEFA_v_pair = 154 x f_legacy/f_new = 1599.6 at 67 000 Gbg/voxel (L2243-2260) | ok |
| k_A,on, k_A,off | 7.6e-4, 0.1 | k_gefRon x 0.1/0.5 = 7.6e-4 | ok |
| k_S | 1 | k_sca1Fb = 1 | ok |
| k_C,on, k_C,off | 2, 0.02 | r8_3, r8_4 | ok |
| k_CG,on/off/hyd | 0.05, 0.005, 150 | r8_4b/c/d | ok |
| k_T,on, k_T,off | 50, 0.5 | L2280-2290 | ok |
| k_P,on, k_P,off | 30, 0.5 | L2668-2670 | ok |
| k_Cr | 10 | k_cracOn | ok |
| k_Cr,P, k_Cr,b, k_Cr,off | 1, 0.2, 0.1 | k_cracPIP3 = nM_to_uM(1e-3) = 1 uM^-1 s^-1 (assuming the helper converts 1e-3 /nM/s), k_cracBas, k_cracOff | ok |
| k_Cr,dp | 1.5 (a = 567, f = 0.44) | 10 x 567/4817.6 x 0.56/0.44 = 1.498 | ok |
| k_X, N_CRAC | 200 molec/s, 2002 -> 0.1 /s | k_m8X_max = 200, CRAC_TOT = 154 x 13 = 2002 | ok |
| K_P, n_P | 5000, 8 | L2797, L2920 | ok |
| k_PKA | 0.938 | 30 x (0.6/60)/(1540/4817.6) = 0.938 | value ok, **attribution corrected** (see below) |
| pools | RasC 385, CGAP/TORC2/CRAC 154, PKBR1 77 /voxel | L2993-3005 | ok |
| D | 0.1 membrane, 10 cytosolic, GEF_A* 2 | M8_DIFFUSION | ok |
| excluded edges 8.7/8.8, 8.13, 8.11 (k = 0) | as in the text | not in MODULE_8_REACTIONS; k_acaOn = 0 | ok |

## Papers

| Statement | Paper | Result |
|---|---|---|
| RasC ACA, RasG chemotaxis; RasGEFA activates RasC not RasG | Kae 2007, abstract | ok |
| "both are G-beta-gamma-driven RasGEFs" [kae2007] | Kae 2007 shows that G proteins activate separate Ras pathways through specific RasGEFs; it does not establish Gbeta-gamma specifically (Kortholt 2013 does that for the initial Ras response) | **wording changed to "G-protein-driven"** in the text and in the k_A,on row |
| RasC upstream of TORC2, PKBR1 independent of PIP3, prolongs ACA activation | Cai 2010 abstract (ACA half-time of decline 1.8 -> 5 min with active RasC) | ok |
| TORC2 activates PkbA and PKBR1 within seconds, also without PI3K activity | Kamimura 2008 abstract | ok |
| PKBR1 has no PH domain, is myristoylated, needs no PI3K | Meili 2000 abstract | ok |
| PKB/PKBR1 phosphorylate Sca1 | Charest 2010 abstract; RasC not adapted at 40 s (p. 8); Sca1 phosphorylation peak 5-10 s (p. 7) | ok |
| PKA feedback on Ras, Rap1, TORC2 | Scavello 2017 abstract | ok |
| CRAC required for ACA activation; binds PI3K products; translocation reporter | Insall 1994; Comer 2005; Parent 1998 abstracts | ok |
| PKBR1 HM peaks at 30-60 s, half-life 40 s, baseline after 2-3 min | Cai 2010 p. 2-3 | ok |
| other substrates peak at 20-60 s | Kamimura 2008 p. 2, p. 4 | ok |
| Laub-Loomis coefficient 0.6 /min as the PKA feedback coefficient | Laub 1998 p. 4: k5 = 0.6 /min is the CAR1 -> ERK2 activation coefficient, k6 = 0.8 /min the PKA-inhibition-of-ERK2 coefficient (the latter is the one used in Module 10, L3125) | **corrected**: the 0.6 /min of the PKBR1 feedback is not a PKA coefficient of Laub-Loomis; row retagged from D to E and the text states this |
