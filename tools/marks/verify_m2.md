# Module 2: verification against `dicty_reactions.jl` and the cited papers

## Julia (lines 889-1143)

| Item | tex | Julia | Result |
|---|---|---|---|
| G:R, G_tot | 60, 184 620 /voxel (38 uM) | L893-894: 3077 x 60 = 184 620; /4817.6 = 38.3 uM | ok (comment "18462 = 3.83 uM" in the code belongs to the old ratio 6) |
| t1/2,diss, k_E | 3.2 s, 0.339 /uM/s | L898, L915-918: ln2/(3.2 x 0.6387) = 0.339 | ok |
| eps_E | 1 | L930 | ok (the comment at r2_2, L1001, still says "5 %": stale) |
| k_E0 | 2.18e-7 /s | L936 | ok |
| k_hyd | 0.162 = 7 ln2/30 | L985-986 | ok |
| k_re | 6.09 = 3500/(30 x [G]tot/2) | L992-995 | ok |
| D | 0.2 for all three forms | L1020-1035 | ok |
| Reaction 2.4 (RGS) | "defined but switched off" | L1006 commented out | ok |
| contrast 0.22 -> 0.77 (60 nM), SNR doubled; 0.85 -> 1.65 (5 nM); dissociation 3.0 -> 2.0 s, reassociation 33 -> 7.5 s; cost x3 | as in tex | L939-981 (0.224 -> 0.766, SNR 15.3 -> 31.7; 0.853 -> 1.654, SNR 34.8 -> 30.9) | ok; **SNR wording at 5 nM corrected** ("no loss" -> 34.8 to 30.9, i.e. -11 %) |

## Papers

| Statement | Paper | Result |
|---|---|---|
| FRET: rapid dissociation/reassociation; steady state, no decline; catalysis whether or not phosphorylated | Janetopoulos 2001, abstract (all three verbatim) | ok |
| FRET-loss half-time 3.4 +/- 1.4 s; reassociation 32.3 +/- 5.3 and 28.7 +/- 6.2 s (about 30 s) | Tang 2014, Fig. 4b and text | ok |
| Biswas reproduce 3.0-3.3 s; k_E0 = 2.18e-7 at ratio 3; ratio 3; D = 0.2; all occupied classes weighted equally | Biswas Table 2 (row t1/2,diss 3.00-3.29 s, ratio-3 column 2.18), p. 8 (ratio 3), Table 3 caption | ok |
| initial Ras response requires Gbeta not Galpha2 | Kortholt 2013, abstract | ok |
| symmetry breaking and amplification between G protein and Ras | Kataria 2013, abstract | ok |
| G proteins cycle between cytosol and membrane; equal FRAP kinetics of the three subunits | Elzie 2009, abstract of results and Fig. 3 | ok. Elzie gives no diffusion constant, so 0.2 uM2/s stays an estimate (as tagged E) |
