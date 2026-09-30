# Module 1: verification against `dicty_reactions.jl` and the cited papers

## Reactions and constants (Julia, lines 521-889)

| Item | tex | Julia | Result |
|---|---|---|---|
| 1.1-1.4, 1.5, 1.5b, 1.6, 1.7, 1.8a/b | reaction table | L705-727 | ok (same educts/products and rate constants) |
| R_tot | 3077 /voxel = 40 001 /cell = 0.64 uM | L560-561; 3077/4817.6 = 0.639 | ok |
| k_on, k_off | 7.5, 0.45 (Kd 60 nM) | L577-578 | ok; equal to Table II of Van Haastert 1984 (H: k1 = 7.5, k-1 = 0.45) and Table 1 rows 1/2 of Biswas |
| eps_p, k_on,p, k_off,p | 0.125, 0.938, 0.45 (Kd,p 480 nM) | L601-603 | ok; the code comment (L590-598) gives the real reason for 1/8: 480 nM puts RpC at 30 % occupancy at the simulated wave peak (204 nM at the sender). Now stated in the table |
| k_ph, k_ph0, k_dp, k_b | 0.02, 5e-4, 5e-3, 5.78e-3 | L699-702 | ok; t1/2 = 34.7 s, 139 s, 120 s |
| k_A, k_A,off, [A]tot | 0.07, 0.02, 154 | L721-723, L1093 | ok; adapter is inert because k_acaOn = 0 (L2588-2617) |
| D_R, D_A | 0.024, 10 | L1036-1043 | ok |
| Optional Biswas table (rows 1-26) | 26 constants | L784-817 and Biswas Table 1 (parsed from the PDF text) | ok, every value equals the paper's Table 1 |

## Statements against the papers

| Statement | Paper | Result |
|---|---|---|
| 4 % slow (Kd 12.5 nM), 40 % fast high (60 nM), 60 % fast low (450 nM), conversion first order t1/2 ~ 9 s | Van Haastert 1984, abstract | ok (the three fractions sum to 104 % in the paper) |
| Kd ratio of the two fast classes 7.5 | Van Haastert 1984, Discussion p. 8: "a = 7.5" | ok |
| both fast classes dissociate with t1/2 ~ 1 s | Van Haastert 1984, abstract (H 1.5 s, L 0.7 s) | ok |
| phosphorylation t1/2 45 s at saturating cAMP, basal ligand-independent component, dephosphorylation t1/2 2 min (22 degC) | Vaughan 1988, abstract and methods | ok |
| serine clusters | Hereld 1994, abstract | ok |
| loss of ligand binding correlated with phosphorylation; affinity falls, receptor stays on the surface; 3-5-fold | Caterina 1995a/b | ok |
| ">80 % loss of high-affinity sites" [Johnson, Van Haastert] | Johnson 1991: >80 % loss of surface binding **sites**; Van Haastert 1984: high-affinity fraction falls from 40 % to about 10 % | **corrected** in text and constants table (two separate statements) |
| about 7e4 binding sites at 4 h [Johnson] | Johnson 1991 Fig. 3 and text | ok. Note: Laub 1998 (Module 0/8-10 papers) quotes about 40 000 cAR1 per cell at 4 h, and the code comment cites Saxe 1996 / Pupillo 1992 for 40 000; the tex calls 40 001 a design choice. Both statements can stand, the literature is not consistent |
| k_ph basis "L" | measured t1/2 45 s (k = 0.0154 /s), model 35 s (0.02 /s) | **changed to D**; text says "of the order of" |
| D = 0.024 (n = 27), 0.027 | Takebayashi 2023; Ueda 2001 (2.7e-10 cm2/s = 0.027 um2/s); Biswas Table 1 caption | ok |
| phospho-deficient receptors still adapt | Kim 1997, abstract | ok |
| fold-change response of the relay | Kamino 2017, abstract | ok |
| Biswas reduced model rows 27/28 | Biswas Table 1: R + cAMP <-> R:C with k = 0.885, k- = 0.524 | topology as stated; the constants of the model come from rows 1/2, not from rows 27/28 |
