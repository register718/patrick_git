# Module 4: verification against `dicty_reactions.jl` and the cited papers

## Julia (lines 1518-2036, 4.14 at L1840-1872)

| Item | tex | Julia | Result |
|---|---|---|---|
| s4 | 20 on 4.1-4.8 only | M4_SPEED = 20 (L1518+), `_m4s` on 4.1-4.8, not on 4.10-4.14 | ok |
| pools PI3K, PIP2, PTEN, PkbA | 200 000, 769 231, 10 000, 100 000 /voxel | L1526-1527, L1875, L1935-1937 | ok |
| k_3K,on, k_3K,off | 0.4, 10 | 0.02 x 20, 0.5 x 20 | ok |
| k_syn | 0.0955 | 6.2 x 20/(200000/154) = 0.09548 | ok |
| k_PT,on, k_PT,off | 0.6, 4 | 30 x 20/1000.3 = 0.5998; 0.2 x 20 | ok |
| k_ev, K_ev, n_ev | 800, 1.31e5 (17 % of the lipid pool), 2 | 40 x 20; 0.17 x 769 231 = 130 769; N_DISP = 2 | ok |
| k_PT | 1.39 | 4.5 x 20/(10000/154) = 1.386 | ok |
| k_SH | 1.28 with 77 SHIP/voxel | 4 x 20 x 77/4817.6 = 1.2786 (k_ship_eff, L1657-1660) | ok |
| k_K,on | 1.32e-3 (rest 0.002 /s at 7298 PIP3) | 0.002/(7298/4817.6) = 1.320e-3 | ok |
| k_K,act | 0.507 (rest 8.5e-3 /s at 81 TORC2*) | 0.00853/(81/4817.6) = 0.5073 | ok |
| k_K,off, k_K,deact, k_K,fb | 0.05, 0.02, 160 | L1846-1855 | ok |
| reaction 4.14b (PKBR1* removes PI3K) | not in the table | comment only (L1857-1872); no `r4_14b` is defined | ok, the table is complete |
| reaction 4.9 (PI(3,4)P2 -> PIP2) | not in the table | collapsed into 4.8, PI34P2 has no producer | ok, as stated in the text |
| "flux about 1.7e3 molecules/voxel/s against a pool of about 6e4" | tex | not stated anywhere in the Julia file | not verifiable from the code; check against the calibration notes |

## Papers

| Statement | Paper | Result |
|---|---|---|
| PI3K at the front, PTEN reciprocal at the back | Funamoto 2002 abstract; Iijima 2002 (PTEN-GFP at the rear) | ok |
| bulk PIP3 transient, PI3K activation within 5 s then declines; larger and longer in pten- cells | Huang 2003 abstract (the PI3K profile is identical in wild type and pten-, the PIP3 changes are larger) | ok |
| CRAC translocation reflects activation of the G-protein system | Parent 1998 abstract | ok |
| PHCrac-GFP: maximum at 6-8 s, gone after 20 s, patches from 30 s | Postma 2003 abstract | ok |
| Ras activates PI3K | Sasaki 2004, Fukushima 2019 abstracts | ok |
| PIP3 suppresses PTEN, bistability | Matsuoka 2018 abstract | ok |
| PKB activation delayed, delayed negative feedback | Miao 2017 p. 3; Charest 2010 p. 8 (RasC not adapted in pkb- cells) | ok |
| PTEN removes the 3-phosphate of PIP3 | Iijima 2002 (PI 3-phosphatase PTEN), Huang 2003 | ok |
| PI(3,4)P2 part of the Ras excitable network | Li 2018 abstract | ok |
| PkbA transient and PI3K dependent; PkbA has a PH domain, PKBR1 has none | Meili 1999 abstract; Meili 2000 abstract | ok |
| TORC2 acts through the hydrophobic motif, at the front; PdkA/PdkB are the activation-loop kinases | Kamimura 2008 abstract; Kamimura 2010 abstract | ok. Kamimura 2008 shows the TORC2-PKB activation is PIP3-independent, the model keeps the PIP3-dependent recruitment of Meili 1999 and adds TORC2 on top (as described in the text) |
| Sca1/RasGEF/PP2A feedback | Charest 2010 abstract | ok |
| TORC2 activated by RasC | Cai 2010 (TORC2 binds active RasC; "RasC acts through TORC2") | ok |
