# Module 9: verification against `dicty_reactions.jl` and the cited papers

## Julia (lines 3022-3078)

| Item | tex | Julia | Result |
|---|---|---|---|
| k_9, k_11, k_9+k_11 | 8, 12, 20 /s | k9_eff, k11_eff (L3040-3041) | ok |
| alpha_H | 0.5 | alpha_H (L98), used in r9_1 | ok |
| k_x | 0.5 /s | k_abcB3 | ok |
| k_2' | 0.0469 | min_to_s(0.9)/(1540/4817.6) = 0.015/0.3197 = 0.0469 | ok |
| k_Ao | 6.4e-3 /s, t1/2 = 108 s | k_acaOff = 0.0064; ln2/0.0064 = 108.3 s | ok |
| ACA pool | 865 /voxel | ACA_v (L516) | ok |
| exported molecules per pulse | about 9e6 | 9.32e6 in the acceptance note (L1919) | ok |
| exported fraction 0.05 at rest, 0.83 at the peak; factor 250; factor 1.0 | as in the text | comment L3049-3058 | ok |
| D | cAMP_i 30, ACA 0.1 | L3732, L3018 | ok |

## Papers

| Statement | Paper | Result |
|---|---|---|
| aca- cells: little cyclase activity, no aggregation; motility and chemotaxis unaffected | Pitt 1992, summary | ok |
| ACA activation requires CRAC | Insall 1994, abstract | ok |
| ERK2 and ACA | Segall 1995 (necessary for receptor-mediated activation) | ok as a report of Segall 1995; wording changed to "was reported to be necessary" with a pointer to Module 0, because Maeda 2004 later showed ERK2 is not essential for ACA activation (see verify_m0.md) |
| ACA peaks at 30 s, decline half-time 1.8 min | Cai 2010, p. 5 | ok |
| AbcB3 is a component of the cAMP export | Miranda 2015, abstract | ok |
| secretion proportional to intracellular cAMP; K_S = 0.34 and 0.94 /min; K_P = 1.73 /min; 16 % (1e-6 M) and 47 % (1e-8 M) secreted | Dinauer 1980, abstract, p. 6, p. 7 | ok |
| about 5 pmol per 1e6 cells (2-min stimulus, 7 h) = 3e6 molecules | Devreotes 1979, abstract; 5e-12 x 6.022e23/1e6 = 3.0e6 | ok |
| PKA-null: ACA and ERK2 rise and stay high; k2 = 0.9 /min | Laub 1998, p. 3 and p. 4 | ok (k2 is the ACA inactivation constant) |
