# Module 0: verification of tables and text against `dicty_reactions.jl`

Checked on the branch state after merging `main` (file `dicty_reactions.jl`, 3989 lines). "ok" = tex value equals the
default value in the Julia code (ENV knobs unset) or follows from it arithmetically.

## Reactions (`tab_m0_rxn.tex`)

| ID | Julia | Result |
|---|---|---|
| 0.8a / 0.8b | L100-103 `k_H_prod = H_max*k_H_deg/M`, `k_H_deg` | ok (zeroth order `k_H H_max`, first order `k_H N_H`) |
| 0.9a | L439-441, gate `Hunger^4/(HC_N+Hunger^4) * KR_N/(KR_N+Refr^4)` | ok |
| 0.9b | L445, `[X]->[X,X]`, rate `13*k1/2`, factor `X-1` | ok (net 2X->3X, propensity k1/2 N(N-1)) |
| 0.9c | L449, `[X]->0`, rate `13^2*k6/6`, factor `(X-1)(X-2)` | ok (net 3X->2X) |
| 0.9d | L452 | ok |
| 0.9y1 / y2 / y3 | L454, L457, L459 | ok |
| 0.9r1 / r2 | L430, L433 | ok |
| 0.9e | L497-499, `[X,ERK2]->[X,ERK2a]`, factor `X-1` | ok (written 2X+ERK2 -> 2X+ERK2* with falling factorial) |
| 0.9e' | L501 | ok |
| 0.9f | L517-519, gate `ERK2a^4/(KE^4+ERK2a^4)` | ok |

## Constants (`tab_m0_const.tex`)

| Symbol | tex | Julia | Result |
|---|---|---|---|
| H_max | 100 | L95 | ok |
| k_H | 2.8e-4 (tau ~ 1 h) | L96; 1/2.8e-4 = 3571 s = 0.99 h | ok |
| H_1/2, n_c | 60, 4 | L159-160; c(100) = 0.8853 | ok |
| k0 | 3.5 | L257 | ok; wording corrected: "at full competence" -> "before the competence gate" (c(H_max) = 0.885, effective seed 3.10 /s) |
| k1, k6, k2 | 0.08, 5e-4, 1 | L258-260 | ok |
| k3, k5, k4 | 6.25e-5, 0.05, 0.01 | L261-263 | ok |
| k_Rp, k_Rd | 5e-3, 5e-3 (tau_R = 200 s) | L356-357 | ok |
| K_R, n_R | 450, 4 | L421-422 | ok |
| k_XE | 1823 uM^-2 s^-1 | L482 | ok |
| k_Ed | 0.02 (tau = 50 s) | L483 | ok |
| k_XA | 65 uM^-1 s^-1 | L490 | ok |
| K_E, n_E | 192 (12.5 % of pool), 4 | L491-492; 192/1540 = 12.47 % | ok (the code comments L468, L492 still say "5000 molecules/cell" and "25 %": stale, they refer to the old pool) |
| N_vox | 13; M_c = 62 629 | L171-172 (62628.8) | ok |
| ERK2_tot | 1540 (0.32 uM) | L86; 1540/4817.6 = 0.320 | ok |
| ACA_tot | 865 | L516 | ok |
| D_X, D_Y, D_Rf, D_ERK2 | 10 | L3762, L3766-3771 | ok |

## Numbers in the text (`m0.tex`)

| Statement | Result |
|---|---|
| mean interval 726 s, 95 % interval 706-747 s, 8 seeds; 2.52 gave 930 s; block adds ~390 s; CV 0.16 | ok (L223-228, L241-247) |
| 1000 cells: 467 +/- 16 s (7.8 min) | ok (L252-256, L402-405) |
| recovery curve: none within 300 s, 1/3 at 360 s, all by 420 s, three seeds | ok (L385-390) |
| resting level / threshold / spike maximum | **corrected**. The code comment (L233-234: 3.97 / 23.21 / 455.82) solves the nullcline without the competence gate. With c(100) = 0.885 (L440) I get 3.43 / 23.78 / 455.79, so the text now says 3.4 / 24 / 456 (recomputed independently from the reactions 0.9a-d) |
| Hunger scales competence, ACA synthesis 9.1, PdsA induction | ok (L440, L3061, L3177/L3242). The header comment L91-92 also claims RegA scaling (10.1); Hunger does not occur there |
| ERK2*: first-order drive left 5 % of pool, brake 87 % off | ok (L473-479) |

## Statements against the papers (corrections)

| Statement | Paper says | Action |
|---|---|---|
| "ERK2 is required for ACA activation" [segall1995] | Segall 1995: ERK2 null cells show little cAMP synthesis, "important for receptor-mediated activation of adenylyl cyclase". Maeda 2004 (p. 3): "ERK2 is not essential for activation of ACA and RegA limits the accumulation of cAMP", circuit modified accordingly | sentence in 0.9f and the ACA-gate row of Table I now state both; the ERK2* gate on ACA is called a design choice |
| gamma distribution (shape 6.6, rate 1.25 /min) together with "about once per 12 min" [ford2023] | Ford 2023 Fig. 1E: fitted gamma has mean mu = 5.3-5.9 min (6.6/1.25 = 5.28 min, mode 4.5 min); the text of Ford quotes 1 pulse per 12 min | tex now gives the mean of 5.3 min and states that the two numbers differ and the model follows 12 min. **Decision for the authors**: which value is the target |
