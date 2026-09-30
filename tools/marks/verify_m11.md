# Module 11: verification against `dicty_reactions.jl` and the cited papers

## Julia (lines 3147-3330)

| Item | tex | Julia | Result |
|---|---|---|---|
| k_pdsa | 0.122 uM^-1 s^-1 (1.2e5 /M/s) | min_to_s(4.9) x 1.5 = 0.1225 | ok |
| k_bas, k_ind | 4.48e-3, 4e-2 | L3187-3188 | ok |
| H_P, n_H | 25, 2 | H_DEV_HALF, N_DEV | ok |
| k_sec, k_shed, k_deg | 0.01, 5e-4, 5e-4 | L3189-3191 | ok (tau_deg = 2000 s) |
| k_pdi, K_pdi, n_pdi | 2.5e-3, 0.05, 2 | L3236-3238 | ok |
| k_trk | 4e-4 (tau = 42 min) | 1/4e-4 = 2500 s = 41.7 min | ok |
| k_i,on, k_i,off, K_d | 1, 1, 1 uM | L3250-3252 | ok |
| k_res | 1.22e-4 = k_pdsa/1000 | k_pdsaResidual | ok |
| D | cAMP 350, PdsA_e and complex 65, PdiA 60, PdsA_i 10, PdsA_m 0.1, cAMP_s 0 | L3779-3781, L3806-3817 | ok |
| G = 16, k_eff 9.6e-4 -> 0.0154 /s, lambda 600 -> 151 um, clearance 65 s | as in the text | measured values quoted in L3181-3196 | ok |
| sqrt(D_PdsA/k_deg) = 361 um | 361 | sqrt(65/5e-4) = 360.6 | ok |
| surface pool 2.5e6 molecules per cell | as in the text | 13 x 0.042127 uM/s x 4817.6 x f_sec x k_sec/(k_shed+k_deg) = 2.51e6 | ok |
| fraction of complexed enzyme "4 % at rho = 7.7e-4, 70 % at 6e-2" | tex | `m11_world_steady` (L3313-3327) gives 3.3 % and 11.3 %; the inhibitor:enzyme synthesis ratio is k_pdi P/((k_bas + k_ind P) f_sec f_shed) = 0.117, so more than 11.7 % of the enzyme can never be complexed | **corrected** to 3 % and 11 % with the upper bound 0.12 |

## Papers

| Statement | Paper | Result |
|---|---|---|
| secreted PDE controls the cAMP level; PDE-null cells cannot move coordinatedly in mounds | Sucgang 1997, abstract | ok |
| extracellular cAMP degraded predominantly by DdPDE1 and DdPDE7; DdPDE1 Km | Bader 2007, abstract and kinetic table (Km 0.75 uM, Hill 0.8) | ok (0.75 rounded to 0.8) |
| pdsA- cells: narrower chemotactic range, no streams at low density | Garcia 2009, abstract | ok |
| three promoters (growth, aggregation, late) | Faure 1990, abstract | ok |
| aggregative promoter down-regulated after aggregation | Weening 2003, abstract | ok |
| inhibitor 1:1, Kd about 1e-10 M, absent in growth, appears on starvation | Franke 1981, abstract (Kd 9e-11 to 1.1e-10 M) | ok |
| inhibitor synthesis repressed by cAMP pulses, independent regulation | Yeh 1978, abstract | ok |
| D(cAMP) = 4.44e-6 cm2/s | Dworkin 1977, abstract | ok |
