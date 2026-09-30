# Module 10: verification against `dicty_reactions.jl` and the cited papers

## Julia (lines 3079-3146)

| Item | tex | Julia | Result |
|---|---|---|---|
| pools RegA, PKA, ERK2 | 1540 /voxel (0.32 uM) | RegA_pool_v, PKA_pool_v, ERK2_pool_v | ok |
| k_B1, k_B2 (resting active fraction 0.9) | 0.05, 5.6e-3 | 0.05/(0.05 + 0.0056) = 0.899 | ok |
| k_B3 | 8 (tau ~0.5 s at the ERK2* peak, knock-down 64x; literal Laub-Loomis value: 1.3x) | L3103-3115 comments | ok |
| k_B4 | 0.0167 (tau = 60 s) | k_regAReact | ok; the code comment says "LL k7 = 1.0 /min", the table of Laub 1998 gives k7 = 2.0 /min mM (k10 = 1.0), the tex quotes 2 /min |
| k_cat | 30 uM^-1 s^-1 = 3e7 /M/s, k_cat = 150 /s at KM = 5 uM | k_regACat | ok |
| phi | 0.05 | REGA_U_FRAC | ok |
| k_D1, k_D2 (Kd 0.1 uM, tau 40 s) | 0.25, 0.025 | k_pkaOn, k_pkaOff | ok |
| k_D3 | 0.0417 | min_to_s(0.8)/(1540/4817.6) = 0.01333/0.3197 = 0.0417 | ok |
| reactions 10.B1-D3 | as in the table | r10_B1 ... r10_D3 | ok |

## Papers

| Statement | Paper | Result |
|---|---|---|
| RegA = DdPDE2 degrades intracellular cAMP | Bader 2007, abstract | ok |
| Asp212 phosphorylation by RdeA, at least 20-fold activation | Thomason 1999, abstract (Thomason 1998 reports "up to 8-fold" with a heterologous donor) | ok; the 20-fold figure is from 1999 |
| Km about 5 uM | Thomason 1998, abstract | ok |
| loss of RegA raises PKA activity | Shaulsky 1998, summary ("inhibition of the phosphodiesterase results in an increase in the activity of PKA") | ok |
| cAMP -> PKA -| ERK2 -> RegA circuit; ERK2 in phase with cAMP | Laub 1998 abstract; Maeda 2004 abstract | ok |
| PKA shapes ERK2 activation and adaptation | Aubry 1997 abstract (pka null: lower, more extended ERK2 activation) | ok |
| ERK2 inhibits RegA | Maeda 2004 p. 3 ("ERK2 may directly inhibit phosphodiesterase activity by phosphorylation of RegA") | ok |
| Laub-Loomis k4 = 1.5 /min, k6 = 0.8, k7 = 2.0, k8 = 1.3 | Laub 1998 Table 2 (p. 4) | ok; k7 is the RegA production coefficient, the text now calls it "replenishment coefficient" instead of "recovery coefficient" |
