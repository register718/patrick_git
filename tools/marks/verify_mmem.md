# Optional memory module: verification against `dicty_reactions.jl` and the cited paper

| Item | tex | Julia (L3607-3690) | Result |
|---|---|---|---|
| [Mem]tot | 5000 /voxel | MEM_V = 5000 | ok |
| k_w | 4e-3 /s | k_memOn = 0.004 | ok |
| K_mem, n_mem | 5000, 3 | K_MEM = 5000, N_MEM = 3; H(1400) = 0.021, H(1600) = 0.032, H(7500) = 0.771 (Table I: 0.03, 0.77) | ok |
| k_e | 5.56e-3 /s (tau = 3 min) | k_memOff = 1/180 | ok |
| D | 10, 0.01 | MEM_DIFFUSION | ok |
| module off by default | yes | MEMORY_ON = (DICTY_MEMORY == "1"), default 0 | ok |

| Statement | Paper | Result |
|---|---|---|
| cells keep their direction in the back of a wave for the natural 6 min period, reverse for longer periods | Skoge 2014, abstract | ok |
| the molecular carrier is not known | Skoge 2014, p. 4 ("the molecular basis of cellular memory remains to be determined") | ok |
| erasure time of 3 min | the code cites Skoge 2014 (about 2 min in the back of 6-10 min waves, at least 5 min after a reversal) and Shi 2013 (not in the library) | a bracket, not a measurement; the constant is tagged E, ok |
