# ============================================================================
# ⚠ 2026-10-02 — MODULE 3 DEFAULTS CHANGED, NOT RE-CALIBRATED (no run was made):
#   * PIP3 feedback back on RasG (3.3b, DICTY_M3_PIP3PATH=ras); 3.3g (PIP3 → Gβγ) is an option.
#   * brake 3.3c read by RasG-GTP (3.3c_R, DICTY_M3_BRAKEREAD=ras; pkb|pip3 are options);
#     BRAKE_HALF_RAS = 6400/voxel is a design anchor to be re-measured.
#   * τ_I 30 → 10 s (Takeda 2012 Supp. Table S1), k_gapOn/k_gapOnG/k_gapBasal ×3.03 so the
#     activated GAP fractions are unchanged.
#   * Sca1 feedback on RasC uses its own species Sca1_cyto/Sca1a_mem (8.2p/8.2q/8.2s).
#   Evidence: sections/m3.tex, sources/M3/SYNTHESE_Rueckkopplung.md, tools/marks/verify_m3.md.
#   Check before use: rest state low branch, step return < 35 s, no post-peak oscillation.
# ============================================================================
# ⚠ STATE AS OF 2026-08-26 — read this before trusting any constant below.
# ============================================================================
# MEASURED ON THE ENGINE AND ACTED ON:
#   * The pre-2026-08-26 committed network FREE-RUNS.  calib/m8_probe.jl rest,
#     1800 s, no stimulus (2.41 spikes expected): with 3.3b+brake OFF, resting
#     PIP3 is 51 419/voxel, Module 8's PIP3 gate reads 1.00 (WIDE OPEN), the
#     relay drive is 110 molec/s against Module 0.9's own seed k0 = 2.52 (44x),
#     and the cell fires 5 times instead of 2.41.  Restoring 3.3b + the brake at
#     their ORIGINAL constants shuts the gate (1.9e-3) and returns the rate to 2.
#     This is the defect the standing K_M8X_PIP3 warning predicted.  3.3b is now
#     LIVE (M3_STEN); DICTY_M3_STEN=0 reverts.
#   * Module 4's own loop gain CANNOT exceed 0.918 over PTEN_v x{4..30} x
#     k_ship/{1..10} x N_DISP{2,3,4,6} x K_DISP — a structural ceiling.  Raising
#     N_DISP makes it WORSE.  The amplifier therefore has to be 3.3b, which is
#     also where Kataria 2013 and Fukushima 2019 put it.
#
# ⚠ MEASURED AND NOT YET FIXED — THE INITIAL CONDITIONS ARE WRONG:
#   calib/m34_rest.jl relaxes the cell from TWO starting points (committed ICs,
#   and the loop variables moved 10x) and compares.  Full cell, 600 s:
#       species          rest    IC        IC/rest
#       RasG_GTP_mem     5 810   69 628     12.0x
#       PI3K_mem         9 387  110 695     11.8x
#       PIP3_mem        48 395  400 946      8.3x
#       RasGEFRa_cyto    5 725      104      0.02x   (IC 55x TOO LOW)
#       Gbg_cyto        26 064      455      0.02x   (IC 57x TOO LOW)
#       RasBrakea_mem    1 905        0      seeds an inert species at 0
#       PKBA_mem/PKBAa_mem 32 565 / 12 715   both seeded at 0
#   i.e. the ICs are wrong in BOTH directions by one to two orders of magnitude.
#   Every driver starts each cell far off the manifold.
#
# ⚠⚠ AND THE FULL CELL HAS NO STATIC REST STATE, so "the initial condition" is
#   not a well-posed quantity as currently derived.  Within-run drift is 168-276 %
#   on RasG_GTP_mem / PIP3_mem / Gbg_cyto — NOT non-convergence: the cell PULSES
#   (Modules 0/9/10 fire every ~750 s), so those species traverse a limit cycle
#   and a window mean depends on where in the cycle it lands.  The A/B agreement
#   is 0.96-1.09, i.e. both starting points reach the SAME cycle.  A correct IC is
#   therefore the INTER-PULSE TROUGH at a defined phase, not a time average — and
#   m34_rest.jl must be extended to report the trough before any IC is rewritten.
#   Deriving one from the window mean would replace a 12x error with a phase error.
#
# ⚠ SCOPE NOTE ON THE 2026-08-25 MODULE 1/2 WORK BELOW (G_RATIO 60, M2_KRE_GAIN,
#   k_hyd = ln2/30, k_phos = 1.0, EPS_GEF = 1.0, M2_GEF_GAIN = 1.463):
#   the VIS deficit those changes were designed to fix was substantially a
#   MEASUREMENT ARTIFACT.  The 0.32/0.36/0.67 figures came from two_cells.jl runs
#   in which both cells ignited spontaneously and flooded their own fields.  With
#   fire_cells.jl (commanded sender, silent receiver) the UNCHANGED network scores
#   VIS(RasG_GTP_mem) = 1.57 / 1.11, bracketing the 1.51 historical benchmark.
#   Paired A/B, 5 seeds/arm, in flight: arm A (these changes) 1.04 ± 0.21 (n=4),
#   arm B (reverted) 1.34 ± 0.23 (n=2) — NO improvement, possibly worse, error
#   bars overlapping.  The upstream effects ARE real and measured (Gbg contrast
#   transfer 0.262 -> 0.330; Ga2GDP_mem 54 315 -> 617/voxel, i.e. the
#   pseudo-first-order recapture the linearisation was built to produce) — they
#   simply do not survive to Ras.  NOT reverted pending the remaining seeds;
#   DICTY_M2_KREGAIN=1 DICTY_M2_KHYD=0.05 DICTY_M1_KPHOS=0.02 DICTY_M2_EPSGEF=0.05
#   DICTY_M2_GEFGAIN=1 reverts all of it.
# ============================================================================
# dicty_reactions.jl
# ============================================================================
# Reaction network of the Dictyostelium discoideum cAMP signaling cell.
#
# Every reaction below is a 1:1 transcription of the corresponding numbered
# reaction in DICTY_CELL.md (the "Reaction." line of each subsection).  The
# section number from DICTY_CELL.md is given in the comment above each channel,
# together with the educt→product form and the literature rate constant.
#
# Two reaction lists are exported for the simulation script:
#   CELL_REACTIONS  — reactions that run *inside* the agent (Modules 1–9,
#                     10.1 RegA, 11).  Their species are agent-internal, except
#                     cAMP_ext which is a shared world species (sensed / secreted).
#   WORLD_REACTIONS — reactions that run *in the world* on agent-free voxels
#                     (Module 10.2–10.6: the extracellular PdsA phosphodiesterase
#                     system that degrades secreted cAMP).
#
# NOTE: this file must be `include`d *after* `using .AgentRD`, because it
#       constructs `Reaction` objects.
# ============================================================================

# ── Unit conversion helpers ────────────────────────────────────────────────
nM_to_µM(k) = Float32(k * 1000)
min_to_s(k) = Float32(k / 60)

const DX_µM_REF    = 2.0                                       
const MOLEC_PER_µM = 6.022e23 * (DX_µM_REF^3 * 1e-15) * 1e-6   
const PKA_pool_v   = 1540  
const ERK2_pool_v  = 1540   

# ============================================================================
# Module 0.8 — Abstract Hunger species (starvation heterogeneity) [NEW Rev 4]
# §0.8: ∅ →[k_H_prod] Hunger →[k_H_deg] ∅
# Hunger encodes per-cell starvation/developmental state; scales ACA synthesis
# (9.1) and RegA effective activity (10.1) so cells near H_max sit closer to
# the Hopf bifurcation.  No diffusion (intracellular signalling proxy).
# ============================================================================
const H_max_vox  = 100          # Hunger molecules/voxel at steady-state maximum
const k_H_deg    = 2.8f-4      
const k_H_prod   = Float32(H_max_vox * k_H_deg / MOLEC_PER_µM)  
const α_H        = 0.5f0       
const β_H        = 0.3f0        
# 0.8a  ∅ → Hunger   (starvation drives ACA-expression cascade)
r0_8a = Reaction(Symbol[], [:Hunger], k_H_prod, nothing)
# 0.8b  Hunger → ∅   (first-order turnover / saturation)
r0_8b = Reaction([:Hunger], Symbol[], k_H_deg,  nothing)

# ============================================================================
# Module 0.9 — INTRACELLULAR EXCITABLE ELEMENT: chemical FitzHugh–Nagumo
# (hunger_test.jl), ported in wholesale, plus a two-step ERK2→ACA output stage.
# [REPLACES the Hill-gated basal cyclase (old r9_0) and the Hill-gated FHN
#  (old r9_0b/c/d, species ACArefr) — 2026-07-30]
# ----------------------------------------------------------------------------
# WHY REPLACED.  r9_0 and r9_0b/c/d approximated an excitable element by
# writing Hill functions directly on ACAa_mem/cAMP_i/ACArefr.  hunger_test.jl
# (this repo, standalone Catalyst/JumpProcesses demo) shows the SAME dynamics —
# a stable rest state, an all-or-none threshold, noise-triggered spiking, and a
# refractory period — emerging from seven ELEMENTARY mass-action reactions among
# two new species with NO Hill function anywhere: every nonlinearity comes from
# reaction stoichiometry (order), exactly the "chemical Schlögl / FHN" construction.
# That file's own header comment (see hunger_test.jl) is the derivation of why
# this particular reaction set is bistable-with-noise-driven-spiking; it is
# reproduced verbatim below as species X (activator, FHN's v) and Y (recovery,
# FHN's w), gated on Hunger through the SAME competence(Hunger) Hill ramp that
# r9_0 used (kept — this is not FHN-specific, it is the developmental-clock gate).
#
# THE SEVEN REACTIONS (hunger_test.jl §SPECIES/§REACTIONS, k0 additionally gated
# by competence(Hunger) here instead of a raw (H/Hmax)^N ratio — same role):
#   (1) k0·comp(Hunger), ∅ → X      spontaneous rare activation (ignition noise)
#   (2) k1,  2X → 3X                autocatalytic +feedback (creates the threshold)
#   (3) k6,  3X → 2X                cubic self-limiting (bounds the spike, keeps
#                                    the SSA finite — see hunger_test.jl reaction 3)
#   (4) k2,  X → ∅                  baseline activator removal
#   (5) k3,  2X → 2X + Y            delayed inhibitor production (quadratic in X,
#                                    so baseline chatter in X does not raise it)
#   (6) k5,  X + Y → Y              inhibitor-mediated spike termination
#   (7) k4,  Y → ∅                  slow recovery decay (sets the refractory floor)
#
# THE OUTPUT STAGE (the delay the user asked for).  X does NOT drive ACA
# directly.  Naively wiring X → ACAa_mem would fire ACA into a fully-armed
# Module-10/11 brake: RegAp_cyto is still 100 % active, so 10.4 clears any
# cAMP_i the spike makes about as fast as 9.1/9.2 can produce it and the spike
# never reaches the world.  ERK2a_cyto is what disarms that brake — 10.2 routes
# active RegAp_cyto → the inactive RegAp_i_cyto pool only in the PRESENCE of
# ERK2a_cyto — so X must open that valve FIRST:
#   0.9e  X + ERK2_cyto → X + ERK2a_cyto      (X-driven, catalytic; fast: this is
#         the "first" step, independent of the receptor-driven r10_1)
#   0.9f  X + ACA_mem → X + ACAa_mem, rate × Hillup(ERK2a_cyto)   ("then": ACA
#         activation is additionally gated on ERK2a_cyto already being up, so it
#         only turns on once 0.9e has had time to build ERK2a_cyto and 10.2 has
#         had time to knock the RegA brake down — a real dynamical delay, not
#         just a rate-constant ordering, because ERK2a_cyto starts at 0 every
#         spike and 0.9f's Hill gate is shut until it rises)
# ⚠ FIRST-PASS PARAMETERS.  k_X_ERK2 / k_X_ACA / KE_ERKGATE below are a
# reasonable first cut (τ_ERK2a ≈ 1/(k_X_ERK2·X_spike + k_erk2Off) ~ sub-second
# once X is spiking, well inside the ~8 s spike hunger_test.jl documents; ACA
# gate half-activates at KE_ERKGATE ≈ 20 % of the ERK2 pool). NOT yet tuned
# against a population run (unlike the constants above them, which carry
# TUNING.md iteration numbers) — treat like the pre-Iteration-14/17 state of
# Module 0.9/0.9b and re-tune once this module is exercised end-to-end.
# All knobs ENV-tunable; DICTY_FHN_K0 = 0 restores X ≡ 0 (module inert).
const H_COMP_HALF = parse(Float32, get(ENV, "DICTY_COMP_HALF", "60"))   # Hunger/voxel at half-competence
const N_COMP      = 4                                                    # competence Hill steepness
const HC_N        = Float32(H_COMP_HALF)^N_COMP                          # baked into the @rate AST

# ── UNITS AND THE COMBINATORIAL / SPATIAL CORRECTION  [CALIBRATED 2026-07-30] ──
# The excitable element is a CELL-level decision device, so its molecule numbers
# are per CELL, not per voxel:  X rests at ~5, its excitation threshold is ~22
# and its spike top ~480 molecules/cell (= 7.7 nM in the 13-voxel, 104 µm³ cell).
# Reading hunger_test.jl's numbers as per-VOXEL counts instead puts the spike top
# at 5900/cell and the barrier at 4.1× rest across a 13× larger reaction volume;
# measured ignition is then ZERO spikes in 24 h against a target of one per
# 12 min (calib/CALIBRATION.md Iteration 1).  Hence M_CELL below.
const N_VOX_CELL = 13
const M_CELL     = Float64(MOLEC_PER_µM) * N_VOX_CELL   # 62628.8 molec / µM / cell
#
# Two engine properties have to be compensated explicitly, and BOTH matter here
# because X sits at only ~5 molecules/cell at rest:
#
#  (1) core_single/propensity.jl builds the RAW product over the educt list, so
#      "2X" gives rate·n² rather than rate·n(n−1) — a molecule is allowed to
#      react with itself.  At n ≈ 5 that is a 23 % error on the rate-limiting
#      step of an exponentially rate-sensitive escape process.
#  (2) The engine evaluates propensities PER VOXEL and X is spread over 13 of
#      them.  For molecules multinomially distributed over v voxels,
#          E[Σ nᵥ(nᵥ−1)]      = N(N−1)/v            ← exact
#          E[Σ nᵥ²]           = N²/v + N(v−1)/v     ← raw product: extra linear term
#      The leftover linear term acts as a spurious FIRST-order X source of
#      6·k1 = 0.48 s⁻¹, which cancels half of the k2 = 1.0 s⁻¹ decay and roughly
#      halves the excitation barrier.
#
# Both are removed by writing the same-species channels with ONE educt and
# carrying the falling factorial in the @rate factor: `[:X] → [:X,:X]` with
# frate (X−1) has the net stoichiometry of 2X→3X and the propensity
# rate·n(n−1).  With the exact form the per-voxel decomposition is EXACT, so a
# well-mixed 13-voxel cell is precisely one compartment with M = M_CELL and the
# calibration in calib/ transfers unchanged.  (n(n−1) and n(n−1)(n−2) self-zero
# at n = 0,1,2, so no guard is needed.)
#
# The constants below are named for hunger_test.jl's molecule-unit values
# (k0 … k6, per CELL); the engine rate each reaction is given is derived from
# them on the reaction line, so the ENV knobs stay comparable to that file.
# [TUNED iter 3/10/12/13]  X production at full competence, molecules per CELL per s.
# Calibrated directly against the ENGINE, which is the deliverable.
#
# Neither reduced model predicts the engine's firing rate well enough to set this:
# the well-mixed mirror is ~2× too slow (ignition nucleates in ONE voxel before
# diffusion levels it out, so a spatially extended excitable medium fires more
# readily), and the 13-voxel spatial mirror is still ~15 % too fast.  Both
# reproduce the PULSE exactly; only the rate needs the engine itself.
# Measured on the engine, 6 seeds × 9 h = 54 h, n = 229 intervals:
#   2.32 molec/s → 833 s ± 30 (95 % CI 775–892)   — 16 % slow
#   2.52 molec/s → see calib/CALIBRATION.md §4.1
# 2.52 molec/s ⇒ X_rest ≈ 2.7, X_threshold ≈ 24, X_top ≈ 454 molecules/cell.
# ── ⚠ [2026-09-16] 2.52 → 3.5.  THE 748 s CALIBRATION PREDATES THE 0.9r CLOCK. ─
# 2.52 was calibrated on the engine to a mean inter-pulse interval of 748 s
# (12.5 min) — the target — but that was measured BEFORE the 0.9r refractory
# clock existed (2026-08-09), and the clock adds a hard ~390 s dead time on top
# of the noise-escape wait.  CLAUDE.md records the consequence and names this
# constant as the fix: "⚠ The spontaneous ISI necessarily lengthens ... k_fhn_k0
# is the compensating knob and was deliberately NOT changed."  It is changed now.
#
# MEASURED on the engine (`calib/engine_check.jl`, the isolated Module 0/9/10
# cell at Hunger = 100 — the same harness and the same statistic as the 748 s
# calibration), 8 seeds x 6 h for the committed value, 4 x 6 h per sweep point:
#     k0 = 2.32  (historical, NO clock)   833.3 s   SE 29.8   CV 0.54   n = 229
#     k0 = 2.52  (historical, NO clock)   748.3 s   SE 23.4   CV 0.50   n = 251
#     k0 = 2.52  (WITH the clock)         930.4 s   SE 23.1   CV 0.33   n = 179
#     k0 = 3.5   (WITH the clock)         726.4 s   SE 10.5   CV 0.16   n = 117
# i.e. 3.5 restores the 12-minute interval the model is supposed to have, with
# the clock live.  95 % CI 706-747 s.
#
# THE PULSE IS UNTOUCHED, AND SO IS THE EXCITABLE CHARACTER.  Solving the
# cell-level X nullcline in the engine's own falling-factorial convention
# (0.9a k0, 0.9b (k1/2)X(X-1), 0.9c (k6/6)X(X-1)(X-2), 0.9d k2 X):
#     k0 = 2.52  rest 2.71  threshold 24.54  top 455.76   barrier 21.8
#     k0 = 3.5   rest 3.97  threshold 23.21  top 455.82   barrier 19.2
# The spike top moves by 0.01 %, the barrier by 12 %, and the fold (where rest
# and threshold merge and the cell would become a free-running OSCILLATOR rather
# than an excitable relay — the failure `problem.md` §1 documents) is far above:
# at k0 = 6.0 the barrier is still 9.95.  This is a change of ignition RATE, not
# of cell type.
#
# ⚠ THE COST IS REGULARITY, AND IT IS NOT A FREE PARAMETER.  CV falls 0.33 →
# 0.16 because the ~390 s block contributes no variance and now makes up more of
# the interval (predicted 0.5·336/726 = 0.23; measured 0.16).  Ford 2023 puts a
# real cell's CV nearer 0.3-0.5, so the model's pacemaker is now more clock-like
# than the literature.  That is inherent to "hard block + 12 min interval" and
# cannot be tuned out from here — it would need the block itself to be
# stochastic.  `DICTY_FHN_K0=2.52` restores the previous value.
#
# ── ✓ [2026-09-17] AND THE COLLECTIVE COUNTERPART, MEASURED ─────────────────
# In a population the cells do not wait out this interval -- the relay entrains
# them, so the period converges onto the refractory block (see the 0.9r block).
# Measured at 1000 cells / 400^2 / 3000 s, seeds 42/43/44: 467 +/- 16 s
# (kymograph) and 480 s (driver), i.e. 7.8-8.0 min against this cell's 12.1 min
# in isolation. Both numbers are wanted and they are not in conflict: 12 min is
# the PACEMAKER rate you asked Module 0 for, ~8 min is the WAVE rate the block
# imposes once relay is working.
const k_fhn_k0 = parse(Float32, get(ENV, "DICTY_FHN_K0", "3.5"))     # molec/s  [2.52 -> 3.5, 2026-09-16]
const k_fhn_k1 = parse(Float32, get(ENV, "DICTY_FHN_K1", "0.08"))    # molec⁻¹s⁻¹  autocatalytic gain (2X→3X)
const k_fhn_k6 = parse(Float32, get(ENV, "DICTY_FHN_K6", "0.0005"))  # molec⁻²s⁻¹  cubic self-limiting (3X→2X)
const k_fhn_k2 = parse(Float32, get(ENV, "DICTY_FHN_K2", "1.0"))     # s⁻¹       baseline X removal
const k_fhn_k3 = parse(Float32, get(ENV, "DICTY_FHN_K3", "6.25e-5")) # molec⁻¹s⁻¹  Y production (2X→2X+Y) [TUNED iter 2]
const k_fhn_k5 = parse(Float32, get(ENV, "DICTY_FHN_K5", "0.05"))    # molec⁻¹s⁻¹  Y-mediated X removal (X+Y→Y)
const k_fhn_k4 = parse(Float32, get(ENV, "DICTY_FHN_K4", "0.01"))    # s⁻¹       Y relaxation (τ=100 s; refractory) [TUNED iter S1]

# ── 0.9r — THE RELAY REFRACTORY CLOCK  [ADDED 2026-08-09, calib/CALIBRATION_M8.md]
# ----------------------------------------------------------------------------
# WHAT IT IS FOR.  Once Module 8 is live, a cell has TWO ways to fire — its own
# noise-driven ignition (0.9a) and a relayed cAMP wave (8.14) — and the
# requirement is a refractory period of ~8 minutes after EITHER.  Nothing in the
# network has that timescale.  Measured on the engine, every candidate is at
# least 3× too fast, and the shortfall is structural, not a matter of retuning:
#   Y (0.9y3)           τ = 100 s, peak ≈ 86/cell — and it CANNOT simply be
#                       slowed or enlarged, because Y is also the variable that
#                       TERMINATES the spike (0.9y2), so any change to it
#                       reshapes the calibrated pulse.
#   PKA* (10.D1/D2)     tracks cAMP_i with K_d = 0.1 µM and τ_off = 40 s, so it
#                       is gone ~2 min after cAMP_i clears.
#   RegA_p reset (10.B4) τ = 60 s.       ACA* (8.12+9.3) effective τ ≈ 80 s.
#   receptor Rp → R (1.7) t½ = 120 s.    FCD adapter (1.8b) τ = 50 s.
# So the refractory state is added explicitly rather than pretended into one of
# these.  Biologically it is the lumped DE-ADAPTATION / RESENSITIZATION state of
# the relay — the slow, post-pulse process that keeps a cell that has just
# signalled from signalling again, and whose 5–10 min duration is what sets the
# natural period of the aggregation-stage cAMP oscillation.  It is deliberately
# ONE lumped species: the model has no evidence to apportion it among receptor
# resensitization, cyclase re-priming and PKA-dependent transcriptional effects.
#
# WHY IT IS CHARGED BY X AND NOT BY cAMP_i.  X is the common final path of both
# firing routes — the whole point of 8.14 is that a relayed pulse is the same X
# spike as a spontaneous one — so charging the clock from X makes the block
# SOURCE-INDEPENDENT by construction.  Charging it from cAMP_i or ACA* would
# work too but ties the refractory length to the output stage's amplitude.
# Production is quadratic in X (same form as Y's 0.9y1, same falling-factorial
# treatment) so it is charged by SPIKES and essentially not at all by the
# resting X leak.
#
# WHAT IT GATES.  Both X sources, through the same Hill factor:
#   0.9a  ∅ → X          the spontaneous ignition seed
#   8.14  CRAC* → X      the relayed injection
# Gating the SEED (0.9a) rather than adding another X sink is what makes the
# block hard: with the seed shut, X drains to 0 through 0.9d and the
# autocatalytic 0.9b (which needs X ≥ 2) has nothing to amplify, so the cell
# cannot ignite by fluctuation either.
#
# ⚠ THIS CLOCK IS LIVE FOR EVERY DRIVER, including calib/engine_check.jl, and it
# gates 0.9a as well as 8.14 — it is a property of the cell, not of Module 8.
# (`DICTY_MODULES` does NOT gate it, or anything else: that variable is read
# nowhere in this file and never has been — it appears only in comments.  What
# makes engine_check.jl "minimal" is that it hands the Agent an `inital_state`
# and `diffusion_constants` covering only the Module 0/9/10 species, so every
# other module's species start and stay at 0.)
# So it CHANGES the calibrated single-cell statistics: measured 5 spikes in
# 5400 s × 2 seeds, i.e. a mean
# inter-pulse interval of ≈1080 s against calib/CALIBRATION.md's 748 s.  That is
# arithmetic, not a bug — a hard 480 s dead time is incompatible with a
# distribution whose measured minimum interval was 231 s, and the 8-minute floor
# was the stated requirement.  `DICTY_REFR_PROD=0` makes Refr identically zero,
# which sets both gates to 1 and restores the pre-2026-08-09 Module 0.9 exactly;
# use it to reproduce any number in calib/CALIBRATION.md.
#
# SIZING (calib/CALIBRATION_M8.md §6.1).  The block lifts when Refr decays below
# K_REFR, i.e. at t = τ_R·ln(R_peak/K_REFR) after the spike.  τ_R = 150 s and a
# peak/threshold ratio of ≈ e^3.2 ≈ 25 puts that at 480 s.  K_REFR is kept at a
# few hundred molecules — not a handful — so the boundary is set by the decay
# and not by shot noise on a near-empty species.
# ── SIZED FROM THE MEASURED RECOVERY CURVE, NOT FROM THE ANALYTIC ESTIMATE ──
# The block does not lift when Refr crosses K_REFR — it lifts when the gated
# relay drive crosses the excitation threshold of the Module 0.9 element, which
# is a different (and lower) point.  `calib/m8_probe.jl window`, six independent
# runs, spontaneous ignition off (DICTY_FHN_K0=0) so the only way to fire is the
# relay and the second pulse is not confounded by a chance spontaneous spike:
#
#   gap [s]   300      360      420      480      540      600
#   pulse 2   blocked  blocked  FIRED    FIRED    FIRED    FIRED
#   Refr@p2     442      290      191      119       73       42
#
# So the boundary is at Refr* ≈ 230/cell, i.e. 2.3× ABOVE K_REFR = 100: the
# relay breaks through while the gate is still only ~7 % open, because 100
# molec/s of drive against a ~24-molecule excitation threshold has that much
# headroom.  Measured R_peak ≈ 2900, τ_R = 150 s ⇒ boundary at
# 150·ln(2900/230) = 380 s post-spike, which is what the table shows and which is
# 6.3 minutes, not 8.
#
# THIS IS WHY THE ANALYTIC SIZING IN THE 0.9r BLOCK ABOVE IS NOT ENOUGH ON ITS
# OWN: it locates K_REFR, and the answer depends on the relay GAIN as well.  The
# calibration has to be closed on the measured curve.
#
# Iteration 2 (k_refrProd 3.2e-3, τ_R 180 s) measured R_peak = 3580 — NOT the
# 4600 the linear-in-k_refrProd estimate predicted — and moved the boundary to
# 446 s, i.e. 480 s was right ON it (blocked at gap 420, fired at 480).  The
# scaling is measurably sub-linear, so the peak has to be read off the run and
# not projected.  Committed: k_refrProd 5.0e-3, τ_R 200 s ⇒ R_peak ≈ 5600 and a
# boundary at 200·ln(5600/300) = 584 s post-spike.  584 and not 480 on purpose —
# the acceptance test is "no pulse before 480 s", the boundary is a stochastic
# one, and a criterion should not be evaluated at the point where it is marginal.
const k_refrProd = parse(Float32, get(ENV, "DICTY_REFR_PROD", "5.0e-3")) # molec⁻¹s⁻¹ (per cell), 2X → 2X + R
const k_refrOff  = parse(Float32, get(ENV, "DICTY_REFR_OFF",  "0.005"))  # s⁻¹  τ_R = 200 s
# ── [2026-08-20] K_REFR 100 -> 450: THE BLOCK WAS 11-12.5 MIN, NOT 8 ─────────
# MEASURED on HEAD with `calib/m8_probe.jl window`, DICTY_FHN_K0=0 (spontaneous
# ignition OFF — mandatory, see below), two seeds; gap at which pulse 2 fires:
#     K_REFR    seed 1    seed 2    block
#       100      720 s    >720 s    ~660-750 s  (11-12.5 min)   <- was committed
#       350      480 s     540 s     420-540 s  (7-9 min, seed-dependent)
#       450      480 s     480 s     420-480 s  (7-8 min)       <- committed
# The 584 s above was ANALYTIC; the engine puts the real boundary 100-170 s
# later still, so the block ran ~50-80 % longer than the 6-7 min intended. 450
# also makes the boundary SEED-REPRODUCIBLE (both seeds blocked@420/fire@480,
# where at 100 one seed fires at 720 and the other still does not) — the
# breakthrough now happens where Refr is larger, so its shot noise is smaller.
#
# WHY THIS IS THE POPULATION PERIOD. Two cells do NOT oscillate: one fires
# spontaneously, the other relays. Oscillation is collective — with many cells
# the probability that SOME cell's 0.9 fires the moment the block lifts goes to
# 1, so the period converges DOWN onto the block. Measured at 250 cells/54 um:
# period 756 s against a 660-750 s block, i.e. the residual wait is only
# ~0-90 s. Predicted period at 450: ~450-540 s.
#
# ── ⚠ [2026-09-17] THE TABLE ABOVE IS THE PRE-REPAIR RELAY GAIN. RE-MEASURED. ─
# This block says it itself: the boundary "depends on the relay GAIN as well".
# The 2026-09-16 Module 8 repair (k_cracDeph 10 -> 1.498, K_M8X_PIP3 13000 ->
# 5000; see calib/CALIBRATION_WAVE_2026-09-16.md) raised CRAC* from 11 % to 44 %
# of pool, so the 420-480 s in the K_REFR = 450 row no longer describes this
# network. Re-measured on the repaired relay, `calib/m8_probe.jl window` with
# DICTY_FHN_K0=0, THREE seeds:
#     gap [s]     300      360      420      480
#     seed 1    blocked  blocked   FIRED    FIRED
#     seed 2    blocked  blocked   FIRED    FIRED
#     seed 3    blocked   FIRED    FIRED    FIRED
# so: no cell fires within 300 s (3/3), 1 in 3 at 360 s, all 3 by 420 s --
# median boundary ~390 s, floor 5 min, full recovery 7 min. Pulse 1 fires 8/8
# and the recovery is monotone; on the pre-repair network pulse 1 fired in only
# 4 of 8 rows and pulse 2 was blocked at EVERY gap out to 720 s, which is a dead
# pathway rather than a refractory period.
# ⚠ `Refr@p2` at the breakthrough spans 614-956, i.e. the scatter is in how much
# clock each cell accumulated as well as in the crossing. ONE SEED IS NOT A
# BOUNDARY -- an n = 1 reading of "360-420 s" taken during that work was
# superseded by these replicates, and K_REFR was deliberately NOT changed, since
# every constant in this block also shifts the spontaneous ISI that k_fhn_k0 is
# set against.
#
# ── ✓ [2026-09-17] "Predicted period at 450: ~450-540 s" IS NOW MEASURED ─────
# 1000 cells / 400^2 / 3000 s on A100, seeds 42/43/44, two independent
# estimators: the driver's per-cell cAMP_ext gives a median 480 s in all three
# seeds, and the kymograph's whole-line mean gives 467 +/- 16 s (n = 15
# intervals, startup excluded). Inside the predicted 450-540 s band, above the
# 300-420 s block, and inside the 5-10 min literature band. 99-100 % of cells
# fire >= 2x.
# ⚠ CV is only 0.034 -- that IS this block's own mechanism ("the period
# converges DOWN onto the block") seen from the collective side, and it is the
# counterpart of the single cell's CV falling 0.33 -> 0.16 when k_fhn_k0 was
# raised to hit 12 minutes with the block live.
#
# ⚠ MEASURE THE BLOCK WITH DICTY_FHN_K0=0. With spontaneous ignition live a
# chance spike lands before the probe pulse and the row reads "pulse 1 no fire"
# (2/8 and 4/8 rows on two seeds) — that looks like an unreliable relay and is
# not one. With it off, pulse 1 fired 48/48 across six runs.
#
# ⚠ K_REFR gates 0.9a (spontaneous seed) as well as 8.14 (relay), so raising it
# shortens the single-cell ISI along with the block. Wanted here, but it means
# the ISI in calib/CALIBRATION.md moves and needs re-measuring.
const K_REFR     = parse(Float32, get(ENV, "DICTY_REFR_K",    "450.0"))  # Refr_cyto/CELL at half-block  [100 -> 450, 2026-08-20]
const N_REFR     = parse(Int,     get(ENV, "DICTY_REFR_N",    "4"))      # block Hill steepness
# Refr_cyto is a per-CELL variable like X and Y, but the engine evaluates the
# gate PER VOXEL, so the threshold has to be expressed per voxel — the species
# diffuses fast (D = 10) and is well mixed over the 13 voxels within ~2.5 s,
# which is what makes the per-voxel count a faithful proxy for the cell total.
const KR_VOX     = Float32(K_REFR / N_VOX_CELL)
const KR_N       = KR_VOX^N_REFR                                          # baked into the @rate AST
# 0.9r1  2X → 2X + Refr   (spike-charged; quadratic in X, exact falling factorial)
r0_9r1 = Reaction([:X], [:X, :Refr_cyto], Float32(N_VOX_CELL * k_refrProd / 2),
                  @rate(x -> x[:X] - 1f0))
# 0.9r2  Refr → ∅   (the recovery clock; τ_R sets how long the block lasts)
r0_9r2 = Reaction([:Refr_cyto], Symbol[], k_refrOff, nothing)

# 0.9a  ∅ → X   (spontaneous ignition, gated by competence(Hunger) AND by the
#   refractory clock — see the 0.9r block above)
#   zeroth order: per-voxel propensity = rate·M_VOX, so the cell total is
#   rate·M_CELL = k_fhn_k0 molecules/s.
r0_9a = Reaction(Symbol[], [:X], Float32(k_fhn_k0 / M_CELL),
                 @rate(x -> (x[:Hunger]^N_COMP / (HC_N + x[:Hunger]^N_COMP)) *
                            (KR_N / (KR_N + x[:Refr_cyto]^N_REFR))))
# 0.9b  2X → 3X   (autocatalytic +feedback — this is what creates the threshold)
#   written as [X] → [X,X] (net +1 X, identical stoichiometry) with frate (X−1);
#   engine rate = N_VOX_CELL · k1/2  so Σ_voxels = (k1/2)·N(N−1) per cell.
r0_9b = Reaction([:X], [:X, :X], Float32(N_VOX_CELL * k_fhn_k1 / 2),
                 @rate(x -> x[:X] - 1f0))
# 0.9c  3X → 2X   (cubic self-limiting; bounds the spike and keeps the SSA finite)
#   [X] → ∅ (net −1 X) with frate (X−1)(X−2); rate = N_VOX_CELL² · k6/6.
r0_9c = Reaction([:X], Symbol[], Float32(N_VOX_CELL^2 * k_fhn_k6 / 6),
                 @rate(x -> (x[:X] - 1f0) * (x[:X] - 2f0)))
# 0.9d  X → ∅   (baseline removal; first order, no correction needed)
r0_9d = Reaction([:X], Symbol[], k_fhn_k2, nothing)
# 0.9y1  2X → 2X + Y   (spike-triggered inhibitor production, quadratic in X)
r0_9y1 = Reaction([:X], [:X, :Y], Float32(N_VOX_CELL * k_fhn_k3 / 2),
                  @rate(x -> x[:X] - 1f0))
# 0.9y2  X + Y → Y   (spike termination; two distinct species ⇒ already exact)
r0_9y2 = Reaction([:X, :Y], [:Y], Float32(k_fhn_k5 * M_CELL), nothing)
# 0.9y3  Y → ∅   (slow recovery decay — sets the refractory floor)
r0_9y3 = Reaction([:Y], Symbol[], k_fhn_k4, nothing)

# ── Output stage: X → ERK2 (first) → ACA (then) ─────────────────────────────
# THE DELAY.  X does not drive ACA directly.  ERK2* has to disarm the RegA brake
# BEFORE ACA reaches full activity, or RegA hydrolyses cAMP_i as fast as ACA can
# make it and nothing reaches the world.  The delay is produced by the two arms
# having very different ERK2* thresholds, so it costs no intermediate species:
#   • RegA disarm (10.B3) needs only ERK2* > k_regAReact·M_CELL/k_regAInh ≈ 130
#     molecules/cell — crossed ~1 s after the X spike starts;
#   • the ACA gate (0.9f) half-opens at KE_ERKGATE·13 ≈ 5000 molecules/cell,
#     38× higher, and is crossed ~3–4 s later.
# Measured: at t = 4 s RegA_p is already 77 % down while ACA* is still 17 % of
# pool; ACA* peaks at 82 % at t ≈ 20 s.
#
# 0.9e is SECOND ORDER in X.  ERK2* is deliberately slow (τ_off = 50 s, Maeda
# 1996), so it integrates X; with a first-order drive the resting X leak left
# ERK2* at ~5 % of pool, which — because RegA integrates ERK2* in turn with
# τ = 60 s — left the brake 87 % OFF before any spike, capping the usable RegA
# knockdown at 11× (it is bounded by the ERK2* peak/rest ratio).  Squaring the
# X-dependence squares that contrast: resting ERK2* falls to 0.05 % of pool and
# the measured knockdown rises to 64×.  Biologically this is the dual (TEY)
# phosphorylation of the MAP kinase, the canonical ultrasensitive step; it also
# keeps Module 0.9's rule that every nonlinearity comes from reaction order.
const k_X_ERK2   = parse(Float32, get(ENV, "DICTY_FHN_XERK2",  "1823"))  # µM⁻²s⁻¹ (cell units) X²-driven ERK2 activation
const k_erk2Off  = parse(Float32, get(ENV, "DICTY_ERK2OFF",    "0.02"))  # s⁻¹     basal ERK2 dephosphorylation, τ = 50 s (Maeda 1996)
# [RECALIBRATED 2026-07-31 — calib/CALIBRATION.md Iteration 14] 32.5 → 65
# µM⁻¹s⁻¹ and the gate threshold 385 → 192/voxel (still 25 % of the — now
# larger — ERK2* pool at half-gate).  Both moved together to compress the ACA
# rise from ~23 s to ~10-12 s to 90 % of peak, which is what the task's "rises
# over 10 s" asks for; see the ACA-pool note at r0_9f below for why the pool
# itself also had to move.
const k_X_ACA    = parse(Float32, get(ENV, "DICTY_FHN_XACA",   "65.0"))  # µM⁻¹s⁻¹ X-driven ACA activation (post-gate)
const N_ERKGATE  = 4                                                       # ERK2a-gate Hill steepness
const KE_ERKGATE = parse(Float32, get(ENV, "DICTY_FHN_ERKGATE", "192"))  # ERK2a_cyto PER VOXEL for a half-open ACA gate (25 % of pool)
const KE_N       = KE_ERKGATE^N_ERKGATE                                    # baked into the @rate AST
# 0.9e  2X + ERK2_cyto → 2X + ERK2a_cyto   ("first": disarm the RegA brake via 10.B3)
#   two educts listed ⇒ engine divides by M_VOX once; the second X and the
#   falling factorial ride in the @rate term.  rate = k_X_ERK2 / M_VOX.
r0_9e = Reaction([:X, :ERK2_cyto], [:X, :ERK2a_cyto],
                 Float32(k_X_ERK2 / MOLEC_PER_µM),
                 @rate(x -> x[:X] - 1f0))
# 0.9e2 ERK2a_cyto → ERK2_cyto  (basal dephosphorylation; PKA adds the rest, 10.D3)
r0_9e2 = Reaction([:ERK2a_cyto], [:ERK2_cyto], k_erk2Off, nothing)
# 0.9f  X + ACA_mem → X + ACAa_mem   ("then": gated on ERK2a_cyto already being up)
# ACA_mem/ACAa_mem pool: 154 → 865/voxel (dicty_simulation.jl:ACA_v).  [ITER 14]
# The task specifies ~10⁷ cAMP_ext molecules exported per pulse at a per-ACA
# turnover of "the order of 20 cAMP/s" (k9_eff+k11_eff below).  ACA's copy
# number was never measured (same status as the RegA pool, Module10.md Part
# VIII); at the old 154/voxel (2002/cell) pool and k9+k11 = 30 s⁻¹ the pulse
# delivered 2.2-2.5×10⁶ molecules — the right order for Devreotes 1979's
# unstimulated relay pulse, but 4-5× short of the task's 10⁷ target. Since
# exported ≈ export_frac · k_cat · ∫ACA*(t) dt and export_frac/shape are set
# elsewhere, ∫ACA* dt has to scale with the target, and it does so linearly in
# the pool at a fixed peak occupancy fraction. Solved for numerically
# (calib/run.jl evoked): 865/voxel (0.18 µM, comparable in scale to the CAR1
# receptor-density literature value of ~10⁴-10⁵/cell) reproduces exactly
# 1.0×10⁷ molecules/cell exported with k9+k11 = 20 s⁻¹.
const ACA_v = parse(Int, get(ENV, "DICTY_ACA_V", "865"))   # ACA_mem+ACAa_mem per voxel
r0_9f = Reaction([:X, :ACA_mem], [:X, :ACAa_mem], k_X_ACA,
                 @rate(x -> x[:ERK2a_cyto]^N_ERKGATE /
                            (KE_N + x[:ERK2a_cyto]^N_ERKGATE)))

# ============================================================================
# Module 1 — cAMP Input & Receptor Binding   [RECALIBRATED 2026-07-31]
# ----------------------------------------------------------------------------
# Calibrated against Module1.pdf = Biswas, Devreotes & Iglesias (2021) PLoS
# Comput Biol 17(7):e1008803, "Three-dimensional stochastic simulation of
# chemoattractant-mediated excitability in cells" — Table 1 (GPCR module).  That
# paper is the closest published analogue of what this file does: a mesoscopic,
# lattice-SSA cAR1 model with explicit copy numbers.  Its reduced-order receptor
# model (Table 1 reactions 27/28: one free R and one occupied R:C state) has
# EXACTLY the topology of 1.1/1.2 here, and its full model supplies the affinity
# classes and the phosphorylated states.  See calib/CALIBRATION_M12.md for the
# full derivation, the timescale/rate-limiting analysis and the iteration log.
#
# Every constant is ENV-overridable so a calibration run can sweep one knob
# without editing this file.  DICTY_M12_LEGACY=1 restores the pre-2026-07-31
# values wholesale (used to reproduce the "what failed" baseline).
# ============================================================================
const M12_LEGACY = get(ENV, "DICTY_M12_LEGACY", "0") == "1"
_m12(env, new, old) = parse(Float32, get(ENV, env, string(M12_LEGACY ? old : new)))

# ── receptor topology: "reduced" (default) vs "biswas" (full multi-affinity) ─
# "reduced" is the pre-existing 2-state R/RC(+Rp/RpC) receptor — Table 1 rxns
# 27/28 of Modules.pdf, the paper's own "reduced-order model". Everything
# downstream of Module 1 in this file was calibrated against it.
# "biswas" switches in the FULL GPCR module (Modules.pdf Table 1, 10 species /
# 26 reactions): three affinity classes — high (H), low (L), slow (S) — each
# free/occupied, with H and L (not S — see Fig 2A) additionally carrying a
# phosphorylated/desensitized shadow pair. Set DICTY_M1_MODEL=biswas to use it;
# see the "Module 1b" block below (after r1_8b) for the reaction set.
const RECEPTOR_MODEL = get(ENV, "DICTY_M1_MODEL", "reduced")
RECEPTOR_MODEL in ("reduced", "biswas") ||
    error("DICTY_M1_MODEL must be \"reduced\" or \"biswas\", got $(repr(RECEPTOR_MODEL))")
const BISWAS_RECEPTOR = RECEPTOR_MODEL == "biswas"

# ── receptor copy number ────────────────────────────────────────────────────
# Aggregation-competent (starved 3–4 h) cells carry ~40,000 cAR1 per cell
# (Saxe et al. 1996; Pupillo et al. 1992) — the hard anchor of DICTY_CELL.md §0.
# Biswas Table 1 uses 70,000 ± 5,000; we keep the lower, directly-measured
# aggregation-stage number.  40,000 / 13 voxels = 3077 per voxel = 0.639 µM.
const R_v          = 3077                               # receptors / voxel = 0.639 µM
const R_TOTAL_CELL = R_v * N_VOX_CELL                   # 40,001 receptors / cell

# ── 1.1/1.2  unphosphorylated receptor: K_d = 60 nM ─────────────────────────
# Biswas Table 1 rxn 1/2 (high-affinity class H, taken from Van Haastert 1984's
# heterogeneous-binding analysis): k_H = 7.5 µM⁻¹s⁻¹, k_-H = 0.45 s⁻¹, so
# K_d = 60 nM — the measured "fast high-affinity" cAR1 state (~40 % of sites).
# WAS k_on = 100 µM⁻¹s⁻¹ / k_off = 1.0 s⁻¹ (K_d = 10 nM), which DICTY_CELL.md
# §1.x itself flags as "tighter than measured (~60 nM); raise k_off or lower
# k_on toward 60 nM if matching biochemistry".  The 6× loosening also moves the
# antenna out of saturation at the cAMP a neighbour actually delivers, which is
# what makes an occupancy GRADIENT survive at all (calib/CALIBRATION_M12.md §4).
# TIMESCALE: occupancy relaxes with τ = 1/(k_on·c + k_off); at c = K_d that is
# 1.1 s — the fastest step in the cascade, at quasi-equilibrium with everything
# downstream.  REGIME: binding-limited, not diffusion-limited — the Smoluchowski
# ceiling for cAMP (D = 350 µm²/s, target radius ~1 nm) is ≈2.6e3 µM⁻¹s⁻¹, some
# 350× above k_on, as expected for a GPCR with a buried orthosteric pocket.
const k_on   = _m12("DICTY_M1_KON",   7.5f0,  100.0f0)   # µM⁻¹s⁻¹
const k_off  = _m12("DICTY_M1_KOFF",  0.45f0,   1.0f0)   # s⁻¹   → K_d = 60 nM

# ── 1.3/1.4  phosphorylated receptor: 30× less sensitive ────────────────────
# The whole 30× is carried by the ON-rate (k_offp ≡ k_off), so the phosphorylated
# receptor stays fully competent to bind and release — it is SLOWER to engage,
# not locked out.  That is the form Biswas measured: their phospho-high state
# (Table 1 rxn 13/14) has k_PHC = 0.04 vs k_H = 7.5 µM⁻¹s⁻¹ — a 187× slower
# on-rate — and they report the phosphorylated states relaxing with t½ = 198 s
# against 10 s for the unphosphorylated ones.  The equilibrium-affinity loss
# alone is only 3–5× (Caterina et al. 1995; Biswas ref [43]), so the requested
# 30× sits inside the bracket [3–5× affinity, 187× kinetics] and is implemented
# on the axis the data actually constrain.  K_d,p = 30 × 60 nM = 1.8 µM.
# ⚠ 1/30 -> 1/8 [2026-08-30].  RpC K_d 1800 -> 480 nM, MATCHED TO THE MEASURED
# WAVE.  fire_cells.jl (commanded sender, 400 s) measures the peak extracellular
# cAMP at 983 molec/voxel = 204 nM at the SENDER's own membrane — the most any
# receptor in the population ever sees — and 39.5 nM at a receiver 60 µm away.
# With +30 % headroom the receptor must work to ~265 nM, so a K_d of 1800 nM sits
# far above anything the cell reaches and RpC never leaves the foot of its curve.
# 480 nM puts RpC at 30 % occupancy at the peak, i.e. on its responsive flank.
# Worst-case transfer over 2-265 nM improves 0.32 -> 0.40 and RpC's share of the
# drive at the peak falls 95 % -> 53 %, so R and RpC actually SPLIT the range.
# ⚠ calib/m8_probe.jl CLAMPS the world field and structurally cannot measure any
# of this — the numbers above are from fire_cells.jl.
const EPS_RP = _m12("DICTY_M1_EPSRP", Float32(1/8), Float32(3.3f-3/0.1f0))
const k_onp  = Float32(k_on * EPS_RP)                    # µM⁻¹s⁻¹
const k_offp = k_off                                     # s⁻¹  (affinity is all on-rate)

# ── 1.5/1.5b/1.6/1.7  the phosphorylation cycle ─────────────────────────────
# k_phos (occupied) : k_phos0 (unoccupied) = 40 : 1.  cAR1 phosphorylation is
# agonist-INDUCED (Vaughan & Devreotes 1988; Hereld et al. 1994) but has a real
# basal component in unstimulated cells, so 1.5b exists rather than being zero —
# that makes "occupied receptors phosphorylate more readily" a measurable ratio
# instead of a division by zero, and it gives the resting ~5 % Rp pool that
# basal phosphorylation implies (k_phos0/k_basal = 0.05).
# k_deph 0.05 → 0.005 s⁻¹ [CHANGED].  Two independent reasons:
#   (i)  DICTY_CELL.md §1.x criterion 2 wants ~80–90 % loss of the high-affinity
#        state under sustained saturating cAMP (Van Haastert 1984; Johnson 1991).
#        The desensitized fraction at steady state is k_phos/(k_phos+k_deph);
#        at 0.02/0.05 that was 29 % — the criterion FAILED, and §1.6's own text
#        ("slower than 1.5, so steady-state favours the desensitized state") was
#        inconsistent with its own number.  0.02/0.005 gives 80 %. ✓
#   (ii) t½ = 139 s lands in Biswas's measured 198 s phospho-state window and in
#        Caterina 1995's "recovery over minutes"; 0.05 s⁻¹ (t½ 14 s) did not.
# k_basal 0.01 → 0.00578 s⁻¹ [RECALIBRATED 2026-08-06, cAR1 two-state task].
# k_basal (1.7, Rp → R) is the physiologically relevant "dephosphorylation"
# rate: cAMP unbinding (k_off/k_offp) is fast (≲1-2 s timescale, essentially
# instantaneous on this clock), so within seconds of cAMP removal essentially
# every receptor sits in RC_mem→R_mem or RpC_mem→Rp_mem's UNBOUND side, and
# k_basal is what then rate-limits the total phosphorylated pool (Rp+RpC, the
# quantity an SDS-PAGE mobility-shift washout assay actually reports) decaying
# back down — k_deph (1.6, dephosphorylation of the still-BOUND RpC pool) is a
# distinct, faster-acting side channel that only matters while cAMP is still
# present (it is what sets the 80–90 % desensitized-fraction steady state
# under SUSTAINED cAMP, unaffected by this change; see (i) above).
# Task target: t½ ≈ 2 min (120 s) at room temperature (22 °C) — explicitly
# NOT the ~20 min figure reported at 4 °C, which the task identifies as a
# cold-arrest assay artifact rather than the physiological rate.
#   k_basal = ln2 / 120 s ≈ 0.00578 s⁻¹  (was 0.01, t½ = 69 s ≈ 1.15 min).
# Downstream check: criterion 4 ("re-arming far inside the pulse period")
# still holds comfortably — 120 s remains >6× shorter than the measured 748 s
# single-cell pulse period (CLAUDE.md), so this does not reopen that criterion.
# ⚠ [2026-08-25] 0.02 -> 1.0.  THE RECEPTOR PHOSPHO-CYCLE IS NOW THE DE-SATURATOR
# for Module 1, which is where the gradient is actually lost.  Read this with the
# EPS_GEF block in Module 2 — the two are ONE change and neither works alone.
#
# THE DEFECT.  After the Module 2 linearisation (M2_KRE_GAIN), the dominant
# compressor was no longer Module 2 but the RECEPTOR's own binding curve: at the
# 67–107 nM working dose against K_d = 60 nM, d ln RC/d ln c = 0.181, i.e. the
# receptor throws away 82 % of the contrast before the G-protein sees it.
#
# THE FIX, and why it was already in the network.  The phosphorylated receptor is
# 30× less sensitive (EPS_RP = 1/30, K_d 60 -> 1800 nM).  A LOW-affinity receptor
# is still on the steep part of its OWN binding curve at doses where the
# high-affinity one is flat — measured directly as a positive control
# (calib/analysis_2026-08-25/m1_phospho.jl §11, matched counts): moving K_d
# 15/60/240/960 nM gives VIS 0.55/1.00/1.56/1.87.  Occupancy pumps receptors into
# that state (1.5), so the pool migrates there exactly as the dose rises: the
# apparent K_d TRACKS the ambient and the stage never saturates.  Measured phospho
# fraction 10.2 % -> 34.9 % at rest and 39 % -> 98 % at the wave peak.
#
# ⚠ THIS ONLY WORKS IF RpC CARRIES THE DRIVE, AND THAT IS WHY IT LOOKED DEAD.
# README_M3_M4_STRUCTURAL_ANALYSIS.md Finding 6 swept k_phos with EPS_GEF pinned at
# 0.05 and concluded "turning desensitisation up does nothing or hurts".  That is
# reproduced here exactly — at EPS_GEF = 0.05 the local transfer DIPS before it
# recovers (0.068 -> 0.056 at k_phos = 0.1) — because with a 5 % leak the phospho
# receptor is a SINK, not a signal path: desensitisation removes drive and returns
# nothing. It is a pure loss until RpC's share of the drive clears ~15 %.  The
# k_phos × EPS_GEF table (§2) is monotone in neither variable alone and rises
# steeply along the diagonal.  PAIR THEM.
#
# WHY 1.0 s⁻¹ AND NOT LESS.  The pool has to re-equilibrate WITHIN the wave, and
# the wave is 40 s FWHM.  At 0.2 s⁻¹ the pump is too slow and VIS is 0.39× (worse
# than committed) at a low ambient; 0.5 gives 0.88×; 1.0 is the first setting that
# wins everywhere.  ⚠ THIS IS A REALISM COST AND IT IS LARGE: cAR1 phosphorylation
# is measured at t½ ~ 30–60 s (Vaughan & Devreotes 1988), so τ = 1 s is 30–60×
# FASTER THAN LITERATURE.  It is the same explicit trade as G_RATIO = 60 —
# function bought with realism, recorded so it is not mistaken for a measurement.
# k_deph/k_basal are deliberately NOT slowed (the user's "slower dephosphorylate"):
# k_basal/4 scores better on VIS (1.45/1.89/3.53 vs 1.20/1.37/2.78) but leaves the
# resting pool 68 % desensitised and stretches 90 % re-sensitisation to 601 s
# against a ~750 s ISI.  At the committed k_basal the cell resets in 352 s —
# FASTER than the committed network's own 460 s — and rests at 34.9 % phospho, so
# low-dose sensitivity (Module 8 relays from 2 nM) is preserved.
# `DICTY_M1_KPHOS=0.02` + `DICTY_M2_EPSGEF=0.05` + `DICTY_M2_GEFGAIN=1` reverts.
# ⚠ 1.0 -> 0.02 [2026-08-30].  AT k_phos = 1.0 THE UNPHOSPHORYLATED RECEPTOR IS
# NOT A SENSOR AT ALL.  `R_free`'s only source is Rp dephosphorylation at
# `k_basal`, so it is RESUPPLY-LIMITED and falls as 1/c, exactly cancelling the
# rise in binding:  RC = k_on*c*R_free/(k_off+k_phos) = k_basal*Rp/(k_off+k_phos),
# INDEPENDENT OF cAMP.  Measured over a 132x range (2 -> 265 nM): R_free
# 1059 -> 12 while RC moves only 11.0 -> 17.3 (1.5x), and RpC carries 95-98 % of
# the drive at every dose.  The cause is that k_phos/(k_phos+k_off) = 69 % of
# bindings drain R into the phospho branch.  At 0.02 the drain is 4 %, RC's range
# is 7.9x, and R does the sensing at the wave's leading edge while RpC takes over
# at the peak — the intended division of labour, which 1.0 makes impossible.
# 0.02 is also the LITERATURE value (cAR1 phosphorylation t½ 30-60 s, Vaughan &
# Devreotes 1988); 1.0 was 30-60x faster and flagged as such where it was set.
# ⚠ THIS IS A DESIGN/REALISM CHANGE, NOT A PERFORMANCE ONE.  Paired over 3 seeds
# on fire_cells.jl it moves gradient visibility by +2.76 ± 3.87 (Gbg) and
# +0.26 ± 0.39 (RasG-GTP) — NOT significant.  A first single-seed run showed a
# 4x Gbg gain and was RETRACTED; the other two seeds go the other way.  Adopt it
# because R becomes a sensor, not because the gradient improves.
const k_phos  = _m12("DICTY_M1_KPHOS",  0.02f0,   0.02f0)  # s⁻¹  RC → RpC, τ 50 s [0.02 -> 1.0 -> 0.02]
const k_phos0 = _m12("DICTY_M1_KPHOS0", 5.0f-4,   0.0f0)   # s⁻¹  R  → Rp   (basal, 40× slower)
const k_deph  = _m12("DICTY_M1_KDEPH",  0.005f0,  0.05f0)  # s⁻¹  RpC → RC, t½ 139 s
const k_basal = _m12("DICTY_M1_KBASAL", 0.00578f0, 0.01f0) # s⁻¹  Rp  → R,  t½ 120 s (2 min, 22°C)

# 1.1 R_mem + cAMP_ext → RC_mem
r1_1 = Reaction([:R_mem, :cAMP_ext], [:RC_mem],     k_on,   nothing)
# 1.2 RC_mem → R_mem + cAMP_ext
r1_2 = Reaction([:RC_mem], [:R_mem, :cAMP_ext],     k_off,  nothing)
# 1.3 Rp_mem + cAMP_ext → RpC_mem                       (30× weaker on-rate)
r1_3 = Reaction([:Rp_mem, :cAMP_ext], [:RpC_mem],   k_onp,  nothing)
# 1.4 RpC_mem → Rp_mem + cAMP_ext
r1_4 = Reaction([:RpC_mem], [:Rp_mem, :cAMP_ext],   k_offp, nothing)
# 1.5  RC_mem → RpC_mem   (agonist-induced desensitization)
r1_5 = Reaction([:RC_mem], [:RpC_mem],              k_phos,  nothing)
# 1.5b R_mem  → Rp_mem    (basal, ligand-independent phosphorylation — 40× slower)
r1_5b = Reaction([:R_mem], [:Rp_mem],               k_phos0, nothing)
# 1.6 RpC_mem → RC_mem   (resensitization, bound)
r1_6 = Reaction([:RpC_mem], [:RC_mem],              k_deph,  nothing)
# 1.7 Rp_mem → R_mem     (basal resensitization, free)
r1_7 = Reaction([:Rp_mem], [:R_mem],                k_basal, nothing)

const k_radOn  = 0.07f0    # [TUNE osc] 0.05→0.07 µM⁻¹s⁻¹: FCD adapter charges a bit faster so it
                           #            engages within one wave (drives the post-pulse reset, with ADAPT_AREF=40)
const k_radOff = 0.02f0    # s⁻¹      adapter relaxation, τ ≈ 50 s (slow vs rise, < 3-min refractory)
# 1.8a RC_mem + RcptAdapt_cyto → RC_mem + RcptAdapta_cyto  (RC catalytic; adapter pool conserved)
r1_8a = Reaction([:RC_mem, :RcptAdapt_cyto], [:RC_mem, :RcptAdapta_cyto], k_radOn,  nothing)
# 1.8b RcptAdapta_cyto → RcptAdapt_cyto                    (slow relaxation back to inactive)
r1_8b = Reaction([:RcptAdapta_cyto], [:RcptAdapt_cyto],    k_radOff,         nothing)

# ============================================================================
# Module 1b — Biswas, Devreotes & Iglesias (2021) full GPCR module  [LITERAL,
# REWRITTEN 2026-08-06]
# ----------------------------------------------------------------------------
# Only built when RECEPTOR_MODEL == "biswas" (DICTY_M1_MODEL=biswas). Eleven
# receptor/phosphate species replace R_mem/RC_mem/Rp_mem/RpC_mem: three
# affinity classes — high (H), low (L), slow (S) — each free (R_?_mem) /
# occupied (RC_?_mem), with H and L (not S — Fig 2A has no S phospho state)
# additionally carrying a phosphorylated/desensitized shadow pair
# (Rp_?_mem / RpC_?_mem), plus a free-phosphate pool Pi_mem (Table 1's "P").
#
# This block reproduces Modules.pdf ("Three-dimensional stochastic simulation
# of chemoattractant-mediated excitability in cells", Biswas/Devreotes/
# Iglesias, PLoS Comp Biol 2021) TABLE 1 LITERALLY — all 26 reactions (rows
# 1-26; rows 27/28 are the paper's own REDUCED model, i.e. this file's
# DEFAULT "reduced" mode above, and are not part of the full network) with
# the paper's own published rate constants, unmodified. This replaces a
# previous version of this block that retargeted the H/L on-rates to a
# task-specific Kd/koff spec and dropped the P-mediated rows 11-22 entirely
# (see git history) — that recalibration is gone; every number below is read
# directly off Table 1 (verified against the table image, not just the PDF's
# text extraction, since pdftotext garbles the superscript-P column headers).
#
# Receptor/G-protein copy numbers also switch to the paper's own literal
# values for this mode: 70,000 receptors/cell (Table 1 caption mean; NOT this
# file's measured 40,001/cell, which remains the "reduced" mode default —
# R_v/R_TOTAL_CELL above are untouched) and G-protein = 3× receptors
# (Table 2/3, same ratio the "reduced" mode already uses). The free-phosphate
# pool P has no number in the paper TEXT; 10,000 molecules/cell is Biswas's
# own value, taken from their published simulation code (param_setup /
# cAMP_dose_response.m in the Cell_Model_Iglesias_2021 repo's Matlab_Code:
# `n_P0_tot = 10000`, against `n_cAR1_free_tot = round(70000+5000*randn)` —
# same repo, same scale, the only available source). D_P = 0.27 µm²/s is the
# same code's `param.D_P = 10*2.7e-2*(1+0.40*randn)`, i.e. 10× the paper's
# GPCR diffusion constant.
#
# Downstream G-protein coupling (Module 2, below) is updated to match: Table
# 2/3's own kE weights ALL occupied classes (H:C, L:C, S:C, ᴾH:C, ᴾL:C)
# EQUALLY (RL = H:C+L:C+S:C+ᴾH:C+ᴾL:C is a single catalytic sum) — the ε=5%
# phospho-discount this file's "reduced"/pre-existing "biswas" modes apply to
# G-protein catalysis is a LOCAL invention, not in the paper, and is dropped
# for biswas mode's Gabg_mem coupling (see Module 2 below).
#
# The 1.8a/1.8b FCD-adapter reactions at the end of this block are NOT in the
# paper — they are this project's own bridge from receptor occupancy into
# Module 3 and are kept unchanged so downstream modules still receive a
# signal in biswas mode.
# ============================================================================
_m1b(env, val) = parse(Float32, get(ENV, env, string(val)))

# ── receptor copy number: paper's own literal value (biswas mode only) ──────
const R_TOTAL_CELL_BISWAS = _m1b("DICTY_M1_RTOT_BISWAS", 70000f0)  # Table 1 caption mean
const R_v_BISWAS = round(Int, R_TOTAL_CELL_BISWAS / N_VOX_CELL)    # 5385 / voxel

# ── H/L/S class cAMP binding — Table 1 rows 1-6, literal ─────────────────────
const k_H  = _m1b("DICTY_M1_KH",  7.5f0)     # µM⁻¹s⁻¹  (rxn 1)   H+cAMP -> H:C
const k_mH = _m1b("DICTY_M1_KMH", 0.45f0)    # s⁻¹      (rxn 2)   H:C -> H+cAMP,  Kd = 60 nM
const k_L  = _m1b("DICTY_M1_KL",  2.2f0)     # µM⁻¹s⁻¹  (rxn 3)   L+cAMP -> L:C
const k_mL = _m1b("DICTY_M1_KML", 1.0f0)     # s⁻¹      (rxn 4)   L:C -> L+cAMP,  Kd = 454 nM
const k_S  = _m1b("DICTY_M1_KS",  4.0f0)     # µM⁻¹s⁻¹  (rxn 5)   S+cAMP -> S:C
const k_mS = _m1b("DICTY_M1_KMS", 0.05f0)    # s⁻¹      (rxn 6)   S:C -> S+cAMP,  Kd = 12.5 nM

# ── H<->L class interconversion, free and occupied — Table 1 rows 7-10 ──────
const k_HL   = _m1b("DICTY_M1_KHL",   0.080f0)   # s⁻¹  (rxn 7)   H -> L, free
const k_mHL  = _m1b("DICTY_M1_KMHL",  0.0536f0)  # s⁻¹  (rxn 8)   L -> H, free
const k_HLC  = _m1b("DICTY_M1_KHLC",  0.080f0)   # s⁻¹  (rxn 9)   H:C -> L:C
const k_mHLC = _m1b("DICTY_M1_KMHLC", 0.0071f0)  # s⁻¹  (rxn 10)  L:C -> H:C

# ── phosphate-mediated phosphorylation/desensitization — Table 1 rows 11-22,
# literal. Pi_mem is the paper's free-phosphate pool "P" (see module note) ──
const k_HCP  = _m1b("DICTY_M1_KHCP",  5.3f-4)  # µM⁻¹s⁻¹ (rxn 11)  H:C+P   -> ᴾH:C
const k_mHCP = _m1b("DICTY_M1_KMHCP", 8.0f-4)  # s⁻¹     (rxn 12)  ᴾH:C    -> H:C+P
const k_PHC  = _m1b("DICTY_M1_KPHC",  4.0f-2)  # µM⁻¹s⁻¹ (rxn 13)  ᴾH+cAMP -> ᴾH:C
const k_mPHC = _m1b("DICTY_M1_KMPHC", 2.42f-3) # s⁻¹     (rxn 14)  ᴾH:C    -> ᴾH+cAMP
const k_PH   = _m1b("DICTY_M1_KPH",   5.3f-4)  # µM⁻¹s⁻¹ (rxn 15)  H+P     -> ᴾH
const k_mPH  = _m1b("DICTY_M1_KMPH",  8.0f-4)  # s⁻¹     (rxn 16)  ᴾH      -> H+P
const k_LCP  = _m1b("DICTY_M1_KLCP",  5.3f-4)  # µM⁻¹s⁻¹ (rxn 17)  L:C+P   -> ᴾL:C
const k_mLCP = _m1b("DICTY_M1_KMLCP", 8.0f-4)  # s⁻¹     (rxn 18)  ᴾL:C    -> L:C+P
const k_PLC  = _m1b("DICTY_M1_KPLC",  5.0f0)   # µM⁻¹s⁻¹ (rxn 19)  ᴾL+cAMP -> ᴾL:C
const k_mPLC = _m1b("DICTY_M1_KMPLC", 1.5f0)   # s⁻¹     (rxn 20)  ᴾL:C    -> ᴾL+cAMP
const k_PL   = _m1b("DICTY_M1_KPL",   5.3f-4)  # µM⁻¹s⁻¹ (rxn 21)  L+P     -> ᴾL
const k_mPL  = _m1b("DICTY_M1_KMPL",  8.0f-4)  # s⁻¹     (rxn 22)  ᴾL      -> L+P

# ── phospho-class interconversion — Table 1 rows 23-26, literal (unchanged
# from the pre-existing version of this block) ──────────────────────────────
const k_PHL   = _m1b("DICTY_M1_KPHL",   4.3f-4)  # s⁻¹      (rxn 23)  PH -> PL, free
const k_mPHL  = _m1b("DICTY_M1_KMPHL",  2.9f-4)  # s⁻¹      (rxn 24)  PL -> PH, free
const k_PHLC  = _m1b("DICTY_M1_KPHLC",  4.3f-4)  # s⁻¹      (rxn 25)  PH:C -> PL:C
const k_mPHLC = _m1b("DICTY_M1_KMPHLC", 3.8f-5)  # s⁻¹      (rxn 26)  PL:C -> PH:C

# ── resting population split (param_setup.m / cAMP_dose_response.m) ─────────
const F_H = _m1b("DICTY_M1_FH", 0.385f0)
const F_L = _m1b("DICTY_M1_FL", 0.575f0)
const F_S = _m1b("DICTY_M1_FS", 0.04f0)
const R_H_v = round(Int, R_v_BISWAS * F_H)
const R_L_v = round(Int, R_v_BISWAS * F_L)
const R_S_v = R_v_BISWAS - R_H_v - R_L_v   # remainder: the three always sum to R_v_BISWAS exactly

# ── free-phosphate pool P: 10,000 molecules/cell (Biswas's own simulation
# code — see module note above; not stated in the paper text) ───────────────
const Pi_TOTAL_CELL = _m1b("DICTY_M1_PTOT", 10000f0)
const Pi_v = round(Int, Pi_TOTAL_CELL / N_VOX_CELL)   # 769 / voxel

# 1/2   R_H_mem + cAMP_ext ⇌ RC_H_mem
r1b_1H = Reaction([:R_H_mem, :cAMP_ext], [:RC_H_mem], k_H,  nothing)
r1b_2H = Reaction([:RC_H_mem], [:R_H_mem, :cAMP_ext], k_mH, nothing)
# 3/4   R_L_mem + cAMP_ext ⇌ RC_L_mem
r1b_1L = Reaction([:R_L_mem, :cAMP_ext], [:RC_L_mem], k_L,  nothing)
r1b_2L = Reaction([:RC_L_mem], [:R_L_mem, :cAMP_ext], k_mL, nothing)
# 5/6   R_S_mem + cAMP_ext ⇌ RC_S_mem   (tight AND slow: never interconverts with H/L)
r1b_1S = Reaction([:R_S_mem, :cAMP_ext], [:RC_S_mem], k_S,  nothing)
r1b_2S = Reaction([:RC_S_mem], [:R_S_mem, :cAMP_ext], k_mS, nothing)
# 7-10  H <-> L class interconversion, free and occupied (S excluded)
r1b_HL   = Reaction([:R_H_mem],  [:R_L_mem],  k_HL,   nothing)
r1b_mHL  = Reaction([:R_L_mem],  [:R_H_mem],  k_mHL,  nothing)
r1b_HLC  = Reaction([:RC_H_mem], [:RC_L_mem], k_HLC,  nothing)
r1b_mHLC = Reaction([:RC_L_mem], [:RC_H_mem], k_mHLC, nothing)

# 11/12  H:C + Pi ⇌ RpC_H_mem   (occupied H phosphorylates WITHOUT losing cAMP)
r1b_HCP  = Reaction([:RC_H_mem, :Pi_mem], [:RpC_H_mem], k_HCP,  nothing)
r1b_mHCP = Reaction([:RpC_H_mem], [:RC_H_mem, :Pi_mem], k_mHCP, nothing)
# 13/14  RpC_H_mem ⇌ Rp_H_mem + cAMP_ext   (occupied phospho-H binds/releases cAMP)
r1b_PHC  = Reaction([:Rp_H_mem, :cAMP_ext], [:RpC_H_mem], k_PHC,  nothing)
r1b_mPHC = Reaction([:RpC_H_mem], [:Rp_H_mem, :cAMP_ext], k_mPHC, nothing)
# 15/16  R_H_mem + Pi ⇌ Rp_H_mem   (free H phosphorylates)
r1b_PH  = Reaction([:R_H_mem, :Pi_mem], [:Rp_H_mem], k_PH,  nothing)
r1b_mPH = Reaction([:Rp_H_mem], [:R_H_mem, :Pi_mem], k_mPH, nothing)
# 17/18  L:C + Pi ⇌ RpC_L_mem
r1b_LCP  = Reaction([:RC_L_mem, :Pi_mem], [:RpC_L_mem], k_LCP,  nothing)
r1b_mLCP = Reaction([:RpC_L_mem], [:RC_L_mem, :Pi_mem], k_mLCP, nothing)
# 19/20  RpC_L_mem ⇌ Rp_L_mem + cAMP_ext
r1b_PLC  = Reaction([:Rp_L_mem, :cAMP_ext], [:RpC_L_mem], k_PLC,  nothing)
r1b_mPLC = Reaction([:RpC_L_mem], [:Rp_L_mem, :cAMP_ext], k_mPLC, nothing)
# 21/22  R_L_mem + Pi ⇌ Rp_L_mem
r1b_PL  = Reaction([:R_L_mem, :Pi_mem], [:Rp_L_mem], k_PL,  nothing)
r1b_mPL = Reaction([:Rp_L_mem], [:R_L_mem, :Pi_mem], k_mPL, nothing)

# 23-26  phospho-class interconversion (Table 1, raw)
r1b_PHL   = Reaction([:Rp_H_mem],  [:Rp_L_mem],  k_PHL,   nothing)
r1b_mPHL  = Reaction([:Rp_L_mem],  [:Rp_H_mem],  k_mPHL,  nothing)
r1b_PHLC  = Reaction([:RpC_H_mem], [:RpC_L_mem], k_PHLC,  nothing)
r1b_mPHLC = Reaction([:RpC_L_mem], [:RpC_H_mem], k_mPHLC, nothing)

# ── 1.8a/1.8b FCD adapter — NOT part of Table 1; this project's own bridge
# into Module 3, kept unchanged from the pre-existing biswas block so the
# rest of the cascade still receives an occupancy signal (same catalytic role
# as the reduced model's) ────────────────────────────────────────────────────
r1b_8aH = Reaction([:RC_H_mem, :RcptAdapt_cyto], [:RC_H_mem, :RcptAdapta_cyto], k_radOn, nothing)
r1b_8aL = Reaction([:RC_L_mem, :RcptAdapt_cyto], [:RC_L_mem, :RcptAdapta_cyto], k_radOn, nothing)
r1b_8aS = Reaction([:RC_S_mem, :RcptAdapt_cyto], [:RC_S_mem, :RcptAdapta_cyto], k_radOn, nothing)

const MODULE_1_REACTIONS_BISWAS = [
    r1b_1H, r1b_2H, r1b_1L, r1b_2L, r1b_1S, r1b_2S,
    r1b_HL, r1b_mHL, r1b_HLC, r1b_mHLC,
    r1b_HCP, r1b_mHCP, r1b_PHC, r1b_mPHC, r1b_PH, r1b_mPH,
    r1b_LCP, r1b_mLCP, r1b_PLC, r1b_mPLC, r1b_PL, r1b_mPL,
    r1b_PHL, r1b_mPHL, r1b_PHLC, r1b_mPHLC,
    r1b_8aH, r1b_8aL, r1b_8aS, r1_8b,
]

# ============================================================================
# Module 2 — G-Protein Cycle   [RECALIBRATED 2026-07-31, re-tuned 2026-08-07]
# ----------------------------------------------------------------------------

const G_RATIO   = _m12("DICTY_M2_GRATIO", 60.0f0, 3.0f0)   # [6 -> 60, 2026-08-24]
const Gp_v      = round(Int, R_v * G_RATIO)              # 18462 Gαβγ / voxel = 3.83 µM
const Gp_v_BISWAS = round(Int, R_v_BISWAS * G_RATIO)     # 32310 Gαβγ / voxel (biswas mode, 70,000-receptor scale)
const RGS_v     = 385                                    # RGS / voxel (~5,000 / cell)

const t_half_diss = 3.2f0                                # s   Janetopoulos 2001

# ⚠ 1.463 -> 1.0 [2026-08-30].  1.463 VIOLATED THE ANCHOR IT IS DOCUMENTED AS
# BEING PINNED BY.  `k_gef` is defined below as GAIN·ln2/(t_half_diss·R_v/M), so
# GAIN = 1 is exactly what makes the Janetopoulos 3.2 s G-protein dissociation
# half-time hold; 1.463 gives 2.46 s — 30 % too fast.  (Measured on the
# deterministic receptor+G-protein cycle at saturating cAMP; this was already
# true BEFORE the 2026-08-30 Module 1 change, which is nearly anchor-neutral:
# 2.46 -> 2.28 s at unchanged gain.)
# ⚠ AND THE EXTRA GAIN WAS COSTING CONTRAST, NOT BUYING IT.  It over-drove the
# G-protein into saturation, and a saturated stage cannot report a difference:
#     d ln Gbg / d ln c   2 nM    20 nM    40 nM   204 nM
#       gain 1.463        0.673    0.158    0.092    0.048   (Gbg 29/72/78/86 % of pool)
#       gain 1.0          0.738    0.204    0.122    0.065   (Gbg 22/63/71/81 % of pool)
# Better at EVERY dose, for a 3 % loss of saturating Gbg counts.  The original
# reason for the raise — restoring molecules thrown away by desensitisation at
# k_phos = 1.0 (Gbg 82k -> 35k) — no longer exists at k_phos = 0.02.
const M2_GEF_GAIN = _m12("DICTY_M2_GEFGAIN", 1.0f0, 1.0f0)
const k_gef  = _m12("DICTY_M2_KGEF",
                    Float32(M2_GEF_GAIN * log(2f0) / (t_half_diss * (R_v / Float32(MOLEC_PER_µM)))),
                    2505.0f0)                            # µM⁻¹s⁻¹

# ⚠ KEPT AT 1.0 — RE-CHECKED 2026-08-30 AGAINST k_phos = 0.02, NOT INHERITED.
# The original argument for 1.0 (at k_phos = 1.0 a 5 % leak made the phospho
# receptor a SINK rather than a signal path) no longer applies, because far less
# receptor is phosphorylated now: EPS_GEF = 0.05 costs only 8 % of saturating
# Gbg here, against the gutting it caused before.  It is kept for a DIFFERENT
# and stronger reason — the R/RpC division of labour requires RpC to CARRY the
# peak, and at 0.05 it carries 7.2 % of the drive and transfer at the measured
# 204 nM peak collapses to 0.014 (vs 0.065 at 1.0).  A phospho receptor that
# neither binds (EPS_RP) nor couples (EPS_GEF) is a pure sink, and then nothing
# senses the top of the wave.
const EPS_GEF = _m12("DICTY_M2_EPSGEF", 1.0f0, Float32(0.25f0 / 2505.0f0))
const k_gefp  = Float32(k_gef * EPS_GEF)                 # µM⁻¹s⁻¹
# ── 2.0  basal, receptor-independent exchange ──────────────────────────────
# Biswas Table 3 rxn 1: k_E0 = 2.18×10⁻⁷ s⁻¹ at ratio 3.  Negligible in flux
# (~0.03 molecules/s per cell) but it is what sets the resting Gβγ floor, so the
# receiver has a defined pre-stimulus baseline rather than a hard zero.
const k_gef0 = _m12("DICTY_M2_KGEF0", 2.18f-7, 0.0f0)    # s⁻¹


# ⚠ 2026-09-12 CONTRAST RETUNE.  `k_hyd = ln2/30` (the Janetopoulos anchor)
# leaves M2 sitting deep in saturation at any dose near/above the 60 nM
# receptor K_d: u = drive/k_hyd is large there even though c itself is not
# "saturating cAMP", so d ln Gβγ/d ln c collapses to 0.10-0.12 at the 60-107 nM
# working-dose range (mean field, `calib/analysis_2026-08-25/m2_linear.jl`;
# committed engine measurement 0.224±0.027 at 60 nM, `calib/m12_probe.jl grad`).
# The lever the user specified is exactly right: k_hyd sets u (recombination
# clears the standing Gβγ pool faster, de-saturating it), and k_reassoc must
# be raised in step or the system falls back under the sqrt-law regime
# (r = k_reassoc·[Gβγ]/(M·k_hyd) collapses toward ≤1).  Both DECISIONS below
# are ENGINE measurements (`calib/m12_probe.jl grad`, n=8 replicate cells),
# not the mean-field proxy, because the two disagree sharply:
#   config       tr@60nM(mean-field)  tr@60nM(ENGINE)  SNR@60nM(ENGINE)
#   committed          0.122               0.224 ± 0.027     15.3
#   khyd×7/kreG×35     0.31-0.32           0.766 ± 0.087     31.7   <- WINS BOTH
#   khyd×100/kreG×1000 0.41                1.056 ± 0.111     13.4
# and at 5 nM the gap is even starker (committed 0.853±0.019 SNR 34.8;
# khyd×7 1.654±0.066 SNR 30.9 — essentially free; khyd×100 1.781±0.102 SNR 9.3,
# a real 3.7× SNR loss).  A naive `transfer·√(molecules)` mean-field visibility
# proxy predicts the OPPOSITE ordering (it says khyd×100 should be closest to
# a wash and khyd×7 should be worse at 5 nM) — it is wrong, and the lesson
# generalises: score M2 retuning on `calib/m12_probe.jl grad`'s own SNR/snapshot,
# never on a Poisson/mean-field stand-in (same class of error as
# `score-with-the-acceptance-estimator`).  Committed: `M2_KHYD_GAIN = 7`,
# `M2_KRE_GAIN` 100 → 3500 (=100×35, holding the old layer's meaning intact).
# ⚠ COSTS, measured, not yet paid down anywhere else in the file:
#  * Janetopoulos anchors move from PASS/PASS to PASS/FAIL: t½ diss 3.0→2.0 s
#    (still inside the 1.5-7 s window), t½ reassoc 33→7.5 s (target 15-60 s,
#    now a real ~4× break — far gentler than khyd×100's 100×-off 2.5 s).
#  * Gbg_v_rest (c=0) collapses 35 → 6/voxel, and every Module-3 constant
#    derived from it moves: `RasGTP_v_rest` 5356 → 2200/cell (measured on
#    `calib/m3_probe.jl pulse`: still structurally sound — peak/basal actually
#    RISES 4.02× → 9.39×, timing and reset criteria still pass — but M3's own
#    absolute counts are now ~41% of what its existing calibration assumed).
#    NOT re-validated beyond that one probe: `calib/CALIBRATION_M3*.md`'s
#    numbers were measured against the old Gbg_v_rest and should be re-run
#    before this is treated as acceptance-complete for the M1-M4 cascade.
#  * a0 (SSA cost) at the 60 nM working point rises ~3× (measured via the
#    steady-state identity a0_M2 = 3·drive(c)·G; rest-state cost is unchanged
#    since drive(0)≈0). Not the "free" paired-scaling trick from `G_RATIO`=60
#    — that held the ACTIVATED FRACTION fixed; this deliberately moves it, and
#    that is exactly what buys the contrast, so cost and contrast are coupled
#    here in a way they were not there.
# `DICTY_M2_KHYDGAIN=1` + `DICTY_M2_KREGAIN=100` restores the pre-2026-09-12
# state; `DICTY_M2_KHYD=0.05` (+ `DICTY_M2_KREGAIN=1`) goes all the way back
# to the pre-linearisation √-regime.
const M2_KHYD_GAIN = _m12("DICTY_M2_KHYDGAIN", 7.0f0, 1.0f0)
const k_hyd  = _m12("DICTY_M2_KHYD", Float32(M2_KHYD_GAIN * log(2f0)/30f0), 0.05f0)   # s⁻¹
const k_RGS  = _m12("DICTY_M2_KRGS", 5.0f0,  5.0f0)      # µM⁻¹s⁻¹  (r2_4 is OFF)

const t_half_reassoc = 30.0f0                            # s   Janetopoulos 2001
#

const M2_KRE_GAIN = _m12("DICTY_M2_KREGAIN", 3500.0f0, 100.0f0)
const k_reassoc = _m12("DICTY_M2_KREASSOC",
                       Float32(M2_KRE_GAIN / (t_half_reassoc * (Gp_v / 2) / Float32(MOLEC_PER_µM))),
                       10.0f0)                           # µM⁻¹s⁻¹

# 2.0 Gabg_mem → Ga2_GTP_mem + Gbg_cyto                     (basal exchange)
r2_0 = Reaction([:Gabg_mem], [:Ga2GTP_mem, :Gbg_cyto],  k_gef0,   nothing)
# 2.1 Gabg_mem + RC_mem → Ga2_GTP_mem + Gbg_cyto + RC_mem
r2_1 = Reaction([:Gabg_mem, :RC_mem],  [:Ga2GTP_mem, :Gbg_cyto, :RC_mem],  k_gef,  nothing)
# 2.2 Gabg_mem + RpC_mem → Ga2_GTP_mem + Gbg_cyto + RpC_mem (ε = 5 % of 2.1)
r2_2 = Reaction([:Gabg_mem, :RpC_mem], [:Ga2GTP_mem, :Gbg_cyto, :RpC_mem], k_gefp, nothing)
# 2.3 Ga2_GTP_mem → Ga2_GDP_mem   (intrinsic hydrolysis)
r2_3 = Reaction([:Ga2GTP_mem], [:Ga2GDP_mem],          k_hyd,  nothing)
# 2.4 Ga2_GTP_mem + RGS_mem → Ga2_GDP_mem + RGS_mem         (RGS-accelerated GTPase)
#r2_4 = Reaction([:Ga2GTP_mem, :RGS_mem], [:Ga2GDP_mem, :RGS_mem],  k_RGS,  nothing)
# 2.5 Ga2_GDP_mem + Gbg_cyto → Gabg_mem  (reassociation)
r2_5 = Reaction([:Ga2GDP_mem, :Gbg_cyto], [:Gabg_mem], k_reassoc, nothing)

const k_gef_biswas = _m1b("DICTY_M2_KE_BISWAS", 0.81f-2)   # µM⁻¹s⁻¹  Table 2/3, ratio 3, literal

r2_1H_lit  = Reaction([:Gabg_mem, :RC_H_mem],  [:Ga2GTP_mem, :Gbg_cyto, :RC_H_mem],  k_gef_biswas, nothing)
r2_1L_lit  = Reaction([:Gabg_mem, :RC_L_mem],  [:Ga2GTP_mem, :Gbg_cyto, :RC_L_mem],  k_gef_biswas, nothing)
r2_1S_lit  = Reaction([:Gabg_mem, :RC_S_mem],  [:Ga2GTP_mem, :Gbg_cyto, :RC_S_mem],  k_gef_biswas, nothing)
r2_1PH_lit = Reaction([:Gabg_mem, :RpC_H_mem], [:Ga2GTP_mem, :Gbg_cyto, :RpC_H_mem], k_gef_biswas, nothing)
r2_1PL_lit = Reaction([:Gabg_mem, :RpC_L_mem], [:Ga2GTP_mem, :Gbg_cyto, :RpC_L_mem], k_gef_biswas, nothing)

const MODULE_2_REACTIONS_BISWAS = [r2_0, r2_1H_lit, r2_1L_lit, r2_1S_lit, r2_1PH_lit, r2_1PL_lit, r2_3, r2_5]

const D_GA2  = _m12("DICTY_M2_DGA2",  0.2f0, 0.2f0)   # free Gα2, C14 myristoyl
const D_GABG = _m12("DICTY_M2_DGABG", 0.2f0, 0.2f0)   # heterotrimer, Biswas T3
# Gβγ 0.1 → 0.05 µm²/s [2026-08-07].  Still inside Diffusion_constants.md §4.6's
# 0.05–0.1 peripheral-membrane bracket — this is the slow end of the SAME
# bracket, not a new claim — and monotonically better at every operating point
# tested (`tune2`: VIS 5.59 / 5.17 / 4.56 at D = 0.05 / 0.10 / 0.20).
# ⚠⚠ 0.05 → 0.2 [2026-09-29], DEFAULT: EQUAL MOBILITY OF ALL THREE G-PROTEIN
# FORMS. The 0.05/0.2 mismatch strands free Gβγ behind the Gα it came from and
# makes Gβγ a spurious post-wave amplifier (centroid gap 0.74 vs 0.004 with equal
# D; post-wave polarity 0.52 vs 0.02 — calib/PERIODIC_MEMORY_2026-09-26.md §1).
# Elzie et al. 2009 measured SIMILAR FRAP recovery for Gα2 and Gβγ (recovery by
# cytosolic exchange, no D quoted), so the VALUE 0.2 is an [ESTIMATE] and the
# EQUALITY is the sourced part. Population effect alone is small (factorial §8:
# D2 ≈ R on displacement); it is required for the combined gain.
# `DICTY_M2_DGBG=0.05` restores the previous network.
const D_GBG  = _m12("DICTY_M2_DGBG",  0.2f0, 10.0f0) # geranylgeranyl CAAX dimer
const M12_DIFFUSION_REDUCED = Dict{Symbol, Float32}(
    :R_mem   => 0.024f0, :RC_mem  => 0.024f0,       # integral TM, Takebayashi 2023
    :Rp_mem  => 0.024f0, :RpC_mem => 0.024f0,
    :Gabg_mem   => D_GABG,
    :Ga2GTP_mem => D_GA2, :Ga2GDP_mem => D_GA2,
    :Gbg_cyto   => D_GBG,
    #:RGS_mem => 0.024f0,                            # membrane-anchored GAP
    :RcptAdapt_cyto => 10.0f0, :RcptAdapta_cyto => 10.0f0,   # genuinely cytosolic
)
# biswas mode: same numbers, spread across the 10 affinity-class species — all
# GPCR states diffuse at the same 0.024 µm²/s (Table 1 uses one D for all of them).
# Pi_mem (Table 1's free-phosphate pool P) has no diffusion constant in the
# paper text; 0.27 µm²/s is Biswas's own simulation code (param.D_P =
# 10*2.7e-2 = 10× the GPCR constant — see the Module 1b block note above).
const M12_DIFFUSION_BISWAS = Dict{Symbol, Float32}(
    :R_H_mem => 0.024f0, :R_L_mem => 0.024f0, :R_S_mem => 0.024f0,
    :RC_H_mem => 0.024f0, :RC_L_mem => 0.024f0, :RC_S_mem => 0.024f0,
    :Rp_H_mem => 0.024f0, :Rp_L_mem => 0.024f0,
    :RpC_H_mem => 0.024f0, :RpC_L_mem => 0.024f0,
    :Pi_mem => 0.27f0,
    :Gabg_mem   => D_GABG,
    :Ga2GTP_mem => D_GA2, :Ga2GDP_mem => D_GA2,
    :Gbg_cyto   => D_GBG,
   #:RGS_mem => 0.024f0,
    :RcptAdapt_cyto => 10.0f0, :RcptAdapta_cyto => 10.0f0,
)
const M12_DIFFUSION = BISWAS_RECEPTOR ? M12_DIFFUSION_BISWAS : M12_DIFFUSION_REDUCED

const _rp_frac  = k_phos0 / (k_phos0 + k_basal)
const Rp_v_rest = round(Int, R_v * _rp_frac)
const R_v_rest  = R_v - Rp_v_rest

"Zero-cAMP (Ga2GTP, Ga2GDP, Gβγ) floor set by the 2.0/2.3/2.5 balance, for a pool of `gp`/voxel."
function _m2_rest(gp)
    M   = Float64(MOLEC_PER_µM)
    p   = Float64(k_gef0); kh = Float64(k_hyd); kr = Float64(k_reassoc)
    gpf = Float64(gp)
    resid(B) = (G = gpf - B; p*G/kh + p*G*M/(kr*B) - B)
    lo, hi = 1e-12, gpf*(1 - 1e-12)
    for _ in 1:200
        m = (lo + hi)/2
        resid(m) > 0 ? (lo = m) : (hi = m)
    end
    B = (lo + hi)/2; T = p*(gpf - B)/kh
    # ⚠ ROUND ONTO THE MANIFOLD, NOT INDEPENDENTLY.  `Gβγ ≡ Ga2GTP + Ga2GDP` is
    # an exact conservation law of this topology, so three separately-rounded
    # integers put the IC OFF it (1.74/32.52/34.26 → 2/33/34, and 2+33 ≠ 34).
    # The engine cannot restore a conservation law it is started off, so Gβγ is
    # derived from the two rounded Gα states rather than rounded itself.
    Ti = round(Int, T); Pi = round(Int, B - T)
    (Ti, Pi, Ti + Pi)
end
_gbg_rest(gp) = _m2_rest(gp)[3]
const Ga2GTP_v_rest, Ga2GDP_v_rest = _m2_rest(Gp_v)[1], _m2_rest(Gp_v)[2]
const Gbg_v_rest = _gbg_rest(Gp_v)
const M12_INITIAL_REDUCED = Dict{Symbol, Int}(
    :R_mem => R_v_rest, :RC_mem => 0, :Rp_mem => Rp_v_rest, :RpC_mem => 0,
    :RcptAdapt_cyto => 154, :RcptAdapta_cyto => 0,
    # [2026-08-25] Ga2GTP is no longer pinned at 0: with k_hyd = ln2/30 the
    # GTP-loaded state is the one that HOLDS the free pool at rest, not Ga2GDP.
    :Gabg_mem   => Gp_v - Ga2GTP_v_rest - Ga2GDP_v_rest,
    :Ga2GTP_mem => Ga2GTP_v_rest, :Ga2GDP_mem => Ga2GDP_v_rest, :Gbg_cyto => Gbg_v_rest,
    #:RGS_mem => RGS_v,
)
# biswas mode: same discipline, applied to the ten-species receptor.  The
# resting split is NOT a free choice here either — at c = 0 the free receptor
# runs a four-state cycle (H ⇌ L via 7/8, H ⇌ ᴾH via 15/16, L ⇌ ᴾL via 21/22,
# ᴾH ⇌ ᴾL via 23/24) that consumes free phosphate, and `Pi_mem` is a CONSERVED
# pool, so the phospho fraction and the free phosphate left over determine each
# other.  Starting every receptor free with `Pi_mem` at its full 769/voxel is
# again the "everything off" corner: it drains ~40 % of the phosphate pool over
# the first several hundred seconds while ᴾH/ᴾL fill up.
#
# `R_S_mem` is the exception and genuinely does stay put: rows 5/6 need cAMP and
# S takes no part in the H⇌L or phospho interconversions (Fig 2A has no S
# phospho state), so with c = 0 it has no exit at all.
#
# Solved numerically rather than in closed form: the cycle is not in detailed
# balance (Table 1's four rate pairs do not satisfy the Kolmogorov condition
# around the H→L→ᴾL→ᴾH→H loop), so the stationary state is a genuine
# steady-state flux solution.  The mass-action ODE is relaxed to that fixed
# point below — Pi is conserved by construction, so the iteration cannot drift
# off the pool.
function _biswas_rest(N, pi_tot; dt = 1.0, nstep = 400_000)
    H, L, PH, PL = Float64(N), 0.0, 0.0, 0.0          # start anywhere on the manifold
    for _ in 1:nstep
        p  = max(pi_tot - PH - PL, 0.0) / Float64(MOLEC_PER_µM)   # free phosphate [µM]
        dH  = -k_HL * H + k_mHL * L - k_PH * p * H + k_mPH * PH
        dL  =  k_HL * H - k_mHL * L - k_PL * p * L + k_mPL * PL
        dPH =  k_PH * p * H - k_mPH * PH - k_PHL * PH + k_mPHL * PL
        dPL =  k_PL * p * L - k_mPL * PL + k_PHL * PH - k_mPHL * PL
        H += dt*dH; L += dt*dL; PH += dt*dPH; PL += dt*dPL
    end
    (round(Int, H), round(Int, L), round(Int, PH), round(Int, PL))
end
const _RH_r, _RL_r, _RPH_r, _RPL_r = _biswas_rest(R_H_v + R_L_v, Pi_v)
const Gbg_v_rest_BISWAS = _gbg_rest(Gp_v_BISWAS)
const M12_INITIAL_BISWAS = Dict{Symbol, Int}(
    :R_H_mem => _RH_r, :R_L_mem => _RL_r, :R_S_mem => R_S_v,
    :RC_H_mem => 0, :RC_L_mem => 0, :RC_S_mem => 0,
    :Rp_H_mem => _RPH_r, :Rp_L_mem => _RPL_r, :RpC_H_mem => 0, :RpC_L_mem => 0,
    :Pi_mem => Pi_v - _RPH_r - _RPL_r,
    :RcptAdapt_cyto => 154, :RcptAdapta_cyto => 0,
    :Gabg_mem   => Gp_v_BISWAS - Gbg_v_rest_BISWAS,
    :Ga2GTP_mem => 0, :Ga2GDP_mem => Gbg_v_rest_BISWAS, :Gbg_cyto => Gbg_v_rest_BISWAS,
   # :RGS_mem => RGS_v,
)
const M12_INITIAL = BISWAS_RECEPTOR ? M12_INITIAL_BISWAS : M12_INITIAL_REDUCED

# ============================================================================
# Module 3 — RasG Activation + Adaptation (incoherent feed-forward / LEGI)
# ----------------------------------------------------------------------------

const M3_LEGACY = get(ENV, "DICTY_M3_LEGACY", "0") == "1"
_m3(env, new, old) = parse(Float32, get(ENV, env, string(M3_LEGACY ? old : new)))

const RasGv  = parse(Int, get(ENV, "DICTY_M3_RASGV",  "38500")) # RasG per voxel  (5.0e5/cell) [385 -> 38500, 2026-08-08]
const GEFR_v = parse(Int, get(ENV, "DICTY_M3_GEFRV", "138600")) # RasGEFR_cyto per voxel (1.8e6/cell) [154 -> 4620 -> 46200 -> 138600, 2026-08-24]

const GAP_v  = parse(Int, get(ENV, "DICTY_M3_GAPV",   "46200")) # RasGAP_cyto per voxel   (6.0e5/cell) [154 -> 4620 -> 46200, 2026-08-24]

# ⚠⚠ 154 → 1540 [2026-09-15]. A PURE SHOT-NOISE RAISE, and the single largest
# lever found in the whole 2026-09 campaign. Committed together with
# BRAKE_HALF 2.17 → 1.5 and BRAKE_STRENGTH 0.3 → 0.6 (see their block).
#
# WHY. At 154/voxel the brake RESTS AT ~22 MOLECULES and regulates a 38500/voxel
# Ras pool — two to three orders below everything it competes with (RasG 38500,
# GEFR_v/GAP_v 46200, PTEN 10000). Its own 1/sqrt(22) = 21 % fluctuation is
# injected straight into the hydrolysis rate by 3.3e.
#
# ⚠ THE MEAN BRAKING ACTION IS EXACTLY INVARIANT, so this is noise-only:
# k_brakeHydro is defined proportional to 1/RasBrake_v while the braking flux
# carries [RasBrakea] = f·RasBrake_v, so the pool cancels. VERIFIED at a FIXED
# BRAKE_STRENGTH: k_brakeHydro·RasBrake_v = 30639.9 at 154, 1540, 15400 and
# 46200 alike (measured at BRAKE_STRENGTH = 0.3). The product is proportional to
# BRAKE_STRENGTH by construction, so at the committed 0.6 it reads 61279.9 —
# that is the STRENGTH change, not a break in the pool invariance.
#
# ⚠ THE EFFECT SATURATES AT ~1540. p300_tuned (46200, seed 701) measured
# over-dispersion 42.2 / alignment 0.670 against pool10's (1540, seed 700)
# 43.1 / 0.670 — identical to two digits across a 30× pool AND two seeds. Bigger
# pools only cost runtime: the 46200 arm ran ~2× slower. ⚠ That also retired my
# own "~2.7 % cost" estimate, which used a stale a0 ≈ 1e5 from CLAUDE.md and was
# wrong by ~40× — a cost figure belongs to a configuration, not to a code base.
#
# MEASURED, at MATCHED input contrast (the unmatched comparison is confounded:
# the brake halves the accepted frames' input contrast and doubles acceptance,
# which alone moves over-dispersion — see analysis/matched_drive.jl), bin
# 0.04-0.06, against the previous default:
#   PIP3 over-dispersion   234.4 -> 43.1      RasG_GTP over-dispersion 29.5 -> 11.0
#   PIP3 alignment         0.591 -> 0.670     derived reliability      +62 %
# and the over-dispersion STOPS GROWING WITH DRIVE (previous default 173->379
# across the four bins, this 40->46 flat), with alignment holding at strong
# drive instead of collapsing (0.374 -> 0.665) — the signature of a loop that is
# damped rather than running away.
#
# VERIFIED IN A POPULATION, seed-PAIRED over six seeds (700-705), 4000 cells,
# 400 k steps, only these three constants differing:
#   straightness  0.2401 (0.2369-0.2440) -> 0.2925 (0.2866-0.3052)
#                 paired +0.0535 ± 0.0050, t = 26.0, 6/6, ranges DO NOT OVERLAP
#   clustered     0.5490 -> 0.6781, paired +0.131 ± 0.018, t = 18.0, 6/6
#   Rg ratio      0.99935 -> 0.99854, 6/6 lower, t = 3.6
#   mean-NN       0.6398 -> 0.6411  (no change)
#   largest_frac  0.0523 -> 0.0550  (no change; its control spread is 0.033-0.093)
# ⚠ Two n=1 claims made during the sweep were WRONG and are retracted by this:
# "mean-NN gets worse" (0.635 -> 0.664 on one seed; paired it is -0.0003) and
# "largest_frac +89 %" (paired it is zero). Both were seed noise.
const RasBrake_v = parse(Int, get(ENV, "DICTY_M3_BRAKEV", "1540")) # RasBrake per voxel (~20,000/cell) — see 3.3c/d/e below

const k_gefRon  = _m3("DICTY_M3_KGEFRON",  0.0038f0, 1.0f0)  # µM⁻¹s⁻¹  [3.0 -> 0.038 -> 0.0038, 2026-08-24]
const k_gefRoff = _m3("DICTY_M3_KGEFROFF", 0.5f0, 0.1f0)     # s⁻¹  τ_E = 2 s  [2.0 -> 0.5]

const k_rasGon  = _m3("DICTY_M3_KRASGON",  0.333325f0, 15.0f0)  # µM⁻¹s⁻¹  [200 -> 0.275 -> 4.0 -> 1.3333, 2026-08-24]

# ── τ_I = 30 s → 10 s [2026-10-02], ALL THREE ARMING RATES SCALED WITH IT ────────
# Takeda et al. 2012 (Supp. Table S1, p. 22) fitted k−GAP = 0.1 s⁻¹ (τ = 10 s) and
# k−GEF = 0.4 s⁻¹; "k-GAP determined the time scale of the return to basal amount"
# (Supp. p. 6).  The model had τ_I = 30 s.  k_gapOff is set to Takeda's value and the
# arming rates of 3.4, 3.4e and 3.4b are multiplied by the SAME factor 0.1/0.033, so
# every activated GAP fraction (rest, adapted plateau, _fI_rest, _fI_rest_G and hence
# k_nf1Cat) is unchanged and only the time constant moves — the pairing Takeda used
# himself ("we fixed kGAP=0.1k-GAP", Supp. p. 5).  ⚠ NOT RE-CALIBRATED: the faster
# inhibitor shortens and lowers the Ras transient; re-check the step return (<35 s,
# Takeda p. 2) and the downstream modules.  `DICTY_M3_KGAPOFF=0.033
# DICTY_M3_KGAPONG=0.000253 DICTY_M3_KGAPBASAL=5e-4 DICTY_M3_KGAPON=0.08` restores τ_I = 30 s.
const _M3_TAUI_PAIR = Float32(0.1 / 0.033)                    # 3.0303, k_off(new)/k_off(old)
const k_gapOn  = _m3("DICTY_M3_KGAPON",  0.08f0 * _M3_TAUI_PAIR, 0.05f0)  # µM⁻¹s⁻¹  [0.08 -> 0.2424, 2026-10-02]
const k_gapOff = _m3("DICTY_M3_KGAPOFF", 0.1f0, 0.005f0)     # s⁻¹  τ_I = 10 s (Takeda 2012)  [0.0015 -> 0.033 -> 0.1]

const M3_INHIB = get(ENV, "DICTY_M3_INHIB", "gbg")
M3_INHIB in ("rc", "gbg", "both") ||
    error("DICTY_M3_INHIB must be rc|gbg|both, got $(repr(M3_INHIB))")

const k_gapOnG = _m3("DICTY_M3_KGAPONG", 0.000253f0 * _M3_TAUI_PAIR, 0.0f0) # µM⁻¹s⁻¹  (3.4e) [0.00253 -> 0.000253 -> 0.000767 (τ_I pairing), 2026-10-02]

const k_rasGoff = _m3("DICTY_M3_KRASGOFF", 5.5f0, 150.0f0)  # µM⁻¹s⁻¹  [2200 -> 1.56 -> 220 -> 22, 2026-08-24]

# ⚠ 0 -> 1.5e-3 [2026-08-27].  THE NETWORK HAS NO OTHER RasG-GTP OFF-TERM.
# 3.7 is the ONLY route RasG_GTP -> RasG_GDP and it is second order in
# `RasGAPa_cyto`; with `k_gapBasal` = 0 and `k_gapOnG` driven purely by Gbg, a
# cell at ZERO ambient cAMP has `RasGAPa_cyto` ~ 3/voxel and therefore
# essentially no hydrolysis at all, while 3.3b keeps supplying activation off
# the basal PIP3.  Measured on the rest-state solver: RasG-GTP sits at 99.88 %
# of pool and PIP3 at 38.9 % of the lipid pool WITH NO STIMULUS — the network
# latches, and `analyse/cell_rest_state.jl` correctly refuses to install an IC.
# There is no intrinsic-GTPase reaction to carry this instead, so the basal arm
# 3.4b has to.  Restoring it is a RATE change to a reaction that is already in
# `MODULE_3_REACTIONS`, not a new edge.
# ⚠ THE REST STATE IS A HARD SWITCH, NOT A GRADED LEVEL, because the Ras/PIP3
# loop is bistable: 5e-4 leaves RasG-GTP at 72.5 %, 1.5e-3 puts it at 0.17 %,
# and there is no setting in between that rests at, say, 10 %.  Pick the OFF
# branch with margin rather than trying to land on a middle that does not exist.
# ⚠ IT IS PAIRED WITH `k_pip3Ras` (3.3b) AND MUST BE RAISED WITH IT.  Raising the
# feedback alone re-opens exactly the rest-state latch this constant closes:
# measured, resting RasG-GTP 47.3 % of pool at k_pip3Ras = 55 with k_gapBasal
# left at 1.5e-3.  2.5e-3 is the smallest value restoring a clean rest state at
# that gain.  ⚠ It is NOT free to raise further "to be safe": being
# stimulus-INDEPENDENT it suppresses the DRIVEN response exactly as hard as it
# suppresses rest, and at k_pip3Ras = 70 every value that fixed rest also
# collapsed the driven gradient (F/B down to 1.2-1.3).
# The 55/0.22/2.5e-3 triple was TRIED AND REVERTED — see the k_pip3Ras note for
# the four-seed measurement that rejected it.
#
# ⚠⚠ 1.5e-3 → 0 [2026-09-06], TOGETHER WITH `DICTY_M3_BRAKE` DEFAULTING ON.
# The note above says the pre-2026-08-27 value of 0 "is the BROKEN one (no rest
# state)", and that WAS true of a network with no brake — reproduced exactly:
# at `k_gapBasal` = 0 with the brake off, resting RasG-GTP is 98.22 % of pool.
# But this term was a PATCH for the Ras/PIP3 latch that appeared when the
# antithetic (AIF) arm was deleted on 2026-08-23, and it worked by crushing the
# resting state 150-fold, which cost the module its adaptation:
#
#   `k_gapBasal` IS INPUT-INDEPENDENT, so the inhibitor stops being proportional
#   to the input, and a LEGI whose inhibitor does not track its input cannot
#   subtract it. That is the whole mechanism of adaptation.
#
# Measured, resting RasG-GTP as a fraction of the Ras pool:
#     k_gapBasal      brake OFF        brake ON
#       0             98.22 % (latch)   1.62 %
#       1.5e-4        21.82 %           1.22 %
#       1.5e-3         0.09 %           0.08 %
# The brake suppresses the SAME latch 60x without an input-independent term, and
# leaves a rest state that is off-but-populated (622 RasG-GTP and 1734 PIP3 per
# voxel, against 36 and 199 with the patch) rather than empty.
#
# Acceptance tests, paired, committed → this pair (calib/m3_kinetics_probe.jl,
# 3 seeds, and calib/m3_probe.jl):
#     PEAK time            11.67 s FAIL  →   6.50 s PASS   (window 3–8 s)
#     % amplitude @30 s     76.1 % FAIL  →  24.0 % FAIL    (t20 34 s; near miss)
#     % amplitude @60 s     59.4 % FAIL  →   7.7 % PASS
#     SINGLE (no 2nd peak)    NaN  FAIL  →  29.8 % PASS
#     kinetics criteria       1/6        →  5/6
#     step adapted plateau  54.4 % of peak (no adaptation, "check")
#                                        →   3.5 % PASS (near-perfect)
#     peak/basal            1135×        →   3.92×  (calibrated value: 3.84×)
#     resting baseline      112          →  67 809 molecules/cell
#     rest state           no pool inverted, both configs
# ⚠ ONE REGRESSION: `calib/m3_probe.jl pulse` locates the peak at 6.50 s, just
# outside ITS 3–5 s window (it is inside the kinetics probe's 3–8 s window; the
# two probes disagree on the window). Residual −2.7 % PASS and decay monotone.
# `DICTY_M3_BRAKE=0 DICTY_M3_KGAPBASAL=0.0015` restores the previous pair.
# ⚠⚠ 0 → 5e-4 [2026-09-30], TOGETHER WITH THE PKB-READ BRAKE BECOMING DEFAULT. The 0 was
# paired with the PIP3-read brake, which (being active at rest) was the only Ras off-switch
# at zero ambient: both IFFL arms scale with Gβγ and vanish there, so with the brake read by
# PKB (≈0 at rest) nothing opposes 3.3b and the rest-state solve lands on the LATCHED branch
# (RasG-GTP 90 % of pool, PIP3 236 k/voxel). A stimulus-independent RasGAP at rest is what
# the literature reports: C2GAP1 "localized on the membrane of resting cells, suggesting its
# potential role in inhibiting the basal Ras activity", and c2gapA⁻ cells have "enhanced
# basal Ras activity" (Xu et al. 2021 Front Cell Dev Biol, PMC8362602); NF1 is constitutive
# (Zhang et al. 2008). [CALIBRATED] minimum: 5e-5 latches, 2e-4 gives the low rest state
# (RasG-GTP ≈ 1, PIP3 ≈ 8 /voxel); 5e-4 is 2.5× above it (1.5 % of the GAP pool armed at
# rest). `DICTY_M3_KGAPBASAL=0` restores the previous value.
const k_gapBasal = _m3("DICTY_M3_KGAPBASAL", 5.0f-4 * _M3_TAUI_PAIR, 0.0f0)  # s⁻¹  (3.4b, basal I→I*)  [8e-5 -> 0 -> 1.5e-3 -> 0 -> 5e-4 -> 1.515e-3 (τ_I pairing, 2026-10-02)]

# 3.1 Gbg_cyto + RasGEFR_cyto → RasGEFR*_cyto + Gbg_cyto
r3_1 = Reaction([:Gbg_cyto, :RasGEFR_cyto], [:RasGEFRa_cyto, :Gbg_cyto],  k_gefRon,   nothing)
# 3.2 RasGEFR*_cyto → RasGEFR_cyto
r3_2 = Reaction([:RasGEFRa_cyto], [:RasGEFR_cyto],     k_gefRoff,             nothing)
# 3.3 RasGEFR*_cyto + RasG_GDP_mem → RasG_GTP_mem + RasGEFR*_cyto
r3_3 = Reaction([:RasGEFRa_cyto, :RasG_GDP_mem], [:RasG_GTP_mem, :RasGEFRa_cyto], k_rasGon, nothing)
# 3.4 RC_mem + RasGAP_cyto → RC_mem + RasGAP*_cyto  (delayed inhibitor)
r3_4 = Reaction([:RC_mem, :RasGAP_cyto], [:RasGAPa_cyto, :RC_mem],        k_gapOn,           nothing)

# 3.4b RasGAP_cyto → RasGAP*_cyto  (basal, receptor-independent inhibitor production)
r3_4b = Reaction([:RasGAP_cyto], [:RasGAPa_cyto], k_gapBasal, nothing)
# 3.4e Gbg_cyto + RasGAP_cyto → RasGAPa_cyto + Gbg_cyto  (LEGI inhibitor, same
#      input as the 3.1 excitor — see the k_gapOnG block above)
r3_4e = Reaction([:Gbg_cyto, :RasGAP_cyto], [:RasGAPa_cyto, :Gbg_cyto], k_gapOnG, nothing)
# 3.5 RasGAP*_cyto → RasGAP_cyto
r3_5 = Reaction([:RasGAPa_cyto], [:RasGAP_cyto],       k_gapOff,           nothing)
# 3.7 RasGAP*_cyto + RasG_GTP_mem → RasG_GDP_mem + RasGAP*_cyto
r3_7 = Reaction([:RasGAPa_cyto, :RasG_GTP_mem], [:RasG_GDP_mem, :RasGAPa_cyto], k_rasGoff, nothing)

_m3frac(pool, frac) = round(Int, frac * pool)

_gbg_uM   = Gbg_v_rest / Float32(MOLEC_PER_µM)
_fE_rest  = k_gefRon * _gbg_uM / (k_gefRon * _gbg_uM + k_gefRoff)
_fI_rest  = (k_gapOnG * _gbg_uM + k_gapBasal) /
            (k_gapOnG * _gbg_uM + k_gapBasal + k_gapOff)
const GEFRa_v_rest = max(_m3frac(GEFR_v, _fE_rest), 1)
const GAPa_v_rest  = max(_m3frac(GAP_v,  _fI_rest), 1)
_ras_on   = k_rasGon  * (GEFRa_v_rest / Float32(MOLEC_PER_µM))
_ras_off  = k_rasGoff * (GAPa_v_rest  / Float32(MOLEC_PER_µM))

# ── 3.7z — NF1 AS A SATURATED (ZERO-ORDER) GAP, written 2026-09-27 ──────────
# ⚠ OPT-IN (`DICTY_M3_ZOU=1`); replaces 3.7 in MODULE_3_REACTIONS when on.
#
# WHY.  3.7 is mass action, k·[GAP*]·[RasGTP], so the Ras cycle is a GRADED
# (hyperbolic) readout of E/I: it neither amplifies a spatial difference nor
# suppresses one.  The wave probe (calib/wave_probe.jl, 2026-09-26) measured
# exactly that — the post-wave random polarity that steers cells after every
# wave is born at Gβγ and passed through by Module 3 at 0.5–0.8× (per-cell
# correlation Gβγ↔RasGTP 0.98–1.00) instead of being rejected.
#
# THE BIOLOGY.
#   * Nakajima A, Ishihara S, Imoto D, Sawai S (2014) Nat Commun 5:5367:
#     Dictyostelium Ras activation at the leading edge is suppressed while the
#     cAMP concentration FALLS ("rectification"); it occurs at or upstream of
#     Ras, survives PI3K inhibition (LY294002) and latrunculin, and "arises
#     naturally in a single-layered incoherent feedforward circuit with ZERO-
#     ORDER ULTRASENSITIVITY" — Michaelis–Menten GEF/GAP legs, G(R) =
#     k_I·I·R/(K_I+R), with K_I = 0.01 of the Ras pool enabling it.
#   * The GAP that runs this leg is NF1 (DdNF1/nfaA): uniformly distributed,
#     essential for directional sensing, and without it Ras is not shut off
#     (Zhang S, Charest PG, Firtel RA (2008) Curr Biol 18:1587).
#   * NF1-GAP: Km = 0.3 µM, kcat = 1.4 s⁻¹ (p120GAP: Km 9.7 µM, kcat 19 s⁻¹) —
#     Wiesmüller L, Wittinghofer A (1992) J Biol Chem 267:10207, as tabulated
#     in Coyle SM, Lim WA (2016) eLife 5:e12435 (verified 2026-09-27; an
#     earlier version of this note said "kcat 5–10 s⁻¹, Scheffzek 1998" — that
#     paper does not contain the constants and the kcat was wrong).  Our membrane
#     Ras pool is RasGv/MOLEC_PER_µM ≈ 8 µM per voxel, so Km/Ras_tot ≈ 0.04:
#     the zero-order regime Nakajima needs is where the real enzyme sits, not a
#     tuning choice.  Mammalian GRD in vitro; no Dictyostelium NF1 kinetics exist.
#
# THE CONSTANTS.
#   K_NF1_µM        [SOURCED, order of magnitude] 0.2 µM (measured 0.3 µM; a
#                   secondary source quotes 0.13 µM, not verified).
#   ZOU_THETA_REST  [CALIBRATED] θ = a·Ras_tot/(kcat·[GAP*]) at rest: the
#                   resting distance to the switch (θ = 1 is threshold; below
#                   it Ras-GTP ≈ K·θ/(1−θ), above it it runs to the pool).
#   k_nf1Cat        [DERIVED] from θ_rest and the IFFL's own resting E/I, which
#                   is ambient-independent to first order (both arms ∝ Gβγ):
#                   kcat = k_rasGon/M · (E/I)_rest · RasGv / θ_rest ≈ 7.9/θ s⁻¹.
#                   ⚠ Only the PRODUCT kcat·[GAP*] enters the model, and GAP* is
#                   the IFFL's activated pool (not a measured NF1 abundance). At
#                   θ_rest = 0.5 the derived 15.8 s⁻¹ is ~11× NF1's 1.4 s⁻¹, i.e.
#                   the model's active-NF1 pool stands in for ~11× more enzyme.
#                   Treat kcat·[GAP*] as [CALIBRATED], not as literature-anchored.
# ⚠⚠ DEFAULT ON [2026-09-29] (was opt-in). Rectification and near-perfect full-cell
# adaptation (step plateau 20 % → 0.2 %, m3_probe; factorial §8). `DICTY_M3_ZOU=0`
# restores the graded 3.7.
const M3_ZOU         = get(ENV, "DICTY_M3_ZOU", "1") == "1"
const K_NF1_µM       = parse(Float64, get(ENV, "DICTY_M3_NF1K", "0.2"))
const K_NF1_MOLEC    = Float32(K_NF1_µM * MOLEC_PER_µM)
const ZOU_THETA_REST = parse(Float64, get(ENV, "DICTY_M3_ZOUTHETA", "0.5"))
# ⚠ [2026-09-30] E/I from the Gβγ-DRIVEN inhibitor fraction only. With the basal arm in
# the denominator, any k_gapBasal was cancelled exactly by a proportionally smaller kcat
# (only kcat·[GAP*] enters), so a stimulus-independent GAP could never add resting
# capacity. Identical to before at the default k_gapBasal = 0.
_fI_rest_G = k_gapOnG * _gbg_uM / (k_gapOnG * _gbg_uM + k_gapOff)
const _EI_rest       = (GEFR_v * Float64(_fE_rest)) / (GAP_v * Float64(_fI_rest_G))
const k_nf1Cat       = Float32(Float64(k_rasGon) / MOLEC_PER_µM * _EI_rest * RasGv / ZOU_THETA_REST)  # s⁻¹
# 3.7z RasGAP*_cyto + RasG_GTP_mem → RasG_GDP_mem + RasGAP*_cyto,
#      propensity kcat·[GAP*]·R/(K+R).  The engine divides a two-educt rate by
#      MOLEC_PER_µM, so the rate handed over is kcat·MOLEC_PER_µM.
r3_7z = Reaction([:RasGAPa_cyto, :RasG_GTP_mem], [:RasG_GDP_mem, :RasGAPa_cyto],
                 Float32(Float64(k_nf1Cat) * MOLEC_PER_µM),
                 @rate(x -> 1f0 / (K_NF1_MOLEC + x[:RasG_GTP_mem])))
# resting Ras-GTP with 3.7z and no feedback: a(Rt−R)(K+R) = V·R, solved exactly
_ras_rest_zou = let a = Float64(k_rasGon) / MOLEC_PER_µM * GEFR_v * Float64(_fE_rest),
                    V = Float64(k_nf1Cat) * GAP_v * Float64(_fI_rest), K = Float64(K_NF1_MOLEC),
                    b = a * (RasGv - K) - V
    (b + sqrt(b^2 + 4a^2 * RasGv * K)) / (2a)
end

const RasGTP_v_rest = M3_ZOU ? round(Int, _ras_rest_zou) :
                      _m3frac(RasGv, _ras_on / (_ras_on + _ras_off))

const M3_INITIAL = Dict{Symbol, Int}(
    :RasGEFR_cyto  => GEFR_v - GEFRa_v_rest, :RasGEFRa_cyto => GEFRa_v_rest,
    :RasG_GDP_mem  => RasGv - RasGTP_v_rest,
    :RasG_GTP_mem  => RasGTP_v_rest,
    :RasGAP_cyto   => GAP_v - GAPa_v_rest,
    :RasGAPa_cyto  => GAPa_v_rest,
    # ⚠ THE STEN BRAKE'S SUBSTRATE, ADDED 2026-09-06 WITH THE 3.3c/d/e EDGES.
    # `RasBrake_cyto` had a diffusion constant and a pool constant but NO INITIAL
    # CONDITION anywhere, so it started at 0 — and nothing in the network
    # produces it. A brake written without this line can never activate, and the
    # comparison "brake on vs brake off" would then have come out as "no effect"
    # for a reason that has nothing to do with the mechanism. That is this
    # repo's own recorded trap (a switch whose gate species has no producer and
    # a zero IC is a no-op that still looks like a tested option), and it was
    # caught here only because the scan printed `RasBrakea = 0.00` and the
    # difference between "not recruited" and "no substrate exists" was checked.
    # Resting state is ALL cytosolic: the brake is recruited by PIP3 (3.3c),
    # and a resting cell has essentially none.
    :RasBrake_cyto => RasBrake_v,
    :RasBrakea_mem => 0,
)

const D_M3_GEFR  = _m3("DICTY_M3_DGEFR",  0.1f0, 0.1f0)   # µm²/s — LOCAL excitor
const D_M3_GAP   = _m3("DICTY_M3_DGAP",   20.0f0, 10.0f0) # µm²/s — GLOBAL inhibitor
# ⚠⚠ 0.4 → 20.0 [2026-09-10]: THE BRAKE IS NOW GLOBAL, AND THAT REVERSES THE
# DESIGN INTENT RECORDED IMMEDIATELY BELOW. The old value was chosen so the
# brake "cannot act as a second LEGI-style global inhibitor" — it was meant to
# be LOCAL, terminating the patch it sits in. MEASURED CONSEQUENCE (mutual
# information about the true gradient direction, 4 tuned vs 2 legacy seeds,
# 1000 cells, bias-corrected, all frames):
#
#   transition                     tuned(local brake)   legacy(no brake)
#   M2 -> RasGEFRa_cyto                  0.819               0.820
#   RasGEFRa_cyto -> RasG_GTP_mem        0.263               0.629   <-- ONLY
#   RasG_GTP_mem -> PIP3_mem             0.388               0.312       DIFF
#
# Every other stage is identical to three digits; the entire information loss
# sits on the one edge the brake acts on (3.3e). The mechanism is the one this
# file already states for `RasGAPa_cyto` at D = 10 — a LOCAL inhibitor "cannot
# create a near/far asymmetry", and read the other way, a local inhibitor
# FOLLOWS the asymmetry and cancels it: `RasBrakea_mem` is activated by PIP3,
# PIP3 is the OUTPUT of the 3.3b loop, so the brake is strongest exactly where
# the signal is strongest, i.e. at the front. It buys adaptation by deleting the
# spatial contrast. A globally-mixed brake subtracts the same amount everywhere
# and leaves the contrast standing — which is what a LEGI inhibitor is for.
#
# ⚠ At REST this changes nothing: the cell is uniform, so local and global
# diffusion give the identical rest state and the identical latch suppression
# (1.62 % of pool). The two differ ONLY inside a patch, which is exactly the
# regime the mean-field derivation below could not see, and the block below
# admits as much ("the mean field cannot see the thing a brake is FOR").
#
# ⚠ PHYSICAL CAVEAT, NOT RESOLVED: `RasBrakea_mem` is named as a MEMBRANE
# species, and 20 µm²/s is not a membrane protein's diffusion constant (those
# run 0.1–1). The physically honest form of this change is a CYTOSOLIC active
# brake (`RasBrakea_cyto`), which is how `RasGAPa_cyto` legitimately gets D = 20.
# The rename is deliberately NOT done here so that this commit changes exactly
# one number and stays A/B comparable; if the mechanism survives its tests, the
# rename is the cleanup.
# `DICTY_M3_DBRAKE=0.4` restores the local brake. Old rationale, still the
# reason the constant exists: Modules.pdf's own D_PKB*s = 0.4 µm²/s (Table 4).
const D_M3_BRAKE = _m3("DICTY_M3_DBRAKE", 20.0f0, 0.4f0)  # µm²/s — GLOBAL (was 0.4, LOCAL)
const M3_DIFFUSION = Dict{Symbol, Float32}(
    :RasGEFR_cyto  => D_M3_GEFR, :RasGEFRa_cyto => D_M3_GEFR,  # LOCAL excitor
    :RasGAP_cyto   => D_M3_GAP,  :RasGAPa_cyto  => D_M3_GAP,   # GLOBAL inhibitor (LEGI)
    :RasG_GDP_mem  => 0.1f0,     :RasG_GTP_mem  => 0.1f0,
    :RasBrake_cyto => D_M3_BRAKE, :RasBrakea_mem => D_M3_BRAKE, # dedicated STEN brake, NOT the LEGI arm
)

# ============================================================================
# Module 4 — PIP3 Generation and Spatial Confinement
# ============================================================================

# regression can be bisected to "Module 4 retune" vs. anything else.

const M4_SPEED = parse(Float64, get(ENV, "DICTY_M4_SPEED", "20.0"))
# ⚠ `get(ENV, env, new)` returns a STRING when the variable IS set, so the old
# one-liner `M4_SPEED * get(ENV, env, new)` threw MethodError (Float64 * String)
# for EVERY `DICTY_M4_*` rate override — i.e. none of them had ever been usable.
_m4s(env, new) = M4_SPEED * (haskey(ENV, env) ? parse(Float64, ENV[env]) : Float64(new))

const PIP2_SCALE = parse(Float64, get(ENV, "DICTY_M4_PIP2V", "769231")) / 769.0
const PI3K_SCALE = parse(Float64, get(ENV, "DICTY_M4_PI3KV", "200000")) / 154.0
const PTEN_SCALE = parse(Float64, get(ENV, "DICTY_M4_PTENV", "10000"))  / 154.0

const k_pi3kOn  = _m4s("DICTY_M4_KPI3KON",  0.02f0) # µM⁻¹s⁻¹ (4.1) [2.5 -> 0.02, 2026-08-08]
const k_pi3kOff = _m4s("DICTY_M4_KPI3KOFF", 0.5f0)           # s⁻¹     (4.2)

const k_pip3 = _m4s("DICTY_M4_KPIP3", 6.2f0) / Float32(PI3K_SCALE)  # µM⁻¹s⁻¹ (4.3), catalyst-compensated  [1.0 -> 1.6 -> 6.2]

const k_ptenOn   = _m4s("DICTY_M4_KPTENON",   30.0f0) / Float32(PIP2_SCALE)  # µM⁻¹s⁻¹ (4.4), catalyst-compensated
const k_ptenOff  = _m4s("DICTY_M4_KPTENOFF",  0.2f0)           # s⁻¹     (4.5)
const k_ptenDisp = _m4s("DICTY_M4_KPTENDISP", 1200.0f0) / Float32(PIP2_SCALE)  # µM⁻¹s⁻¹ (4.6), catalyst-compensated

const k_pten = _m4s("DICTY_M4_KPTEN", 4.5f0) / Float32(PTEN_SCALE)  # µM⁻¹s⁻¹ (4.7), catalyst-compensated

const k_ship  = _m4s("DICTY_M4_KSHIP",  4.0f0)          # µM⁻¹s⁻¹ (4.8) [1.0 -> 4.0, 2026-08-08]
const k_inpp4 = _m4s("DICTY_M4_KINPP4", 2.0f0)          # µM⁻¹s⁻¹ (4.9)

# [REMOVED 2026-08-27] `SHIP_like`/`INPP4_like` as explicit species. Each was
# its own trivial one-species conservation law (nothing ever produced or
# consumed it — see the old `77 = SHIP_like` / `77 = INPP4_like` moiety laws),
# so in 4.8/4.9 it appeared as an UNCHANGED catalyst on both sides of a mass-
# action reaction: propensity = k_bimolecular · [PIP3] · [enzyme], and since
# [enzyme] never moves, `k_bimolecular · [enzyme]` is a constant, not a
# species-dependent quantity — collapsing it into a single pseudo-first-order
# rate constant is an EXACT algebraic simplification, not an approximation.
# Was already the reasoning `PI34P2_v_rest` used below (`k_ship * _uM(SHIP_v) /
# (k_inpp4 * _uM(INPP4_v))`), just not yet applied to the reactions themselves.
# `SHIP_v`/`INPP4_v` are kept as the constants that fix the rate, not as species.
const SHIP_v  = parse(Int, get(ENV, "DICTY_M4_SHIPV",  "77"))   # constant pool, molec/voxel
const INPP4_v = parse(Int, get(ENV, "DICTY_M4_INPP4V", "77"))   # constant pool, molec/voxel
const k_ship_eff  = k_ship  * Float32(SHIP_v)  / Float32(MOLEC_PER_µM)  # s⁻¹ (4.8), PIP3_mem   → PI34P2_mem
const k_inpp4_eff = k_inpp4 * Float32(INPP4_v) / Float32(MOLEC_PER_µM)  # s⁻¹ (4.9), PI34P2_mem → PIP2_mem

# 4.1 RasG_GTP_mem + PI3K_cyto → RasG_GTP_mem + PI3K_mem
r4_1 = Reaction([:RasG_GTP_mem, :PI3K_cyto], [:PI3K_mem, :RasG_GTP_mem],  k_pi3kOn,   nothing)
# 4.2 PI3K_mem → PI3K_cyto
r4_2 = Reaction([:PI3K_mem], [:PI3K_cyto],                                k_pi3kOff,  nothing)
# 4.3 PIP2_mem + PI3K_mem → PIP3_mem + PI3K_mem
r4_3 = Reaction([:PIP2_mem, :PI3K_mem], [:PIP3_mem, :PI3K_mem],           k_pip3,     nothing)
# 4.4 PTEN_cyto + PIP2_mem → PTEN_mem + PIP2_mem  (PIP2 = docking site)
r4_4 = Reaction([:PTEN_cyto, :PIP2_mem], [:PTEN_mem, :PIP2_mem],          k_ptenOn,   nothing)
# 4.5 PTEN_mem → PTEN_cyto
r4_5 = Reaction([:PTEN_mem], [:PTEN_cyto],                                k_ptenOff,  nothing)
# 4.6 PTEN_mem + PIP3_mem → PTEN_cyto + PIP3_mem  (mutual inhibition)
r4_6 = Reaction([:PTEN_mem, :PIP3_mem], [:PTEN_cyto, :PIP3_mem],          k_ptenDisp, nothing)
# 4.6c COOPERATIVE eviction — the Hill-gated alternative to 4.6.

const M4_COOP  = get(ENV, "DICTY_M4_COOP", "1") == "1"

const N_DISP   = parse(Int, get(ENV, "DICTY_M4_NDISP", "2"))   # [4 -> 2, 2026-08-20]
# ⚠ 0.05 -> 0.17 [2026-08-27].  THE 4.6c EVICTION GATE WAS SATURATED, so
# Module 4's mutual-inhibition amplifier contributed no contrast at all: at
# K = 0.05 the gate reads 0.853 (front) / 0.833 (back) against a PIP3 that
# differs by 7.8 %, i.e. the switch was pinned near 1 and could not report the
# difference.  Same rule as "check the activated fraction of every switch
# before tuning its rates", applied to a Hill gate.  K = 0.17 puts the gate at
# ~0.18, on its steep flank, and the front/back PIP3 ratio goes 1.078 -> 3.58.
# ⚠ THIS PARAMETER HAS A CLIFF AT ~0.195 AND MUST NOT BE RAISED PAST IT.
# Measured on calib/tune_2026-08-27/fb.jl, driven at 60 nM, solving from a
# right-way AND a wrong-way initial split:
#     0.14  1.47 / 1.47   unique
#     0.17  3.58 / 3.58   unique   <-- committed
#     0.19  4.19 / 4.19   unique
#     0.20  0.40 / 0.63   HYSTERETIC and sign-INVERTED
#     0.30  4.72 / 0.38   HYSTERETIC (holds F/B 3.9 with NO input at all)
# Past the cliff the PIP3 domain is a winner-take-all latch whose orientation is
# set by the initial condition rather than by the gradient — it still LOOKS like
# a strong gradient from a single solve, which is exactly why the acceptance
# test has to be run from both starts.
# ⚠ 0.05 -> 0.17 [2026-08-27].  At 0.05 the 4.6c eviction gate was SATURATED
# (0.853 front / 0.833 back against a PIP3 that differs by 7.8 %), so Module 4's
# mutual-inhibition amplifier could not report the difference at all — the same
# rule as "check the activated fraction of every switch before tuning its rates",
# applied to a Hill gate.  0.17 puts the gate at ~0.20, on its steep flank.
# ⚠ HARD CLIFF AT ~0.195 — DO NOT RAISE PAST IT.  Measured on
# calib/tune_2026-08-27/fb.jl (driven, solved from a right-way AND a wrong-way
# initial split; a single solve cannot see this):
#     0.14  1.47 / 1.47   unique
#     0.17  3.58 / 3.58   unique   <-- committed
#     0.19  4.19 / 4.19   unique
#     0.20  0.40 / 0.63   HYSTERETIC and sign-INVERTED
#     0.30  4.72 / 0.38   HYSTERETIC (holds F/B 3.9 with NO input at all)
# Past the cliff the PIP3 domain is a winner-take-all latch oriented by the
# initial condition rather than by the gradient — and it still LOOKS like a
# strong gradient from one solve.
# A 0.22 setting was TRIED AND REVERTED with the 3.3b raise; see k_pip3Ras.
const K_DISP_F = parse(Float64, get(ENV, "DICTY_M4_KDISP", "0.17"))
# PIP2_v itself is declared below the reaction block, so read its ENV source
# here the same way PIP2_SCALE does rather than reordering the file.
const K_DISP_N = Float32((K_DISP_F *
    parse(Int, get(ENV, "DICTY_M4_PIP2V", "769231")))^N_DISP)   # K^n, molecules^n
# Rate is the SATURATED eviction rate; the gate supplies the PIP3 dependence
# that 4.6 carried in its second educt, so the two forms are the same process.
const k_ptenDispC = _m4s("DICTY_M4_KPTENDISPC", 40.0f0)   # s⁻¹ (4.6c)
r4_6c = Reaction([:PTEN_mem], [:PTEN_cyto], k_ptenDispC,
                 @rate(x -> x[:PIP3_mem]^N_DISP / (K_DISP_N + x[:PIP3_mem]^N_DISP)))
# 4.7 PIP3_mem + PTEN_mem → PIP2_mem + PTEN_mem
r4_7 = Reaction([:PIP3_mem, :PTEN_mem], [:PIP2_mem, :PTEN_mem],           k_pten,     nothing)

# ── 4.8  PIP3_mem → PIP2_mem   (PTEN-INDEPENDENT 5-phosphatase removal) ─────
# [ADDED 2026-08-27]  ⚠ THIS REACTION WAS MISSING FROM THE NETWORK ENTIRELY.
# The 2026-08-27 "fold-in to pseudo-first-order" note above `SHIP_v`/`INPP4_v`
# describes 4.8/4.9 as PERMANENT and derives `k_ship_eff`/`k_inpp4_eff` and
# `PI34P2_v_rest`/`M4_INITIAL` from them — but the two `Reaction` objects were
# never written, and `MODULE_4_REACTIONS` never listed them.  So PIP3's
# dominant sink (CLAUDE.md: "SHIP does 87 % of PIP3 removal") did not exist,
# leaving PTEN as the ONLY route out of PIP3 — and PTEN is EVICTED by PIP3
# (4.6c).  That is an unopposed runaway, and it is what the rest-state solver
# reports: with the 3.3b amplifier fully OFF the resting cell still converts
# 60.6 % of the lipid pool to PIP3 and strips PTEN to 4.2 % of its pool.  The
# mean-field balance reproduces that number exactly when `k_ship_eff` is
# omitted (prod/PIP2 = 0.184 s⁻¹ vs PTEN-only removal 0.120 s⁻¹ ⇒ 60.5 %), so
# the missing reaction is the whole of the discrepancy.  Same defect class as
# the dead `r3_3b`, the commented-out `r2_4` and the `r3_3c-e` the Module 3
# list names but does not define: prose describing a live mechanism that the
# reaction list does not contain.
#
# COLLAPSED TO ONE STEP, by decision: the PI(3,4)P2 intermediate is skipped and
# PIP3 returns straight to PIP2.  In the two-step chain the flux OUT of PIP3 is
# set by the FIRST step alone, so a single first-order `k_ship_eff` sink
# preserves PIP3's removal rate exactly; only `PI34P2_mem`'s own level is lost.
# ⚠ `PI34P2_mem` therefore has no producer and stays at 0.  It is deliberately
# KEPT as a species (M4_DIFFUSION/CELL_DIFF) because eight probes and drivers
# track it by name (calib/m4_probe.jl, calib/m4_polarity_probe.jl,
# dicty_simulation.jl, …) and would error on an unknown species; it is simply
# no longer given any molecules.  This also returns the ~61 700 molec/voxel
# that `M4_INITIAL` used to park in it to the usable PIP2 pool, so the lipid
# conservation law becomes the full `PIP2_v` = 769 231 rather than 707 547.
const k_pip3Dec = haskey(ENV, "DICTY_M4_KPIP3DEC") ?
    parse(Float64, ENV["DICTY_M4_KPIP3DEC"]) : Float64(k_ship_eff)   # s⁻¹
r4_8 = Reaction([:PIP3_mem], [:PIP2_mem], Float32(k_pip3Dec), nothing)


# ============================================================================
# 4.10-4.14 — the DELAYED NEGATIVE FEEDBACK that makes Module 4 EXCITABLE
# rather than either bistable-and-latching or monostable-and-soft.
# [2026-08-08]
# ----------------------------------------------------------------------------
# WHY THIS ARM EXISTS. Mutual inhibition alone (4.3-4.7) gives exactly two
# unsatisfactory regimes, both measured (calib/CALIBRATION_M4.md §6):
#   loop gain > 1  sharp domain (pol_x 0.89) but it LATCHES — self-polarises with
#                  no stimulus, invents a domain along y, never resets (1.175)
#   loop gain < 1  resets (0.14) and cannot break symmetry (pol_y 0.00) but is
#                  soft (pol_x 0.42), so the back keeps ~40 % of the front's PIP3
# An EXCITABLE network gets both: locally supra-threshold positive feedback for
# the sharp all-or-none patch, plus a DELAYED negative feedback that terminates
# it and imposes a refractory period, so it cannot latch.
#
# THE BIOLOGY (this is the documented Dictyostelium STEN mechanism, not a
# generic relaxation oscillator):
#   * The STEN's positive feedback is the Ras/Rap ↔ PIP2 mutual inhibition
#     already in 4.3-4.7. Its negative feedback is "delayed PKB activation by
#     PIP3" acting back on the active state, and the refractory state carries
#     LOWER Ras/Rap activity together with HIGHER PI(3,4)P2 and PKB activity
#     (Devreotes/Iglesias STEN reviews; Miao et al. 2019 Wave patterns;
#     Fukushima, Matsuoka & Ueda 2019 J Cell Sci 132:jcs224121, already cited
#     for the 3.3b PIP3→Ras arm).
#   * PKB activation is a COINCIDENCE DETECTOR, which is where the delay comes
#     from: PKBA needs (i) PIP3 for PH-domain membrane recruitment and (ii)
#     TORC2-mediated hydrophobic-motif phosphorylation, and TORC2 is itself
#     activated by RasC at the front. Two sequential steps ⇒ a genuine lag
#     behind PIP3 rather than an exponential that tracks it.
#   * The feedback target is the RasGEF machinery: Charest PG, Shen Z, Lakoduk A,
#     Sasaki AT, Briggs SP, Firtel RA (2010) "A Ras signaling complex controls
#     the RasC-TORC2 pathway and directed cell migration", Dev Cell 18:737-749 —
#     PKB and PKBR1 PHOSPHORYLATE Sca1, the scaffold of the Sca1/RasGEF(Aimless)/
#     PP2A complex, and thereby control that complex's MEMBRANE LOCALISATION and
#     so RasC activity, explicitly "in a negative feedback fashion".
#
# WHAT IS LUMPED, stated plainly. The literature target is Sca1/RasGEF → RasC →
# TORC2. This model has no Sca1 and its RasGEF lives in Module 3, which is
# separately calibrated and adapts correctly; wiring the feedback there would
# re-open that calibration. So the feedback is applied to Module 4's own
# PI3K_mem, which RasG-GTP recruits (4.1) — i.e. "PKB ⊣ Sca1/RasGEF ⇒ less
# Ras-GTP ⇒ less membrane PI3K" is collapsed into its Module 4 endpoint. The
# SIGN, the LOCALITY and the DELAY are all as published; only the intermediate
# is elided. Restoring the faithful route means a PKBAa_mem ⊣ RasGEFRa_cyto edge
# in Module 3 and a re-run of CALIBRATION_M3.md.
#
# SPECIES NAMING [CORRECTED 2026-08-27]. This arm is **PKBA** (PkbA / "AKT",
# DDB_G0287449) and nothing else, so it is now spelled `PKBA_cyto` / `PKBA_mem` /
# `PKBAa_mem` (was `PKB_cyto` / `PKB_mem` / `PKBa_mem`).  Module 8b's
# `PKBR1_mem`/`PKBR1a_mem` is the OTHER gene product, PkbR1 (DDB_G0282133) — a
# genuinely different protein, not a second parameterisation of this one:
#   * PkbA carries an N-terminal PIP3-specific PH domain; its membrane
#     recruitment AND its subsequent activation-loop phosphorylation both require
#     PIP3 ("inhibition of PtdIns(3,4,5)P₃ signalling blocks AKT recruitment to
#     membrane compartments and consequent phosphorylation by both PDK1 and
#     TORC2" — Kamimura, Tang & Devreotes 2010, J Cell Sci 123:983).  That is
#     4.10 + 4.12: recruitment first, phosphorylation second — the coincidence
#     detector, and where this arm's delay comes from.
#   * PkbR1 is myristoylated, constitutively membrane-anchored, and activated
#     independently of PIP3 (Meili et al. 2000 Curr Biol 10:708; Kamimura et al.
#     2008 Curr Biol 18:1034).  That is 8b.1: a single TORC2-gated step with no
#     recruitment stage.
# The two therefore differ in EXACTLY one place — whether a PIP3 recruitment step
# precedes the TORC2 step — and having different pool sizes (10⁴/voxel here vs
# `PKBR1_v` = 77/voxel) is expected of two different proteins, not a conflict.
#
# ⚠ THE OLD NOTE HERE WAS STALE AND SAID THE OPPOSITE ("the SAME molecules
# Module 8's 8.7/8.8 use — deliberately, since it is the same kinase", plus
# "Module 8 is currently commented out").  Neither has been true since
# 2026-08-09: Module 8 is ON, and 8.7/8.8/8.13 are explicitly EXCLUDED from
# `MODULE_8_REACTIONS` precisely because they were a second, PIP3-INDEPENDENT
# activation route for this pool — i.e. they modelled PkbR1's mechanism while
# consuming PkbA's species.  Module 8b already supplies PkbR1 properly, so 8.7
# is redundant as well as mis-typed and stays out.  Reading the old note is what
# makes a reviewer conclude "Module 4 and Module 8 model the two paralogs and
# collide"; the split is already done — only the names were lying.
#
# NOT MODELLED, stated plainly: PDK1.  Kamimura 2010's central result is that
# TORC2 hydrophobic-motif phosphorylation is "wholly insufficient" for either
# kinase and that PDK1-site (activation-loop) phosphorylation is also required.
# 4.12 and 8b.1 both lump AL+HM into one step.  The lump is safe for PkbA
# (PDK1 access is itself PIP3-gated, and 4.10 already carries that gate) and is
# a genuine simplification for PkbR1 (whose PDK1 step is PIP3-independent).
# ⚠ 10^4 -> 10^5 [2026-08-29].  THE BRAKE ONLY EVER RECRUITS ~2.5 % OF ITS POOL
# (PKBAa_mem peaked at 250/voxel out of 10 000), so ~40x of its capacity was
# idle and the arm could not pull PIP3 down however large `k_pkb4Fb` was made.
# Authority is the PRODUCT k_pkb4Fb * PKBA_v, and `k_pkb4Fb` is NOT compensated
# by any PKBA scale — so the pool is a real, independent lever, and the better
# one: it raises PKBAa itself instead of a per-molecule rate that was already
# dominant over k_pi3kOff.  Measured (calib/m8_probe.jl single, seed 1):
#     10k pool, fb 400 : export 9.75e6, PIP3 @480 s = 4.8x basal
#    100k pool, fb 160 : export 9.32e6, PIP3 @480 s = 3.1x basal   <-- committed
# Same authority at a 2.5x SMALLER and far more defensible bimolecular constant,
# with 5.5x more molecules carrying it (less shot noise in the brake itself).
const PKBA_v = parse(Int, get(ENV, "DICTY_M4_PKBV", "100000"))   # PKBA per voxel
#
# ── THE RESTING-DRIVE ANCHOR ("the RasG-GTP / GEF problem") [2026-08-27] ─────
# 4.10 and 4.12 are the only two edges in Module 4 whose rate is multiplied by a
# species OWNED BY ANOTHER, SEPARATELY-CALIBRATED MODULE (`PIP3_mem`, and
# `RasG_GTP_mem` from Module 3).  Every Module 3 re-pairing since 2026-08-08 was
# built to hold Module 3's OWN products invariant — `k_rasGon·GEFR_v`,
# `k_rasGoff·GAP_v`, `f_E`, `f_I`, τ_E, τ_I — and NONE of them held
# `k_pkb4Act·[RasG_GTP]`.  So this arm was silently re-scaled by work that was
# explicitly designed not to change any dynamics:
#     RasGv 385 → 38500, GEFR_v 154 → 4620 → 46200 → 138600, GAP_v 154 → 46200,
#     k_rasGon 200 → 0.275 → 4.0 → 1.3333 → …, k_rasGoff 2200 → …
#   ⇒ `RasGTP_v_rest`  ~854–956  →  5356 molec/voxel   (≈5.6×)
# and `k_pip3` 1.6 → 6.2 (2026-08-21) did the same to 4.10 through `PIP3_mem`.
# CONSEQUENCE at the committed constants (k = 0.026 / 0.043 fixed), computed
# from this file's own numbers against the MEASURED resting `PIP3_mem` (7298/vox)
# and `RasGTP_v_rest` (5356/vox): the resting activation rate
# `k_pkb4Act·[RasG_GTP]` is 0.0478 s⁻¹ against `k_pkb4Off` = 0.05 s⁻¹, so ~49 %
# of recruited PKBA is activated AT REST, and the steady state of the 4.10–4.13
# cycle puts `PKBAa_mem` at **40.7 % of pool (4070/voxel)** — against the
# 455/voxel (4.6 %) that CALIBRATION_M4.md measured and CALIBRATION_M8.md §4
# quotes.  An **8.9× over-charge of a DELAYED NEGATIVE FEEDBACK, at rest**, in
# the module whose recorded defects are "does not return to rest" and
# "M4_INITIAL is still a startup transient".  r4_14's resting PI3K-eviction rate
# goes 0.104 → 0.845 s⁻¹ with it, i.e. 1.0 % → 8.4 % of `k_pi3kOff` = 10 s⁻¹.
# (Reproduce either side with `DICTY_M4_PKBANCHOR=0` / `=1`.)
#
# FIX, same discipline as `k_reassoc` and the ICs: HOLD THE MEASUREMENT, DERIVE
# THE CONSTANT.  The calibrated quantity is the resting per-molecule rate, not
# the bimolecular constant, so the two rates are now solved from the anchors
# below against whatever rest state Modules 3/4 currently have.  Back-solved from
# the 2026-08-08 calibration (PKBAa_mem 4.6 % of pool at rest, PIP3_mem 1526/vox,
# RasG_GTP_mem ~854–956/vox):
const PKB4_ANCHOR   = get(ENV, "DICTY_M4_PKBANCHOR", "1") == "1"
# ⚠ 0.00824 -> 0.002 [2026-08-29].  DELAYING THE BRAKE IS WHAT PROTECTS THE
# cAMP PULSE.  Strength and timing are separate knobs and must be set as a pair:
# a strong brake with the original recruitment rate bites DURING the pulse and
# clips the export (fb150 suppressed PIP3 71 % by t=140 and cut export to
# 7.64e6), while the same strength with 4x slower recruitment acts after the
# pulse and keeps export at 9.3-9.8e6.  This is the "inert for the first ~30 s,
# then reset" requirement expressed as a rate.
const PKB4_ON_REST  = parse(Float64, get(ENV, "DICTY_M4_PKBONREST",  "0.002"))  # s⁻¹, 4.10 at rest
const PKB4_ACT_REST = parse(Float64, get(ENV, "DICTY_M4_PKBACTREST", "0.00853")) # s⁻¹, 4.12 at rest
# ⚠ `PKB4_ANCHOR=0` restores the literal 2026-08-08 constants (0.026 / 0.043) and
# with them the 13× over-charge — it exists so the pre-fix network is
# reproducible, not because it is a supported configuration.
# ⚠ NOT ACCEPTANCE-MEASURED.  The anchor values are back-solved from the
# published rest fraction, not re-measured on the engine.  Run
# `calib/m4_polarity_probe.jl` + `calib/m4_kinetics_probe.jl` and re-check
# `PKBAa_mem` at rest before quoting any Module 4 number taken after this change.
#
# The 4.10 reference level is the MEASURED resting `PIP3_mem`, not this file's
# derived `PIP3_v_rest` — (i) `PIP3_v_rest` is declared below the reaction block
# (it needs `PIP2_v`), and (ii) it is a known startup transient: it evaluates to
# ~30 800/voxel against the 7 298/voxel `CALIBRATION_M4_PIP3.md` measures on a
# no-stimulus run, because it is derived from the PTEN docking equilibrium while
# `k_ptenDispC·H` ≫ `k_ptenOff` actually holds the variable
# (README_M3_M4_STRUCTURAL_ANALYSIS.md).  Anchoring against a number that is
# itself 4× wrong would defeat the point.  ⚠ Re-measure after any `k_pip3`,
# `k_ship` or PTEN change.
const PIP3_v_rest_ref = parse(Int, get(ENV, "DICTY_M4_PIP3REST", "7298"))  # molec/voxel, MEASURED
# 4.10  PIP3_mem + PKBA_cyto → PKBA_mem + PIP3_mem   (PH-domain recruitment; LOCAL
#       — this is what makes the feedback terminate the patch where the patch is)
const k_pkb4On    = haskey(ENV, "DICTY_M4_KPKBON") ?
    parse(Float64, ENV["DICTY_M4_KPKBON"]) :
    (PKB4_ANCHOR ? PKB4_ON_REST * MOLEC_PER_µM / max(PIP3_v_rest_ref, 1) : 0.026)  # µM⁻¹s⁻¹
# 4.11  PKBA_mem → PKBA_cyto                        (release)
const k_pkb4Off   = parse(Float64, get(ENV, "DICTY_M4_KPKBOFF",   "0.05"))   # s⁻¹
# 4.12  PKBA_mem + <gate> → PKBAa_mem + <gate>       (the TORC2 hydrophobic-motif
#       step — the second, kinase-gated half of the coincidence detection, and
#       the second stage of the delay)
#
# ⚠ WHICH SPECIES IS THE GATE.  The literature is unambiguous that PKBA's HM
# kinase is TORC2 and that TORC2's Ras activator is **RasC, not RasG**: "RasG is
# the main mediator of PI3K activation in response to cAMP, while RasC is the
# main activator of TORC2" (Cai et al. 2010, J Cell Biol 190:233; Charest et al.
# 2010, Dev Cell 18:737).  Module 8 gets this right (8.5 reads `RasC_GTP_mem`);
# this edge does not — it reads `RasG_GTP_mem`, which the old comment defended
# with "TORC2 is Ras-activated", true only of the wrong Ras.
#

const PKB4_GATE_REST    = 81
const k_pkb4Act   = haskey(ENV, "DICTY_M4_KPKBACT") ?
    parse(Float64, ENV["DICTY_M4_KPKBACT"]) :
    (PKB4_ANCHOR ? PKB4_ACT_REST * MOLEC_PER_µM / max(PKB4_GATE_REST, 1) : 0.043)  # µM⁻¹s⁻¹
# 4.13  PKBAa_mem → PKBA_cyto                       (sets the REFRACTORY duration)
const k_pkb4Deact = parse(Float64, get(ENV, "DICTY_M4_KPKBDEACT", "0.02"))   # s⁻¹
# 4.14  PKBAa_mem + PI3K_mem → PKBAa_mem + PI3K_cyto  (the Charest 2010 feedback,
#       lumped onto PI3K as explained above)
# ⚠ 1.0 -> 160 [2026-08-29].  At 1.0 the arm was STRUCTURALLY UNABLE to act:
# its ceiling is k_pkb4Fb*PKBA_v/M = 2.08 /s against k_pi3kOff = 10 /s, i.e. a
# 21 % bump to PI3K's off-rate even with the whole pool activated.  Root cause is
# a missed pairing — `k_pi3kOff` IS multiplied by M4_SPEED (=20) and `k_pkb4Fb`
# is NOT, so speeding Module 4 up left its brake 20x behind.
const k_pkb4Fb    = parse(Float64, get(ENV, "DICTY_M4_KPKBFB", "160.0"))  # µM⁻¹s⁻¹
# 4.14b PKBR1a_mem + PI3K_mem → PKBR1a_mem + PI3K_cyto  [ADDED 2026-08-27]
#       The Charest 2010 feedback is attributed to BOTH kinases, not to PkbA
#       alone: "PKB and PKB-related PKBR1 phosphorylate Sca1 and regulate the
#       membrane localization of the Sca1/RasGEF/PP2A complex, and thereby RasC
#       activity, in a negative feedback fashion", and pkbR1⁻/pkbA⁻ double nulls
#       show reduced Sca1 phosphorylation with elevated basal and post-stimulus
#       RasC (Charest et al. 2010, Dev Cell 18:737-749).  4.14 carried only
#       `PKBAa_mem`, so the larger of the two activities contributed nothing to
#       the loop.  Same rate constant, because Sca1 is a shared substrate and
#       nothing in the source separates the two kinases' gains.
#       MAGNITUDE, so this is not mistaken for a tuning change: `PKBR1a_mem`
#       measures 525 → 516 molec/CELL (≈40/voxel, flat) in CALIBRATION_M8.md §4,
#       so 4.14b adds ~0.008 s⁻¹ against `k_pi3kOff` = 10 s⁻¹ — under 0.1 %, i.e.
#       numerically ~a no-op TODAY.  It is in because it becomes the correct
#       wiring the moment Module 8's RasC arm is fixed (see the 4.12 gate note),
#       and leaving it out would hide that dependency.
#       ⚠ It is OMITTED when Module 8 is off, because `PKBR1a_mem` then has no
#       producer and would enter the species list as a permanent zero.
const M4_PKBR1_FB = get(ENV, "DICTY_M4_PKBR1FB", "1") == "1" &&
                    get(ENV, "DICTY_MODULE8", "1") == "1"
r4_10 = Reaction([:PIP3_mem, :PKBA_cyto], [:PKBA_mem, :PIP3_mem],        k_pkb4On,    nothing)
r4_11 = Reaction([:PKBA_mem], [:PKBA_cyto],                              k_pkb4Off,   nothing)
r4_12 = Reaction([:PKBA_mem, :TORC2a_mem], [:PKBAa_mem, :TORC2a_mem], k_pkb4Act, nothing)
r4_13 = Reaction([:PKBAa_mem], [:PKBA_cyto],                             k_pkb4Deact, nothing)
r4_14 = Reaction([:PKBAa_mem, :PI3K_mem], [:PKBAa_mem, :PI3K_cyto],      k_pkb4Fb,    nothing)
# ⚠ DEFAULT FLIPPED 1 -> 0 [2026-08-27].  THE COMMITTED M3/M4 CALIBRATION IS A
# BRAKE-OFF CALIBRATION, and leaving the default at 1 meant that every driver
# that did NOT set `DICTY_M4_REFRACTORY` loaded a configuration in which nothing
# had been measured.  The whole 2026-08-27 tuning round (k_pip3Ras 55 / K_DISP
# 0.22 / k_gapBasal 2.5e-3, and the solved rest state installed by DICTY_IC=auto)
# was validated with this arm OFF.
# ⚠⚠ THE ORIGINAL JUSTIFICATION FOR THIS FLIP WAS WRONG AND IS RETRACTED
# [2026-08-29].  It claimed the arm was "measured WORSE" because turning it on
# blew up the y-null (PTEN y-null -0.019 -> -0.173).  That measurement was taken
# with `calib/m4_polarity_probe.jl`, whose REDUCED network (M12+M3+M4) does not
# contain Module 8 — and r4_12's gate `TORC2a_mem` is a Module 8 species that
# appears in 4.12 only as a CATALYST, has no producer there, and starts at 0.
# So r4_12's propensity is identically zero, `PKBAa_mem` never forms, 4.14 never
# fires, and THE BRAKE CANNOT ACT IN THAT PROBE AT ALL.  Verified directly:
# brake ON vs OFF give PIP3 peak 93453 vs 92949 and RESET 1.376 vs 1.392 — the
# same run.  The apparent y-null difference was a DIFFERENT RNG REALISATION,
# because adding 5 reactions changes the reaction-list length and hence the SSA
# draw sequence; with the measured seed-to-seed y-null spread (0.012 to 0.28)
# that fully accounts for it.
# ⚠ GENERAL RULE THIS COST US: before attributing any effect to a switch, check
# that the switch's reactions CAN FIRE in the network under test — every gate
# species needs a producer and a non-zero IC.  A switch whose gate is a species
# from an absent module is a no-op that still perturbs the RNG stream.
#
# WHAT IS ACTUALLY TRUE, measured in the FULL network (calib/m8_probe.jl single):
# with the arm OFF the cell does NOT reset after a pulse — PIP3 keeps CLIMBING
# once the stimulus ends (4 747 -> 619 376 -> 1 603 147 /cell at t = 120/150/295,
# stimulus off at t = 160) because r3_3b sustains the Ras/PIP3 loop with Gbg
# already falling.  The arm is REQUIRED for the post-pulse reset; it is off by
# default here only because the committed M3/M4 gradient numbers were taken
# without it, NOT because it is harmful.
# ⚠ THIS IS A TUNING-STAGE STATE, NOT A CLAIM THAT THE ARM IS WRONG.  The arm is
# what returns PIP3 to baseline after a pulse, and Module 8's relay and the
# refractory period were calibrated with Module 4 in a different regime — so
# re-enabling it is expected, but it re-opens this calibration and must be
# re-measured on the y-null, not just on amplitude.  DICTY_M4_REFRACTORY=1
# restores it.
# ⚠ BACK TO 1 [2026-08-29].  The arm is REQUIRED: with it off the full-network
# cell never resets after a pulse (PIP3 4 747 -> 1 603 147 /cell and still
# CLIMBING at t=295 with the stimulus off since t=160), which holds Module 8's
# PIP3 relay gate open and makes the refractory logic meaningless.
# ⚠ This does NOT disturb the committed M3/M4 gradient numbers: those were taken
# on calib/m4_polarity_probe.jl, whose reduced network cannot fire this arm at
# all (r4_12's gate TORC2a_mem is a Module 8 species with no producer there), so
# the arm is a no-op in that probe either way.
# ACCEPTANCE, on the engine (calib/m8_probe.jl, full network):
#   export 9.32e6 molecules/pulse   (target ~1e7)
#   blocked 420 s 2/2 and 480 s 4/4; fires again at 540 s 3/3  -> block = 8 min
#     (9 valid rows over 5 seeds, unanimous; rows where pulse 1 did not fire are
#      excluded - the fixed-onset trap, see calib/m8_probe.jl)
#   PTEN_mem reset 99 %;  PIP3 312x-runaway -> 3.1x basal at 8 min
# ⚠ STILL OPEN: RasG-GTP resets only 77-78 %, and that is NOT brake-limited -
# it is unchanged across a 6.7x authority range, a 10x pool range and every
# k_pip3Ras tried.  Gbg decays with tau = 55 s and is still 231x elevated at
# t=290, so the residual is Module 2's recapture, not Module 4's.
const M4_REFRACTORY = get(ENV, "DICTY_M4_REFRACTORY", "1") == "1"

const MODULE_4_REACTIONS = [r4_1, r4_2, r4_3, r4_4, r4_5,
                            M4_COOP ? r4_6c : r4_6,     # see the 4.6c block above
                            r4_7,
                            r4_8,       # PTEN-independent PIP3 removal — see the 4.8 block
                            (M4_REFRACTORY ? [r4_10, r4_11, r4_12, r4_13, r4_14] : Reaction[])...]

const PI3K_v  = parse(Int, get(ENV, "DICTY_M4_PI3KV",  "200000")) # 154 -> 1e4 (2026-08-05) -> 2e5/voxel (2026-08-19, see note above)
const PIP2_v  = parse(Int, get(ENV, "DICTY_M4_PIP2V",  "769231"))  # 769 -> 1e7 molec/cell total (was ~1e4/cell)
const PTEN_v  = parse(Int, get(ENV, "DICTY_M4_PTENV",  "10000"))  # 154 -> 1e4/voxel, 2026-08-05 (see note above)
#
_m4frac(pool, frac) = round(Int, frac * pool)

_uM(n)      = n / Float32(MOLEC_PER_µM)
_pi3k_f     = k_pi3kOn * _uM(RasGTP_v_rest) / (k_pi3kOn * _uM(RasGTP_v_rest) + k_pi3kOff)
const PI3Kmem_v_rest = max(round(Int, _pi3k_f * PI3K_v), 1)
_pten_f     = k_ptenOn * _uM(PIP2_v) / (k_ptenOn * _uM(PIP2_v) + k_ptenOff)
const PTENmem_v_rest = round(Int, _pten_f * PTEN_v)
_prod       = k_pip3 * _uM(PI3Kmem_v_rest)
# 4.8/4.9 are already pseudo-first-order (k_ship_eff/k_inpp4_eff, defined next
# to k_ship/k_inpp4 above), so the balance below uses them directly rather than
# re-deriving `k_ship * _uM(SHIP_v)` from the now-removed SHIP_v/INPP4_v pair.
_rem        = k_pten * _uM(PTENmem_v_rest) + k_ship_eff
const PIP3_v_rest   = max(round(Int, PIP2_v * _prod / (_prod + _rem)), 1)
# `PI34P2_mem` is an orphan now that 4.8 goes straight to PIP2 (see the 4.8
# block): no producer, so its rest level is 0 and the lipid budget below no
# longer subtracts it.
const PI34P2_v_rest = 0
const M4_INITIAL = Dict{Symbol, Int}(
    :PI3K_cyto  => PI3K_v - PI3Kmem_v_rest, :PI3K_mem => PI3Kmem_v_rest,
    :PIP3_mem   => PIP3_v_rest,
    :PI34P2_mem => PI34P2_v_rest,
    :PIP2_mem   => PIP2_v - PIP3_v_rest,
    :PTEN_cyto  => PTEN_v - PTENmem_v_rest, :PTEN_mem => PTENmem_v_rest,
    # Refractory arm at rest: PIP3 is at its basal 0.2 % so recruitment (4.10)
    # is negligible and the whole PKB pool sits cytosolic — the cell starts
    # NON-refractory, which is what lets the first stimulus fire.
    :PKBA_mem   => 0,  :PKBAa_mem => 0,  :PKBA_cyto => PKBA_v,
)

const M4_DIFFUSION = Dict{Symbol, Float32}(
    :PI3K_cyto  => 10.0f0,  :PI3K_mem   => 0.1f0,
    # PTEN_cyto is the SHUTTLE that makes winner-take-all possible, and its
    # cell-crossing time L²/(4D) = 2.5 s at D = 10 is NOT scaled by M4_SPEED —
    # so it becomes rate-limiting once the sped-up PTEN chemistry runs faster
    # than that. ENV-tunable for exactly that reason; 10-30 µm²/s is the generic
    # in-cell cytosolic range this file uses (CELL_DIFF).
    :PTEN_cyto  => parse(Float32, get(ENV, "DICTY_M4_DPTENCYTO", "10.0")),  :PTEN_mem => 0.1f0,
    # PKB: cytosolic when free, membrane-slow once PIP3-recruited. The ACTIVE
    # form must stay LOCAL or the refractory signal would terminate the whole
    # cell instead of the patch that produced it.
    :PKBA_cyto  => 10.0f0,  :PKBA_mem   => 0.1f0,  :PKBAa_mem => 0.1f0,
    :PIP2_mem   => 0.1f0,   :PIP3_mem   => 0.1f0,  :PI34P2_mem => 0.1f0,
)

# ============================================================================
# Module 3-4 coupling — PIP3/PI(3,4)P2 feedback onto Ras (LEGI-BEN)
# ============================================================================

# ── 3.3b IS A THREE-WAY PAIRING: feedback <-> M4 gate <-> M3 basal inhibitor ──
# [2026-08-27]  40 -> 55, and it CANNOT be moved alone.  The PIP3 -> RasG-GTP
# feedback is the largest available lever on the PIP3 gradient — larger than any
# increase in RasG-GTP production from the Module 2 signal — but each of its two
# partners fails FIRST, in a different place, and each failure looks like a
# different bug:
#   * `K_DISP` held fixed  -> the extra PIP3 pushes the 4.6c eviction gate toward
#     saturation and the gradient COLLAPSES.  Measured at kr = 70: F/B 10.12 with
#     K = 0.28, but 1.13 with K = 0.17.  This one shows up in the DRIVEN state.
#   * `k_gapBasal` held fixed -> the RESTING cell latches ON, because 3.3b now
#     out-drives the only Ras off-term at zero stimulus.  Measured at K_DISP
#     paired correctly: rest RasG-GTP 0.17 % at kr = 40, but 47.3 % at kr = 55
#     and 55.3 % at kr = 70.  This one is INVISIBLE to any driven scan.
# ⚠ AND `k_gapBasal` CANNOT SIMPLY BE RAISED TO MATCH, because it is
# stimulus-INDEPENDENT: it suppresses the response as hard as it suppresses rest.
# At kr = 70 every value that fixed the rest state also killed the drive (driven
# PIP3 0.9-1.7 % of the lipid pool, gate ~0.003, F/B 1.18-1.30).  The window is
# narrow and closes fast — kr = 62 is already HYSTERETIC (0.13 vs 9.59 from the
# two starts) while kr = 55 is unique (8.758 / 8.758).
# ⚠⚠ THE 55 / 0.22 / 2.5e-3 TRIPLE WAS TRIED AND **REVERTED** [2026-08-27].
# It is kept documented here because the reasoning above is sound and the triple
# looked like a clear win right up until the acceptance test was run properly.
# WHAT WENT WRONG, in order:
#  1. It was selected against a PROBE ARTEFACT.  calib/m4_polarity_probe.jl had
#     an asymmetric y-margin that put a real absorbing-boundary cAMP gradient on
#     the very axis it uses as its signal-free null (fixed, see Y_MARGIN there).
#     The low-gain network does not amplify that; a high-gain one does — so the
#     triple's own gain manufactured the evidence against it AND for it.
#  2. The first "PIP3 front/back 1.33 -> 2.08, target met" was ONE SEED. RETRACTED.
#  3. Clean geometry, FOUR seeds each, paired (this is the acceptance test):
#        new 55/0.22/2.5e-3 : 1.23  1.60  4.66  1.64   mean 2.28 +/- 0.80  median 1.62
#        old 40/0.17/1.5e-3 : 1.27  1.30  1.67  1.12   mean 1.34 +/- 0.12  median 1.29
#        paired diff +0.94 +/- 0.69, t = 1.36  -> NOT SIGNIFICANT
#     The whole mean is carried by seed 3; drop it and new averages 1.49 vs 1.23.
#     Seeds reaching the F/B >= 2 target: new 1/4, old 0/4.  BOTH MISS.
#  4. THE COST IS REPRODUCIBLE WHERE THE BENEFIT IS NOT: |y-null| 0.149 +/- 0.094
#     (new) vs 0.010 +/- 0.005 (old) — 15x worse, in every seed.  And the seed
#     that produced the impressive 4.66 is the SAME seed with the worst y-null
#     (-0.413): the big ratio and the off-axis domain are one event.  The network
#     fires a strong patch the gradient only loosely steers; it is not reporting
#     the gradient better.
# CONCLUSION: 40 is committed.  Raising 3.3b buys variance and off-axis
# polarisation, not gradient fidelity.  The gain needed for F/B >= 2 is not
# available inside Module 4 without going excitable — it has to come from
# upstream contrast (Module 2).
# Anything reading these notes should also know what DOES survive from that
# round and is still committed: r4_8 (PIP3's PTEN-independent sink, which had
# never existed), the DICTY_IC=auto solved rest state, k_gapBasal > 0 (required
# at ANY feedback gain), and M4_REFRACTORY defaulting off.
const k_pip3Ras    = _m3("DICTY_M3_KPIP3RAS",   40.0f0,  0.0f0) / Float32(PIP2_SCALE)  # µM⁻¹s⁻¹  (3.3b), catalyst-compensated
const k_pi34p2Gap  = _m3("DICTY_M3_KPI34P2GAP",  0.2f0,  0.0f0) / Float32(PIP2_SCALE)  # µM⁻¹s⁻¹  (3.4c), catalyst-compensated
const k_pip3GapOff = _m3("DICTY_M3_KPIP3GAPOFF", 1.0f0,  0.0f0) / Float32(PIP2_SCALE)  # µM⁻¹s⁻¹  (3.4d), catalyst-compensated

# 3.3b RasG_GDP_mem + PIP3_mem → RasG_GTP_mem + PIP3_mem
r3_3b = Reaction([:RasG_GDP_mem, :PIP3_mem], [:RasG_GTP_mem, :PIP3_mem], k_pip3Ras, nothing)

# ── 3.3b IS THE DEFAULT AGAIN [2026-10-02, later the same day]; 3.3g IS AN OPTION ──
# `DICTY_M3_PIP3PATH = ras (DEFAULT, 3.3b above) | gbg (3.3g below)`.
# WHY 3.3b.  The literature check (sources/M3/SYNTHESE_Rueckkopplung.md, sections/m3.tex)
# found every G-protein measurement against 3.3g: G-protein activation persists while
# downstream responses adapt (Janetopoulos 2001 p. 2408; Elzie 2009 p. 2598), rises
# step-like under two cAMP steps while Ras is transient (Xu 2022 p. 3), and the Gβγ
# domain change is PI3K-independent (van Hemert 2010 p. 2927).  The published Ras/PIP3
# model puts the PIP3 feedback into the Ras GEF term (Fukushima 2019 p. 10, Eqn 11:
# "One term defines the basal activity of Ras and the other defines feedback from
# PIP3"), which is 3.3b.  Its PIP3 share of Ras activity is small (Li 2018 p. E9128;
# Sasaki 2004 p. 511: LY reduces, does not abolish).  k_pip3Ras = 40 is the value that
# was committed for 3.3b before the switch to 3.3g.
# The history of the 3.3g attempt follows; it is kept as the `gbg` option.
# ── (history) 3.3g — THE PIP3 FEEDBACK ENTERS AT Gβγ, NOT AT RasG-GTP [2026-10-02, morning] ─
# WHY.  3.3b made PIP3 convert RasG-GDP to RasG-GTP DIRECTLY.  No source in the
# repo documents such a PIP3 → Ras-GEF edge; what is documented is (i) the
# Ras → PI3K → PIP3 → Ras loop as a whole (Sasaki 2004, Fukushima 2019) and
# (ii) that "the amplification sits between the G protein and Ras" (Kataria
# 2013).  The feedback is therefore moved UPSTREAM of Ras, onto the G protein:
# PIP3 accelerates the dissociation of the heterotrimer, which releases Gβγ,
# and Gβγ then drives everything it already drives (3.1 → RasG, 3.4e → GAP,
# 8.1 → RasGEFA).  PIP3 no longer touches RasG-GTP at all.
#   3.3g  Gabg_mem + PIP3_mem → Ga2GTP_mem + Gbg_cyto + PIP3_mem
# It is the same product set as 2.1 (receptor-catalysed exchange) with PIP3 as a
# second catalyst, and r2_5 (reassociation) removes the extra Gβγ again, so the
# G-protein pool is conserved.
# ⚠ EVIDENCE CHECK 2026-10-02 (sources in the repo): NO source measures a PIP3/PI3K effect
# on the G protein, and what exists points the other way.  Elzie 2009: membrane
# G-protein activation does not decline during continuous stimulation although PI3K
# activity subsides; van Hemert 2010: the cAMP-induced Gβγ mobility change/domain
# formation is PI3K-independent (mobility, NOT dissociation, was measured); Kortholt
# 2013: Ras symmetry breaking needs Gα2/Gβγ but not the PIP3 pathway; Fukushima 2019:
# Ras waves occur without PIP3 (PI3K inhibition only destabilises them).  The measured
# PIP3 → Ras effect has no documented attachment point.  `DICTY_M3_PIP3PATH=ras` is the
# alternative; decide by tuning which one restores the rest state.
# ⚠ WHAT THIS CHANGES FOR THE LEGI STAGE — read before tuning.  Gβγ drives BOTH
# arms (3.1 excitor, 3.4e inhibitor), so the feedback is amplified by the
# excitor and partly CANCELLED by the inhibitor (the GAP is armed by the same
# Gβγ), and it also raises RasGEFA (M8) and thereby the RasC arm.  This is
# deliberate: it is what "feedback onto Gβγ" means, and it is the reason the
# gain below cannot be taken over from 3.3b.  It also means the Gabg pool
# (≈70 % dissociated at the working dose, see Module 2) caps what the loop can
# add: the feedback can only act on the ≈30 % that is still intact.
# ⚠ k_pip3Gbg IS NOT CALIBRATED.  Starting value = the old k_pip3Ras per molecule
# (same units, same PIP2_SCALE catalyst compensation) so that the first run is
# comparable; there is no measurement for it.  TUNING: raise it until the cell
# fires and the PIP3 patch forms at the working dose, lower it until the rest
# state is the low branch (RasG-GTP ≈ 1, PIP3 ≈ 8 /voxel, see k_gapBasal).  The
# check list is in the 3.3c′ block (the brake is what restores the rest state).
const M3_PIP3PATH = get(ENV, "DICTY_M3_PIP3PATH", "ras")   # ras (3.3b, DEFAULT) | gbg (3.3g)
const k_pip3Gbg   = _m3("DICTY_M3_KPIP3GBG", 40.0f0, 0.0f0) / Float32(PIP2_SCALE)  # µM⁻¹s⁻¹ (3.3g), catalyst-compensated, UNTUNED
r3_3g = Reaction([:Gabg_mem, :PIP3_mem], [:Ga2GTP_mem, :Gbg_cyto, :PIP3_mem], k_pip3Gbg, nothing)
#r3_4c removed since only depended on constant species => useless
# 3.4d RasGAPa_cyto + PIP3_mem -> RasGAP_cyto + PIP3_mem
# ⚠ OFF BY DEFAULT (DICTY_M3_PIP3GAP=1 restores it): "biologically not
# validated" — and note it is NOT a brake despite the old MODULE_3_REACTIONS
# comment calling 3.3b/3.4d "STEN amplifier + brake".  PIP3 INACTIVATING the
# GAP removes RasG-GTP's inhibitor, so 3.4d is a SECOND positive-feedback arm
# stacked on 3.3b, not a negative one.
# It stays DEFINED rather than commented out so that switching it off cannot
# make the file unloadable — the exact failure this line caused (MODULE_3_
# REACTIONS still referenced it, so every include died with UndefVarError).
const M3_PIP3GAP = get(ENV, "DICTY_M3_PIP3GAP", "0") == "1"
r3_4d = Reaction([:RasGAPa_cyto, :PIP3_mem], [:RasGAP_cyto, :PIP3_mem], k_pip3GapOff, nothing)

# ── 3.3c/d/e — THE STEN BRAKE, WRITTEN 2026-09-06 ───────────────────────────
# ⚠ 2026-10-02: THE BRAKE IS A PROPOSED TERM WITHOUT LITERATURE SUPPORT FOR RasG,
# and the amplifier it brakes is now 3.3g (PIP3 → Gβγ), not 3.3b.  Read the
# "PROPOSED" block below (3.3c′) for the status and the tuning decision rule.
# The text below is the history of how it was built; where it says "3.3b" read
# "the PIP3 amplifier" (3.3g by default).
# ⚠ THIS IS A NEW REACTION EDGE, OFF BY DEFAULT (`DICTY_M3_BRAKE=1` enables).
#
# WHY IT DID NOT EXIST.  `MODULE_3_REACTIONS`' own comment calls its contents
# "STEN amplifier + brake" and a 2026-08-26 note in this file says "⚠ THE BRAKE
# IS NOT OPTIONAL … 3.3c-e is what makes the result EXCITABLE rather than
# LATCHED … They are restored together and DICTY_M3_STEN=0 removes both" — but
# the list contains ONLY `r3_3b`, and `r3_3c`/`r3_3d`/`r3_3e` were never
# written.  Everything else for them was: `RasBrake_cyto`/`RasBrakea_mem` are
# real species with `CELL_DIFF` entries (`D_M3_BRAKE`), a pool constant
# (`RasBrake_v` = 154/voxel) and an IC that this file's own header describes as
# "seeds an inert species at 0".  Same defect class as the `r4_8`/`r4_9` sink
# that "had never existed" — the prose, the species and the constants are all
# present and only the edges are missing.
#
# WHY IT MATTERS FOR GAIN.  3.3b is autocatalytic (PIP3 → Ras-GTP → PI3K →
# PIP3).  With no delayed inhibitor the loop has exactly two behaviours: too
# little gain to amplify anything, or bistable latching.  That is precisely what
# the 55/0.22/2.5e-3 attempt measured and why it was reverted — its own note
# records the cost as "variance and off-axis polarisation, not gradient
# fidelity", |y-null| 15x worse in every seed, and the impressive 4.66 front/back
# coming from the SAME seed as the worst y-null: "the network fires a strong
# patch the gradient only loosely steers".  A latched patch is what an amplifier
# without a brake produces.  The conclusion drawn there — that the gain "has to
# come from upstream contrast (Module 2)" — is now measurable as a dead end:
# Module 2 sits at 70 % dissociated at the working dose, so (1-f) caps its
# transfer near 0.30, and both routes to lowering that (less drive, bigger pool)
# break a literature anchor (the Janetopoulos 3.2 s dissociation time, and a
# G:R ratio already 24x literature).
#
# THE BIOLOGY.  Dictyostelium's STEN is an EXCITABLE network, not a bistable
# one: it fires, propagates and RESETS, and a gradient biases where and when it
# fires rather than holding it on (Fukushima 2019; Miao 2017 — both already
# cited in this file's Module 3 notes).  Excitability needs fast autocatalysis
# plus a SLOWER negative feedback; 3.3b is the former and this is the latter.
# The structure below is the standard one and mirrors Module 0.9's own
# chemical-FitzHugh-Nagumo pair, which this repo has already validated:
#   3.3c  PIP3 recruits/activates the brake         (slow, follows the patch)
#   3.3d  the brake decays back to cytosol          (sets the refractory time)
#   3.3e  the active brake hydrolyses RasG-GTP      (terminates the patch)
# `RasBrakea_mem` is left on the membrane (D = D_M3_BRAKE) so the brake is
# LOCAL to the patch it terminates — a globally-mixed brake would suppress every
# patch equally and act as a gain reduction rather than a terminator, the same
# distinction CLAUDE.md draws for `RasGAPa_cyto` at D = 10 ("cannot create a
# near/far asymmetry ... which makes 3.3b the spatial knob").
#
# ⚠ THE THREE CONSTANTS ARE NOT CALIBRATED.  They are set to give a brake that
# is ~5x slower than the amplifier and can hydrolyse the Ras pool on the ~10 s
# timescale of a PIP3 patch, which is the regime the structure requires — not
# fitted to any measurement.  Anything quoting a number from a run with
# `DICTY_M3_BRAKE=1` must say so.
# ⚠ DEFAULTS ON since 2026-09-06, PAIRED WITH `k_gapBasal` = 0. The brake is not
# an optional extra: it is what makes an input-independent basal inhibitor
# unnecessary, and that term was costing the module its adaptation. See the
# `k_gapBasal` block for the paired measurement (kinetics 1/6 → 5/6, step
# plateau 54.4 % → 3.5 %, peak/basal 1135× → 3.92×).
const M3_BRAKE      = get(ENV, "DICTY_M3_BRAKE", "1") == "1"
# ⚠ THE TWO GAINS ARE DERIVED FROM ANCHORS, NOT PICKED — because picking them
# produced a brake that could not brake. The first attempt used k_brakeOn = 0.02
# and k_brakeHydro = 2.0 "in the regime the structure needs", and MEASURED:
# `RasBrakea_mem` = 0.07 of a 154 pool (0.05 %) at the driven state, with
# brake-on and brake-off steady states identical to 5 significant figures
# (RasG_GTP 6390.9 vs 6391.0 at gain 200). Two independent reasons, both
# arithmetic, both invisible without doing the arithmetic:
#   * `k_brakeOn` is divided by PIP2_SCALE = 1000.3 (the catalyst compensation
#     every PIP3-driven rate in this block carries), so 0.02 is an EFFECTIVE
#     2.0e-5 µM⁻¹s⁻¹. Against k_brakeOff = 0.1 s⁻¹ at the driven PIP3 of
#     2.17 µM that gives an activated fraction of 4.3e-4 — the 0.07 measured.
#   * `k_brakeHydro` = 2.0 acting on 0.07 molecules is ~10⁴ times weaker than
#     the LEGI GAP arm it has to compete with (`k_rasGoff·[RasGAPa]` = 5.5 ×
#     1.93 µM = 10.6 s⁻¹ at the same point).
# So both are now SOLVED from a stated anchor, in this file's own "hold the
# measurement, derive the constant" style:
#   BRAKE_HALF_µM   PIP3 at which the brake is half-activated. Set to the DRIVEN
#                   PIP3 (2.17 µM at 60 nM ambient, measured on the rest-state
#                   solver) so the brake is off at rest — resting PIP3 is
#                   0.041 µM, giving ~1.9 % activation — and engages inside a
#                   patch. This is what makes it a terminator rather than a
#                   standing gain reduction.
#   BRAKE_STRENGTH  the brake's hydrolysis at half activation, as a FRACTION of
#                   the LEGI GAP arm's hydrolysis at the same point. 0.3 makes
#                   it a real but not dominant terminator; 0 reproduces a network
#                   with no brake at all.
# ⚠ NEITHER ANCHOR IS A MEASUREMENT FROM THE LITERATURE. They are design targets
# for an excitable-with-reset STEN, chosen so the mechanism can act at all.
# Any number quoted from a DICTY_M3_BRAKE=1 run must say so.
#
# ⚠⚠ MEASURED RESULT [2026-09-06]: THE BRAKE DOES NOT BUY GAIN HEADROOM FOR
# CONTRAST, BECAUSE GAIN IS NOT A CONTRAST LEVER. It works — 1.7 % activated at
# rest, ~47 % driven, and it cuts the resting latch tendency 2.5x (at
# k_pip3Ras = 200, resting Ras-GTP 7.26 % of pool → 2.86 %). It is also free on
# the transient, exactly as a 10 s delayed feedback against a 0.25 s PIP3
# response should be. But raising 3.3b behind it does NOT raise the contrast
# transfer d ln PIP3 / d ln c — it LOWERS it, monotonically:
#
#     PIP3 transient transfer   @2 nM    @20 nM   @60 nM
#       gain  40, no brake      0.245    0.101    0.060
#       gain  40, brake         0.259    0.107    0.063   ← best
#       gain 120, brake         0.232    0.100    0.059
#       gain 200, brake         0.205    0.094    0.056
#
# The reason is the one this file already records for `k_pip3` ("A PRODUCTION-
# RATE CONSTANT BUYS COUNTS, NOT CONTRAST"): 3.3b is positive feedback on the
# LEVEL, so it lifts the resting baseline as hard as the driven response, and a
# contrast is a ratio. This CONFIRMS the reverted 55/0.22/2.5e-3 round's
# conclusion — "Raising 3.3b buys variance and off-axis polarisation, not
# gradient fidelity" — and supplies the mechanism it was missing. Do not
# re-attempt a 3.3b gain raise as a contrast fix at any brake setting.
#
# WHAT IS STILL UNTESTED, and the only reason the brake stays in the file: the
# mean field cannot see the thing a brake is FOR. Its job is to terminate a
# patch so it neither latches nor drifts off-axis, and that is a spatial,
# stochastic property measured by the y-null of calib/m4_polarity_probe.jl —
# the very statistic that killed the 55/0.22/2.5e-3 triple (|y-null| 0.149 vs
# 0.010, 15x worse in every seed). The brake at the COMMITTED gain 40 is the
# configuration worth testing there; it costs nothing on contrast (+5 %) and
# might cost much less on the y-null.
# ⚠⚠ 2.17/0.3 → 1.5/0.6 [2026-09-15], TOGETHER WITH `RasBrake_v` 154 → 1540.
# The three move as ONE set; see the RasBrake_v block for the derivation and the
# full measurement. Short version: the M2 contrast retune (31e12a3) moved the
# driven PIP3 from 1.804 to 1.254 µM while the excursions did not follow, so the
# BRAKE_HALF anchor — defined as "the brake is half-activated at the DRIVEN
# PIP3" — pointed at an operating point the network had left, and the 3.3b loop
# ran under-damped (typical cell braked 37 % instead of 45 %).
# `DICTY_M3_BRAKEHALF=2.17 DICTY_M3_BRAKESTR=0.3 DICTY_M3_BRAKEV=154` restores it.
const BRAKE_HALF_µM  = parse(Float64, get(ENV, "DICTY_M3_BRAKEHALF", "1.5"))
const BRAKE_STRENGTH = parse(Float64, get(ENV, "DICTY_M3_BRAKESTR",  "0.6"))
const k_brakeOff    = _m3("DICTY_M3_KBRAKEOFF",  0.10f0, 0.0f0)   # s⁻¹  (3.3d) τ = 10 s
# half-activation at BRAKE_HALF_µM  ⇒  k_on·[PIP3]½ = k_off
const k_brakeOn     = Float32(Float64(k_brakeOff) / BRAKE_HALF_µM) # µM⁻¹s⁻¹ (3.3c), already effective
# at half activation the brake pool contributes RasBrake_v/2 molecules; ask its
# hydrolysis to be BRAKE_STRENGTH × the GAP arm's 10.6 s⁻¹ at the driven point
const k_brakeHydro  = Float32(BRAKE_STRENGTH * 10.6 /
                              ((RasBrake_v / 2) / Float64(MOLEC_PER_µM)))  # µM⁻¹s⁻¹ (3.3e)
# 3.3c PIP3_mem + RasBrake_cyto → RasBrakea_mem + PIP3_mem
r3_3c = Reaction([:PIP3_mem, :RasBrake_cyto], [:RasBrakea_mem, :PIP3_mem], k_brakeOn, nothing)
# 3.3d RasBrakea_mem → RasBrake_cyto
r3_3d = Reaction([:RasBrakea_mem], [:RasBrake_cyto], k_brakeOff, nothing)
# 3.3e RasBrakea_mem + RasG_GTP_mem → RasG_GDP_mem + RasBrakea_mem
r3_3e = Reaction([:RasBrakea_mem, :RasG_GTP_mem], [:RasG_GDP_mem, :RasBrakea_mem], k_brakeHydro, nothing)

# ── 3.3c/d/e — THE BRAKE IS A PROPOSAL, NOT A LITERATURE MECHANISM [2026-10-02] ─
# ⚠⚠ PROPOSED.  NOT TAKEN FROM THE LITERATURE.  KEPT ONLY IF TUNING SHOWS IT IS
# NEEDED TO RETURN THE SYSTEM TO ITS REST STATE.
# WHAT THE SOURCES DO AND DO NOT DOCUMENT (checked against the PDFs in sources/):
#   * RasG's documented brake is the Gβγ-driven RasGAP of the incoherent
#     feed-forward loop (Takeda 2012; the RasG GAP is NF1/NfaA, Zhang 2008,
#     Kortholt 2013) = reactions 3.4/3.4e/3.7 in this file.  That is ADAPTATION,
#     not a feedback from PIP3 or PKB.  Takeda 2012 tested and did not support a
#     GAP activated by Ras-GTP (integral control).
#   * PKB/PKBR1 → Sca1 → RasGEFA ↓ is documented for RasC (Charest 2010, Cai 2010,
#     Scavello 2017) = Module 8, reaction 8.2s.  Cai 2010: RasG-Q61L does NOT
#     change PKB phosphorylation, so the RasC pathway is separate from RasG.
#   * PKB → RasG is NOT documented.  Miao 2017 reports elevated RBD (an
#     isoform-NON-specific Ras sensor) and PHcrac in pkbA⁻/pkbR1⁻ cells and calls
#     the PKBs "candidates"; "the mechanism for the RasG reporter is not
#     identified".
#   * No source documents a PIP3-activated RasG brake either (see 3.3c above).
# So THIS EDGE CROSSES THE TWO BRANCHES: B (a Sca1-like substrate of PKBR1*, which
# is RasC's pathway) is made to hydrolyse RasG-GTP (3.3e).  It is a design
# choice that keeps the behaviour of the previous default (activation by PKBR1*,
# τ = 10 s decay, hydrolysis of RasG-GTP) and nothing more.
# PROPOSED ROLE.  Without any delayed inhibitor the PIP3 → Gβγ → RasG → PI3K →
# PIP3 loop (3.3g) has two regimes: too little gain, or a latched patch / latched
# rest state (measured, see k_gapBasal).  The brake is what removed the latched
# rest state (RasG-GTP 90 % of pool without it, ≈1 with it).
# HOW TO DECIDE WHETHER IT STAYS (do this when tuning 3.3g):
#   1. Run the rest-state solve + a 1800 s no-stimulus run with the brake OFF
#      (`DICTY_M3_BRAKE=0`).  If rest is already the low branch (RasG-GTP ≈ 1,
#      PIP3 ≈ 8 /voxel) and the cell fires at the Module-0 rate (2.41 spikes /
#      1800 s), the brake is NOT NEEDED: delete 3.3c/d/e and the 3.3c′ reader and
#      say so in the text.  The Gβγ route may already damp the loop, because the
#      GAP arm 3.4e is armed by the same Gβγ.
#   2. If it latches with the brake off, switch it on and tune ONLY
#      `DICTY_M3_BRAKEHALFPKB` (PKBR1* level, molecules/voxel, for half
#      activation; default 20), `DICTY_M3_BRAKESTR` (hydrolysis as a fraction of
#      the LEGI GAP arm at half activation; default 0.6) and
#      `DICTY_M3_KBRAKEOFF` (decay, default 0.1 s⁻¹ = τ 10 s).  The brake has to
#      be OFF at rest and after adaptation (Takeda 2012: the adapted level is set
#      by the IFFL) and act only on the transient.
#   3. Alternative reader if the PKB reading does not help: `DICTY_M3_BRAKEREAD=pip3`
#      (3.3c, activation by PIP3, the earlier default).  Report which reader was
#      needed; that is the result to carry back into the TeX (Module 3).
# Any number quoted from a run with the brake ON must say it is a proposed term.
#
# ── 3.3c′ — THE BRAKE'S READER, PKB-BASED (default since 2026-09-30) ────────────
# `DICTY_M3_BRAKEREAD = pkb (DEFAULT) | pip3 (the 3.3c above, previous default)`.
# NOTE ON THE EARLIER COMMENT HERE: it said Miao 2017 "extends it to the RasG
# reporter".  That overstated it — see the list above; Charest 2010 is RasC only
# ("RasG activation is unaffected" by the Sca1 complex).  Kept for the record:
#   * Charest et al. 2010 Dev Cell 18:737: PKB and PKBR1 phosphorylate Sca1
#     ("a peak at ~5–10 sec after stimulation"), releasing the Sca1/RasGEF
#     (Aimless)/PP2A complex from the membrane; in pkbA⁻/pkbr1⁻ cells RasC
#     activation "fails to rapidly adapt as it normally does by 40 s".  RasC.
#   * Weights below: Sca1 phosphorylation is "reduced in cells lacking PKB
#     (pkbA⁻), reduced to a greater extent in pkbr1⁻ cells, and abolished" in
#     both, i.e. both kinases, PKBR1 the larger (again a RasC statement).
# Constants are anchors, see calib/BRAKE_PKB_2026-09-29.md.
# ⚠ [2026-10-02] NO LONGER THE DEFAULT — superseded by the RasG-GTP reader 3.3c_R below.
# `DICTY_M3_BRAKEREAD = ras (DEFAULT) | pkb (this block) | pip3 (3.3c)`.
const M3_BRAKE_READ   = get(ENV, "DICTY_M3_BRAKEREAD", "ras")
M3_BRAKE_READ in ("ras", "pkb", "pip3") ||
    error("DICTY_M3_BRAKEREAD must be ras|pkb|pip3, got $(repr(M3_BRAKE_READ))")
# PKB activity (weighted molecules/voxel) at which the brake is half-activated
const BRAKE_HALF_PKB  = parse(Float64, get(ENV, "DICTY_M3_BRAKEHALFPKB", "20"))
const BRAKE_W_PKBR1   = parse(Float64, get(ENV, "DICTY_M3_BRAKEWR1", "1.0"))
const BRAKE_W_PKBA    = parse(Float64, get(ENV, "DICTY_M3_BRAKEWA",  "0.0"))
const k_brakeOnPKB    = Float32(Float64(k_brakeOff) * Float64(MOLEC_PER_µM) / BRAKE_HALF_PKB)  # µM⁻¹s⁻¹
r3_3cR = Reaction([:PKBR1a_mem, :RasBrake_cyto], [:RasBrakea_mem, :PKBR1a_mem],
                  Float32(BRAKE_W_PKBR1) * k_brakeOnPKB, nothing)
r3_3cA = Reaction([:PKBAa_mem, :RasBrake_cyto], [:RasBrakea_mem, :PKBAa_mem],
                  Float32(BRAKE_W_PKBA) * k_brakeOnPKB, nothing)

# ── 3.3c_R — THE BRAKE IS READ BY RasG-GTP (DEFAULT since 2026-10-02) ────────────
# `DICTY_M3_BRAKEREAD=ras`.  The brake becomes a negative feedback loop with a
# buffer node on RasG itself (NFBLB, Xu et al. 2022 p. 1: "the output is shut down
# by an inhibitor induced by the output itself").  Sources (all in sources/M3,
# marked in sources_marked/M3, see sections/m3.tex 3.3c_R):
#   * the reader can be neither PIP3 nor PKB(R1): Ras is excitable with the PIP3,
#     TorC2, PLA2 and sGC pathways all inhibited and without F-actin (Fukushima 2019
#     p. 2, p. 7), and the cAMP Ras response is unchanged in sgc/pla2(/pkbR1)-null+LY
#     cells, with LatA too (Kortholt 2011 p. 1273-1274);
#   * a RasGAP whose membrane recruitment and activation require Ras exists
#     (C2GAP1, Xu 2022 p. 2), and c2gapA- cells have elevated basal Ras
#     (Xu 2021 p. 4); a RasGAP activated by Ras-GTP reproduces the transient,
#     adaptive Ras response (Xu 2022 p. 4);
#   * k_brakeOff = 0.1 s⁻¹ (τ = 10 s) now has an anchor: the refractory period
#     recovers with t½ ≈ 7 s (Artemenko 2016 p. E7503), τ = 7 s / ln 2 ≈ 10 s.
# It is NOT Takeda's rejected integral controller (Supp. p. 3: GAP produced ∝ Ras-GTP,
# removed at a CONSTANT rate): B* decays in first order and has no setpoint, and the
# adapted level stays the IFFL's.  Takeda Supp. p. 9: the buffer must be slow to
# avoid oscillations — check the post-peak decay for sign changes when tuning.
# ⚠ BRAKE_HALF_RAS IS A DESIGN ANCHOR, NOT A MEASUREMENT: half activation at the
# driven RasG-GTP so the brake is ~off at rest and in the adapted plateau (its
# removal flux ∝ [B*][RasG-GTP] ∝ RasG-GTP² while unsaturated).  6400/voxel is the
# only driven RasG-GTP recorded in this file (6390.9 at gain 200, see the 3.3c block);
# re-measure the driven level at the committed gain and set it there.
# k_brakeHydro (3.3e) and RasBrake_v are unchanged.
const BRAKE_HALF_RAS = parse(Float64, get(ENV, "DICTY_M3_BRAKEHALFRAS", "6400"))  # RasG-GTP molec/voxel
const k_brakeOnRas   = Float32(Float64(k_brakeOff) * Float64(MOLEC_PER_µM) / BRAKE_HALF_RAS)  # µM⁻¹s⁻¹
# 3.3c_R RasG_GTP_mem + RasBrake_cyto → RasBrakea_mem + RasG_GTP_mem
r3_3cG = Reaction([:RasG_GTP_mem, :RasBrake_cyto], [:RasBrakea_mem, :RasG_GTP_mem], k_brakeOnRas, nothing)


# ============================================================================
# Module 8 — RasC → TORC2 → PKB → CRAC → ACA  (PIP3-independent relay cascade)
# ============================================================================
# 8.1  Gbg_cyto + RasGEFA_cyto → RasGEFA*_cyto + Gbg_cyto       (k_gefAon = 1e-3 nM⁻¹s⁻¹)
# ⚠ [2026-09-29] 8.1/8.2 ARE A COPY OF MODULE 3's LEGACY GEF (3.1/3.2 were 1.0 µM⁻¹s⁻¹ /
# 0.1 s⁻¹ — `_m3` legacy values) AND MISSED EVERY RE-PAIRING 3.1 GOT (3.0 → 0.038 → 0.0038,
# the last paired with G_RATIO 6 → 60, 2026-08-24). Measured consequence
# (calib/BRAKE_PKB_2026-09-29.md §4): resting Gβγ 2282/voxel puts RasGEFA at ~83 % active
# AT REST and ~99 % at 60 nM — a saturated excitor — while the RC-driven RasCGAP arm (8.4b)
# still rises, so a cAMP step drives RasC-GTP DOWN 12× (107 → 9/voxel) and TORC2/PKBR1 with
# it. Literature: RasC is rapidly and transiently GTP-loaded by cAMP (Kae et al. 2004) and
# drives TORC2 → PKB (Cai et al. 2010; Charest et al. 2010).
# `DICTY_M8_GEFAPAIR=1`: pair 8.1 with 3.1 (same activated-fraction curve in Gβγ, since both
# are Gβγ-driven RasGEFs — Kae et al. 2007), and size the pool so the DRIVEN activated count
# at the working point (Gβγ = M8_GBG_DRIVEN_V, measured at a 60 nM step) is what the old
# saturated arm delivered there — i.e. the relay drive at the stimulated state is held and
# only the resting level (hence the contrast) moves.
# DEFAULT ON since 2026-09-30 (`DICTY_M8_GEFAPAIR=0` restores the saturated legacy arm).
const M8_GEFA_PAIR    = get(ENV, "DICTY_M8_GEFAPAIR", "1") == "1"
const k_gefAoff       = 0.1f0
const k_gefAon_legacy = nM_to_µM(1f-3)
const k_gefAon        = M8_GEFA_PAIR ? Float32(Float64(k_gefRon) * Float64(k_gefAoff) / Float64(k_gefRoff)) : k_gefAon_legacy
const M8_GBG_DRIVEN_V = parse(Float64, get(ENV, "DICTY_M8_GBGDRIVEN", "67000"))  # molec/voxel, MEASURED (60 nM step)
_fA(k, g) = k * g / (k * g + Float64(k_gefAoff))
const GEFA_v_legacy   = 154
const GEFA_v_pair     = round(Int, GEFA_v_legacy * _fA(Float64(k_gefAon_legacy), M8_GBG_DRIVEN_V / MOLEC_PER_µM) /
                                   _fA(Float64(k_gefAon), M8_GBG_DRIVEN_V / MOLEC_PER_µM))
r8_1  = Reaction([:Gbg_cyto, :RasGEFA_cyto], [:RasGEFAa_cyto, :Gbg_cyto], k_gefAon,   nothing)
# 8.2s  RasBrakea_mem + RasGEFAa_cyto → RasGEFA_cyto + RasBrakea_mem   [2026-09-29, opt-in]
# THE Sca1 FEEDBACK AS CHAREST ET AL. 2010 MEASURED IT: PKB/PKBR1 phosphorylate the scaffold
# Sca1 and release the Sca1/Aimless(RasGEFA)/PP2A complex from the membrane, "and thereby
# RasC activity, in a negative feedback fashion"; without PKB or TORC2, RasC "fails to rapidly
# adapt as it normally does by 40 s" and its basal level is elevated. `RasBrakea_mem` is the
# PKB-phosphorylated Sca1-like substrate (3.3c′, DICTY_M3_BRAKEREAD=pkb), so the same species
# releases active Aimless here. `DICTY_M8_SCA1FB` = k in µM⁻¹s⁻¹ (DEFAULT 1 since
# 2026-09-30; 0 = off). ⚠ k = 1 vs 3 not yet decided (calib/BRAKE_PKB_2026-09-29.md round 3).
# [CALIBRATED] against 'RasC back near baseline by ~40 s' (Charest 2010).
const k_sca1Fb = parse(Float32, get(ENV, "DICTY_M8_SCA1FB", "1"))
# ⚠ [2026-10-02] Sca1 IS NOW ITS OWN SPECIES (`Sca1_cyto` / `Sca1a_mem`), NO LONGER
# Module 3's `RasBrakea_mem`.  The shared species made one substrate brake both RasC
# (documented, Charest 2010) and RasG (not documented: "the Sca1 complex regulates
# RasC activity while RasG activation is unaffected", Charest 2010 p. 740).  The RasG
# brake is now read by RasG-GTP (3.3c_R).  8.2p/8.2q copy the former PKB reader 3.3c′
# and decay 3.3d exactly (same k_brakeOnPKB, weights, k_brakeOff, pool, D), and every
# step is catalytic for Sca1*, so the RasC arm behaves as under the previous default
# (DICTY_M3_BRAKEREAD=pkb), only no longer gated on the M3 brake switches.
const Sca1_v = RasBrake_v                                              # Sca1-like substrate / voxel
# 8.2p PKBR1a_mem + Sca1_cyto → Sca1a_mem + PKBR1a_mem   (PKBR1 phosphorylates Sca1)
r8_2p  = Reaction([:PKBR1a_mem, :Sca1_cyto], [:Sca1a_mem, :PKBR1a_mem],
                  Float32(BRAKE_W_PKBR1) * k_brakeOnPKB, nothing)
# 8.2pA PKBAa_mem + Sca1_cyto → Sca1a_mem + PKBAa_mem   (PKBA, weight 0 by default)
r8_2pA = Reaction([:PKBAa_mem, :Sca1_cyto], [:Sca1a_mem, :PKBAa_mem],
                  Float32(BRAKE_W_PKBA) * k_brakeOnPKB, nothing)
# 8.2q Sca1a_mem → Sca1_cyto   (τ = 10 s, as 3.3d)
r8_2q  = Reaction([:Sca1a_mem], [:Sca1_cyto], k_brakeOff, nothing)
# 8.2s Sca1a_mem + RasGEFAa_cyto → RasGEFA_cyto + Sca1a_mem
r8_2s = Reaction([:Sca1a_mem, :RasGEFAa_cyto], [:RasGEFA_cyto, :Sca1a_mem], k_sca1Fb, nothing)
# 8.2  RasGEFA*_cyto → RasGEFA_cyto                             (k_gefAoff = 0.1  s⁻¹)
r8_2  = Reaction([:RasGEFAa_cyto], [:RasGEFA_cyto],    k_gefAoff,         nothing)
r8_3  = Reaction([:RasGEFAa_cyto, :RasC_GDP_mem], [:RasC_GTP_mem, :RasGEFAa_cyto], 2.0f0, nothing)
# 8.4  RasC_GTP_mem → RasC_GDP_mem  (constitutive, slow)       (k_RasC_off = 0.02 s⁻¹)
r8_4  = Reaction([:RasC_GTP_mem], [:RasC_GDP_mem],     0.02f0,            nothing)
const RasCGAP_v = 154        # RasCGAP pool molecules/voxel (mirrors RasGAP)
# 8.4b RC_mem + RasCGAP_cyto → RC_mem + RasCGAP*_cyto  (delayed receptor-driven activation)
r8_4b = Reaction([:RC_mem, :RasCGAP_cyto], [:RC_mem, :RasCGAPa_cyto], 0.05f0,  nothing)
# 8.4c RasCGAP*_cyto → RasCGAP_cyto  (slow, persistent inhibitor — holds RasC adapted)
r8_4c = Reaction([:RasCGAPa_cyto], [:RasCGAP_cyto],    0.005f0,           nothing)
# 8.4d RasCGAP*_cyto + RasC_GTP_mem → RasC_GDP_mem + RasCGAP*_cyto  (GAP-driven hydrolysis)
r8_4d = Reaction([:RasCGAPa_cyto, :RasC_GTP_mem], [:RasC_GDP_mem, :RasCGAPa_cyto], 150.0f0, nothing)

# ── 8.5 / 8.6 — RasC → TORC2.  [TUNE M8 2026-07-21] k_torOn 10 → 50 µM⁻¹s⁻¹,
# k_torOff 0.1 → 0.5 s⁻¹ (τ 10 s → 2 s).  BOTH scaled by the same factor 5.
#
# Why jointly: TORC2* is a first-order relaxation driven by RasC-GTP, so
#     τ_stage = 1/(k_torOn·[RasC]/M_VOX + k_torOff)   ← sets the LATENCY
#     f_ss    = k_torOn·[RasC]/M_VOX / (…)            ← sets the AMPLITUDE
# Scaling k_torOff alone (the obvious "make it faster" move, tried in
# validation/v6_sensitivity.jl) divides the amplitude by the same factor it
# divides the latency — that is what collapsed the ACA peak 53.6 % → 23.1 % and
# produced the "speed vs. drive" tension of TUNING.md iters 2–5.  Scaling the
# PAIR leaves f_ss invariant, so the stage gets 5× faster at NO cost in drive.
#
# Effect (validation/data/m8_calibration.md, sustained 1000 nM):
#   TORC2* peak 14.52 s → 7.46 s   (into the 5–10 s band; Liu 2017, Kamimura 2008)
#   PKBR1a peak 16.30 s →  9.30 s   (into the 5–10 s band — crit 8b.x #1 / 8.x #1)
#   ACA peak    53.6 %  →  55.2 %   of pool (slightly BETTER, not worse)
#   ACA decline t½ 113.4 s → 107.0 s (Cai 2010: ~108 s — fit improves)
#   ACA peak time 31.6 s → 26.1 s   (Cai 2010: 30 s — the one mild regression)
# This is the fix for discrepancy D1: the binding constraint on PKB timing was
# always Module 8's TORC2 stage, never Module 8b's PKBR1 rates.
const k_torOn  = parse(Float32, get(ENV, "DICTY_M8_KTORON", "50.0"))   # µM⁻¹s⁻¹  RasC-GTP-driven TORC2 activation (0 = PP242, TORC2 inhibited)
const k_torOff = 0.5f0    # s⁻¹      TORC2 deactivation, τ = 2 s
r8_5  = Reaction([:RasC_GTP_mem, :TORC2_mem], [:TORC2a_mem, :RasC_GTP_mem], k_torOn, nothing)
# 8.6  TORC2*_mem → TORC2_mem
r8_6  = Reaction([:TORC2a_mem], [:TORC2_mem],          k_torOff,          nothing)
# 8.7  TORC2*_mem + PKBA_cyto → TORC2*_mem + PKB*_mem           (k_pkbOn  = 10 µM⁻¹s⁻¹)
r8_7  = Reaction([:TORC2a_mem, :PKBA_cyto], [:PKBAa_mem, :TORC2a_mem],    10.0f0,            nothing)  # TUNED iter3 (see 8.5)
# 8.8  PKB*_mem → PKBA_cyto                                     (k_pkbOff = 0.1  s⁻¹, τ = 10 s)
r8_8  = Reaction([:PKBAa_mem], [:PKBA_cyto],           0.1f0,             nothing)
# ── CRAC is a THREE-STATE moiety  [RESTRUCTURED 2026-07-22, iteration 13;
#                                   chain CLOSED + retuned 2026-07-26, iteration 18] ───
#     CRAC_cyto  --(8.9b PIP3 | 8.9c basal)-->  CRAC_mem  --(8.9 PKBa | 8b.3 PKBR1a)-->  CRACa_mem
#                <--------- 8.10b -----------            <--------- 8.10 ----------------
# Conserved: CRAC_cyto + CRAC_mem + CRACa_mem = 154/voxel.
#
# WHY three states.  Two experimental facts are individually solid and jointly
# impossible under a single CRAC*_mem state:
#   (i)  CRAC is a PH-domain protein recruited to the membrane by binding PIP3 —
#        PHcrac-GFP is THE Dd PIP3 biosensor (Insall 1994; Parent 1998);
#   (ii) pkbR1⁻/pkbA⁻ double mutants ABOLISH ACA-relay cAMP synthesis and cannot
#        aggregate (Kamimura 2008; Cai 2010).
# A first attempt wrote (i) as `PIP3 + CRAC_cyto → CRAC*_mem` — fully active CRAC
# straight from PIP3.  That made PKB dispensable and the double null retained 24 %
# of WT ACA, i.e. it broke (ii) (battery 14/19 → 12/19; TUNING.md §13.6).  The two
# are compatible only if RECRUITMENT and ACTIVATION are separate steps: PIP3 puts
# CRAC on the membrane, PKB phosphorylates it there.  That is this structure, and
# it is more faithful than either the pre-13 model (no PIP3 dependence at all) or
# the rejected binary edge (no PKB dependence).
#
# WHY 8.9c, the PIP3-INDEPENDENT recruitment route, is REQUIRED and not padding.
# Module 8 is a **PIP3-independent** relay — that is the module's defining claim
# (DICTY_CELL.md §8 role; pi3k1/2-null cells still relay and aggregate, Funamoto
# 2002), and criteria 8.x#2 / 8b.x#3 test exactly that.  If PIP3 were the ONLY way
# onto the membrane, a PI3K-null would abolish the relay and we would have traded
# criterion 8b.3 for criterion 8b.x#3.  8.9c is the residual membrane association
# that carries the relay when PIP3 is absent — physically the PI(3,4)P2 route
# (§4.8: "CRAC binds PI(3,4)P2 too") plus non-specific membrane affinity.
#
# WHY 8.10 RETURNS CRAC* TO THE MEMBRANE, NOT TO THE CYTOSOL  [FIXED 2026-07-26].
# As first written, the three-state chain was closed by sending CRAC* straight back
# to CRAC_cyto — deactivation ALSO ejected the protein from the membrane.  That
# short-circuits the chain: every activation cycle has to pay the cyto→mem
# recruitment step again, so the CRAC* the cascade can reach is capped at
#     CRAC*/pool  ≤  k_recruit / (k_recruit + k_offA)  =  0.24 / 0.34  =  71 %
# no matter how hard PKB drives it, and the measured peak sat at 32 % of pool (ACA*
# 44.3 %).  It is also not what the protein does: CRAC is held at the membrane by its
# PH domain binding PIP3/PI(3,4)P2 (Parent 1998), and a phosphatase removing the PKB
# phosphate does not strip that lipid anchor.  Closing the chain step-by-step
# (CRAC* --8.10--> CRAC_mem --8.10b--> CRAC_cyto) removes the artificial ceiling,
# lets the activation cycle run entirely on the membrane, and puts the Parent-1998
# translocation constraint on 8.10b — the reaction that actually models release —
# where its "returns within seconds–tens of seconds" τ = 10 s is untouched.
# Measured effect of the topology fix ALONE (mean-field, 1000 nM clamp):
#   CRAC* peak 48.7 → 55.2, ACA* peak 44.3 % → 46.4 % of pool, all timings unmoved.
#
# RATES.  8.9/8b.3 keep their committed k_cracOn = 10 µM⁻¹s⁻¹; 8.9b keeps the
# canonical PIP3-recruitment constant 1e-3 nM⁻¹s⁻¹ (as 5.1/4.1/4.4); 8.10b keeps
# k_cracOff = 0.1 s⁻¹ (Parent 1998).  k_cracBas is DERIVED, not fitted: inserting an
# intermediate state divides the CRAC→CRACa flux by the membrane-pool occupancy, so to
# leave the committed relay strength intact the recruitment rate must balance the
# losses out of CRAC_mem,  k_cracBas ≈ k_cracOff + k_cracOn·(PKBa+PKBR1a)/M_VOX
# ≈ 0.1 + 0.17 ≈ 0.27 at the operating point; 0.2 is that value less the PIP3 term it
# now shares the job with.
#
# NB the PIP3 arm (8.9b) is numerically minor as parameterised: PIP3_mem rests at 0
# and peaks at 63 molec/vox under a saturating clamp, so 8.9b contributes
# 1.0·63/M_VOX ≈ 0.013 s⁻¹ against 8.9c's 0.2 s⁻¹ — ~6 % of recruitment at the peak
# and none at rest.  The edge is real and correctly signed (it is what makes the
# recruitment step stimulus-dependent at all), but it is NOT what carries the relay;
# 8.9c is.  Raising it is not free — the rate is pinned by the "same constant as
# 5.1/4.1/4.4" argument in DICTY_CELL.md §8.9b.
const k_cracOn   = 10.0f0        # µM⁻¹s⁻¹  PKB-driven activation of the membrane pool
const k_cracPIP3 = nM_to_µM(1f-3)  # µM⁻¹s⁻¹  PIP3-facilitated recruitment (8.9b)
const k_cracBas  = 0.2f0         # s⁻¹      PIP3-independent membrane association (8.9c)
const k_cracOff  = 0.1f0         # s⁻¹      CRAC_mem → CRAC_cyto release (8.10b; Parent 1998)
# [TUNED 2026-07-26, iteration 18]  8.10 is a DEPHOSPHORYLATION rate now that it no
# longer doubles as the membrane-release step, so it is no longer pinned by Parent
# 1998's translocation imaging (that constraint moved to k_cracOff/8.10b above).
# 0.1 → 0.03 s⁻¹ (τ 10 → 33 s).  Tuned against the ACA observables of criterion 8.x#5
# (Cai 2010: ACA* peaks ~30 s and declines with t½ ≈ 1.8 min = 108 s), which the model
# had been missing by a factor of 2.6 on the decline.  Sweep at the committed
# k_acaOn = 10 (mean-field, 1000 nM clamp), ACA* peak / t_peak / decline t½:
#   0.100 → 46.4 % / 23.2 s /  42 s      0.040 → 52.5 % / 24.1 s /  71 s
#   0.070 → 49.2 % / 23.6 s /  48 s      0.035 → 53.1 % / 24.2 s /  83 s
#   0.050 → 51.4 % / 23.9 s /  59 s      0.030 → 53.8 % / 24.4 s / 102 s  ← committed
# 0.03 is the point that lands the decline on Cai's 108 s.  Going lower buys almost
# no further t_peak (it asymptotes at ≈24.5 s — see the 8.x#5 note in DICTY_CELL.md)
# and overshoots the decline.  Population-validated: latchfrac 0.228 < 0.30 (PASS),
# 98 % of cells fire ≥2×, period 163 s (validation/data/population_m8_diffusion.csv,
# tag `crac_3state_deph003`; the pre-change control is `crac_ctrl_2state`, 0.194).
# ── [RE-TUNED 2026-08-09] k_cracDeph 0.03 → 10 s⁻¹.  CRAC* WAS SATURATED. ───
# MEASURED on the engine with Module 8 live (calib/m8_probe.jl, first run):
# CRAC* sits at 1923 of the 2002-molecule pool — 96 % — in a cell at REST, and
# reaches 100 % under stimulus.  A stage pinned at its ceiling carries no signal
# (the same defect Module 2's G-protein and Module 3's GAP arm had; CLAUDE.md
# "check the activated fraction of every switch before tuning its rates"), so
# with 8.11 gone Module 8 would have been a constant multiplier and the PIP3
# gate would have been doing all the work by itself.
#
# WHY IT WAS PINNED, and why that is new.  8.9 is driven by `PKBAa_mem`, which
# Module 4's refractory arm now supplies on a PKBA_v = 10^4/voxel pool: PKBAa_mem
# rests at 455/voxel and peaks at 5720, so the ON-rate k_cracOn·[PKBa] runs
# 0.94 → 11.9 s⁻¹ against an OFF-rate of 0.03.  Both ends of that swing are
# ≫ k_cracDeph, so the ratio saturates.  0.03 was calibrated in 2026-07-26
# against the ACA* decline (Cai 2010, t½ ≈ 108 s) at a time when CRAC* DROVE ACA
# through 8.11 and PKBa came from Module 8's own 77/voxel pool.  Neither holds:
# 8.11 is superseded by 8.14, ACA is driven by X through 0.9f and its decline is
# set by k_acaOff/9.3, and PKB is Module 4's.  The constraint that pinned this
# rate no longer exists, so it is free to be set for SIGNAL CONTRAST instead.
#
# CHOICE.  With k_cracDeph = d the activated fraction of the closed three-state
# chain (0.2 on / 0.1 off recruitment, a = k_cracOn·[PKBa]) is
#     f_a = a / (a + d·(1.5 + a/d))  →  rest 5.9 %, peak 44 % at d = 10
# i.e. a 7.5× swing out of a 12.6× drive swing, against 1.04× at d = 0.03.
# Larger d buys a little more contrast at a proportionally lower peak (d = 30:
# 2.1 %/21 %); 10 keeps the peak high enough that 8.14's gain stays modest.
# Physically this makes the CRAC phospho-cycle a fast covalent-modification
# cycle that FOLLOWS PKB (τ = 0.1 s) rather than integrating it — which is what
# it should be, since the integration in this pathway is done by Module 3's LEGI
# and Module 4's STEN upstream, not by the kinase substrate.
# ── ⚠ [RE-ANCHORED 2026-09-16] 10 → 1.5 s⁻¹.  THE RELAY HAD STOPPED WORKING. ──
# 10 s⁻¹ was not a free choice: the block above picks it from ONE measurement,
# `a_peak = k_cracOn·[PKBAa_mem]_peak = 11.9 s⁻¹` (PKBAa_mem peaking at
# 5720/voxel), because what it sets is the ACTIVATED FRACTION of the CRAC
# phospho-cycle, and the fraction is a function of `a/d` alone.  `a` has since
# moved by an order of magnitude and `d` did not follow, so the switch fell off
# its operating point — which is exactly the failure this file keeps recording
# under "check the activated fraction of every switch before tuning its rates",
# arriving this time through a constant that was right when it was written.
#
# MEASURED on HEAD, `calib/m8_probe.jl single` (one cell, field clamped, 0.06 µM
# box at t = 120 s):  PKBAa_mem peaks at 7373/CELL = **567/voxel**, i.e. 10.1×
# below the value 10 s⁻¹ was derived from, so a_peak = 1.18 s⁻¹ and
#     f(CRAC*) = (a/d)/(k_cracOff/r + 1 + a/d)  →  **11 % of pool**, not 44 %.
# The relay drive is `XGAIN·f(CRAC*)·gate(PIP3)`, so an 11 % CRAC* is a 4×
# under-drive, and the engine says that is the difference between firing and
# not: at 10 s⁻¹ the probe reads
#     pulse 1 @ 120 s → NO FIRE   X_pk 18   cAMP_i_pk 0   exported 0
# against a Module 0.9 excitation threshold of ≈24 molecules of X.  X_pk 18 is
# the whole story — the relay was pushing the cell to 75 % of threshold and
# stopping there.  A population of such cells cannot propagate a wave at all:
# every pulse it shows is a SPONTANEOUS one, and what looks like synchrony is
# the refractory clock releasing everybody at once.
#
# WHY `k_cracOn` OR `XGAIN` IS THE WRONG KNOB.  Both would restore the product
# and leave the fraction at 11 %, i.e. leave a switch sitting at a tenth of its
# range, where its contrast (and hence the relay's dose-independence) is gone.
# `d` is the constant that sets the fraction and the constant that was derived
# from `a`; re-deriving it is a restoration, not a new fit.  Inverting
# f = (a/d)/(1 + a/d) at the block's own target f = 0.44 gives
#     d = a_peak·(1 − f)/f = 1.18 · 0.56/0.44 = 1.50 s⁻¹.
#
# MEASURED AFTER, same probe, two independent seeds and two doses two decades
# apart — and these are CALIBRATION_M8.md's OWN published numbers, which a free
# re-tune would have no reason to land on:
#     0.002 µM  FIRED  X_pk 447  cAMP_i 386390  ACA* 94 % of pool  exported 8.9e6
#     0.060 µM  FIRED  X_pk 468  cAMP_i 380992  ACA* 94 % of pool  exported 8.6e6
#     CRAC*_pk 44 % of pool (target 44 %), Refr peak ≈ 6100
# τ = 1/d = 0.67 s still makes the phospho-cycle a FOLLOWER of PKB rather than
# an integrator, which is the block above's physical requirement.
#
# ── ✓ [2026-09-17] CONFIRMED IN A POPULATION, AND END-TO-END AT TWO CELLS ───
# `two_cells.jl` at 60 um, 3 seeds, real diffusive field, no clamp and no
# commanded pulse: relay delay 9-14 s (mean 12.3) => wave speed 4.9 um/s,
# against 18.3 +/- 0.5 s / 3.28 um/s before the repair and a literature ~5 um/s.
# 0/6 co-ignitions and the leader ALTERNATES, i.e. mutual entrainment rather
# than one cell driving a passive partner. The speed gain is mechanistic, not a
# fit: at 44 % CRAC* the gate opens on the 2 nM leading EDGE of the arriving
# pulse instead of needing its 60 nM peak, so the neighbour fires earlier.
# In a population (1000 cells / 400^2 / 3000 s, seeds 42/43/44),
# `analysis/kymograph.jl` classifies 13 of 21 firing episodes as carrying
# propagating fronts, against 0 of 2 -- pure FLASH -- on the pre-repair network.
# ⚠ THAT IS FRONTS, NOT ONE WAVE SWEEPING THE FIELD. The whole 784 um scan line
# lights up in 24-72 s, which as a single front would imply 10-33 um/s -- faster
# than this relay can propagate. Each episode is SEVERAL sources igniting
# near-simultaneously, each driving a short local front of ~4 um/s (the
# per-front regressions and two_cells agree on that figure). At rho = 6.25e-3
# with a 726 s spontaneous interval ~1.4 cells ignite per second field-wide, so
# the pacemaker spacing is ~180 um against an 800 um field. That is geometry and
# rate, NOT a constant to retune here -- and the 12-minute specification fixes
# the rate.
#
# ⚠ `M8_PKBA_PEAK_V` IS A MEASUREMENT AND GOES STALE WHEN MODULE 4 MOVES.  It is
# written out as a named constant, and `k_cracDeph` derived from it, precisely
# so that the next time Module 4's PKBA arm is re-anchored the repair is one
# number taken off `calib/m8_probe.jl single` rather than an archaeology
# exercise.  Re-measure it after ANY change to 4.10-4.14, `PKBA_v`, or the
# `DICTY_M4_PKBANCHOR` anchors.  `DICTY_M8_KCRACDEPH=10` restores the old value.
const M8_PKBA_PEAK_V = parse(Float32, get(ENV, "DICTY_M8_PKBAPEAK", "567.0"))  # molec/voxel, MEASURED
const M8_CRAC_FPEAK  = parse(Float32, get(ENV, "DICTY_M8_CRACFPEAK", "0.44"))  # target activated fraction
const k_cracDeph = parse(Float32, get(ENV, "DICTY_M8_KCRACDEPH",
    string(k_cracOn * M8_PKBA_PEAK_V / Float32(MOLEC_PER_µM) *
           (1f0 - M8_CRAC_FPEAK) / M8_CRAC_FPEAK)))                # s⁻¹  (8.10)  [10 -> 1.50, 2026-09-16]
# 8.9  PKB*_mem + CRAC_mem → PKB*_mem + CRAC*_mem              (k_cracOn = 10 µM⁻¹s⁻¹)
#   PKBa is itself membrane-localised, so it acts on the membrane-recruited pool.
r8_9  = Reaction([:PKBAa_mem, :CRAC_mem], [:CRACa_mem, :PKBAa_mem],       k_cracOn,         nothing)  # TUNED iter3 (see 8.5)
# 8.9b PIP3_mem + CRAC_cyto → PIP3_mem + CRAC_mem   (PH-domain recruitment, Parent 1998)
r8_9b = Reaction([:PIP3_mem, :CRAC_cyto], [:CRAC_mem, :PIP3_mem],         k_cracPIP3,       nothing)
# 8.9c CRAC_cyto → CRAC_mem   (PIP3-independent residual association; keeps M8 PIP3-independent)
r8_9c = Reaction([:CRAC_cyto], [:CRAC_mem],                               k_cracBas,        nothing)
# 8.10b CRAC_mem → CRAC_cyto  (release of the un-phosphorylated membrane pool)
r8_10b = Reaction([:CRAC_mem], [:CRAC_cyto],                              k_cracOff,        nothing)
# 8.10 CRAC*_mem → CRAC_mem   (dephosphorylation; the protein stays lipid-anchored)
# [2026-07-26] product CRAC_cyto → CRAC_mem, rate 0.1 → k_cracDeph = 0.03 s⁻¹.  See the
# "WHY 8.10 RETURNS CRAC* TO THE MEMBRANE" block above.  The old note here warned that
# re-slowing this rate stretched the ACA decline past 300 s — that was measured with the
# cytosol-returning product, where slowing 8.10 also starved the recruitment cycle; with
# the chain closed on the membrane, 0.05 s⁻¹ lands the decline at t½ = 102 s (Cai 2010:
# ~108 s), i.e. the constraint the old note was protecting is now MET, not violated.
r8_10 = Reaction([:CRACa_mem], [:CRAC_mem],            k_cracDeph,        nothing)
# [RE-ANCHORED 2026-07-22] KTH_RC 470 → 180.  This is a change of *observable*, not
# a fit.  The old 470 ⇔ 2.5 nM (the annotation is numerically right: RC_mem measures
# 445.8 at a sustained 2.5 nM, TUNING.md §11.5 §E) was taken from Gross et al. 1976's
# POPULATION RELAY THRESHOLD — the field at which a lawn propagates a wave.  8.11 is a
# SINGLE-CELL firing gate, and the population threshold necessarily sits above the
# concentration at which one cell half-activates its own cyclase, so the gate was
# anchored ~3.4× too high and could not open at physiological occupancy: measured
# spontaneous firing was 10^-150 per cell·s, ~147 decades below the ~10^-2.9 per
# cell·s (≈1 pulse/12 min) seen in pre-wave cells (TUNING.md §11.6).
# 180 ⇔ 0.73 nM is the single-cell half-activation constant:
#   • it reproduces the measured single-cell firing rate (10^-3.0 vs 10^-2.86);
#   • it is an independent match to Gregor et al. 2010, who measure half-maximal
#     SINGLE-CELL cAMP output at 0.5 ± 0.1 nM (rhythmic threshold ~1 nM);
#   • it keeps a quiescent rest state (ACA* ≤ 0.01 % of pool) and sits 17 % above
#     the bisected latch boundary (rest survives to 153.67, gone by 153.65) — the
#     usable window is [154, 194] and 180 is the point inside it that matches the rate.
# NOT free-standing: on its own this change LATCHES the population, because A_REF
# below rides on KTH_RC.  The two must move together — see the ADAPT_AREF note.
const KTH_RC     = 180.0f0   # RC_mem copy-number at half-max firing ≈ 0.73 nM cAMP_ext (un-adapted)
const N_ACA      = 4         # Hill cooperativity of the firing threshold (sharpness)
# [PAIRED WITH KTH_RC 2026-07-22] ADAPT_AREF 40 → 8.  A_REF is NOT independent of
# KTH_RC: the gate is
#     K_eff = KTH_RC·(1 + adapt/A_REF) = KTH_RC + (KTH_RC/A_REF)·adapt
# so the adapter's threshold SWING carries a factor KTH_RC.  Re-anchoring KTH_RC
# 470 → 180 divided that swing by 2.6 as a side effect, and the swing is exactly
# what makes co-firing dense cells adapt and fall quiet so the collective field can
# trough and reset.  Left at 40 the lawn LATCHED outright: population latchfrac
# 1.000 (PASS is < 0.30), 0 % of cells oscillating, ACA* pinned at 31 % of pool
# (validation/data/population_m8_diffusion.csv, tag `dense_v11_s1`).
#
# Measured, not guessed (validation/v11_abcb3.jl §G/§H/§I):
#   • it is NOT a drive problem — cutting secretion 9.2 ×0.5 / ×0.375 leaves
#     latchfrac at 1.000 in mean field AND on GPU (`kth180sec36_s1`), and raising
#     PdsA clearance ×3 / ×10 fails too.  A_REF is the ONLY knob of the three that
#     re-opens the window, which is the §1.8/Kamino-2017 mechanism doing its job.
#   • the naive ∝KTH_RC scaling (40·180/470 = 15.3) is NOT enough — still latched.
#     Restoring the absolute swing is insufficient because the field now sits 2.6×
#     further above threshold; the gate has to close HARDER than it used to.
#     At A_REF = 8 the swing is 22.5·adapt vs the old 11.75·adapt, i.e. stronger
#     fold-change detection — sharper on the rising edge, quieter on a background.
#   • mean-field cliff: 15.3 / 12 / 10 all latch (field 65 / 43 / 26 molec/vox);
#     8 resets (field 0.45) with the ACA* pulse still at 90.7 of pool.  Going lower
#     starts costing the pulse itself — at 4 the peak collapses 95 → 42.
#   • population, dense/both/D_cAMP=50, seed 1: latchfrac 0.000, 26 % of cells fire
#     ≥2×, period 159 s, wave 0.94 µm/s (vs 0.61 before) — `kth180aref8_s1`.
const ADAPT_AREF = 8.0f0     # was 40 (which was itself 60→40 in round 7, same mechanism)
# [EXAMINED AND HELD 2026-07-26, iteration 18]  k_acaOn stays at 10 µM⁻¹s⁻¹.
# DICTY_CELL.md §8 says the chain's on-rates "are scaled … so ACA reaches ~70 % of its
# pool" (criterion 8b.x#6, 65–75 %), and 8.11 is the only one of them that moves the ACA
# amplitude without disturbing the PKB timing that 8.x#1/8b.x#1 pin.  It is also the ONLY
# reachable knob for that criterion: the saturating fraction is
#     ACA*/pool = kA / (kA + k_acaOff + PKA-brake),   kA = k_acaOn·CRAC*/M_VOX
# so at k_acaOn = 10 even ALL 154 CRAC molecules being active gives kA = 0.32 ⇒ 72 %.
# No recruitment, phosphorylation or dephosphorylation rate can reach 65–75 %; only the
# gain can.  It was raised to 20 and MEASURED, and it must be put back, because the
# amplitude and the population reset are the same knob:
#   k_acaOn = 10  → ACA* 53.8 % of pool, GPU latchfrac 0.228 PASS, wave 7.2 µm/s
#   k_acaOn = 20  → ACA* 66.1 % of pool, GPU latchfrac 0.414 FAIL, wave 7.2 µm/s,
#                   sustained ACA* 13 % of pool (vs 6 %) — the lawn stops resetting
# 8b.x#6 (a "should", written against a baseline the harness itself flags as stale) is
# therefore in direct conflict with 8b.x#5 / 8.x#9 (no latch — "the primary constraint").
# The primary constraint wins.  See TUNING.md iteration 18 for the full A/B.
#
# ── [DEACTIVATED 2026-07-30] k_acaOn 10.0 → 0.0 ─────────────────────────────
# 8.11 is the ONLY route from Modules 1–8 into ACA, and it is switched OFF.
#
# WHY.  Modules 1–8 have not been re-calibrated against the recalibrated Module
# 9/10 (calib/CALIBRATION.md).  With 8.11 live, ACA has TWO independent drivers:
# the receptor arm (RC_mem → … → CRAC* → 8.11) and the excitable element
# (0.9f, X-driven and gated on ERK2*).  The receptor arm carries its own,
# uncalibrated gain and its own Hill gate (KTH_RC / ADAPT_AREF, both tuned at
# D_cAMP = 50 against the OLD Module 9/10), so in `DICTY_MODULES=full` it emits
# cAMP pulses that the Module-0/9/10 brake was never calibrated to terminate —
# a second pulse generator running in parallel with the calibrated one, not a
# refinement of it.  Zeroing the rate leaves the whole Module-8 cascade intact
# and observable (CRAC still cycles, PKB still fires) while making Module 0.9
# the sole ACA driver, which is the condition every number in
# calib/CALIBRATION.md was measured under.
#
# NOTE this changes NOTHING in `DICTY_MODULES=minimal` (the calibrated default):
# 8.11 is not in that bundle at all.  It only makes `full` mode behave as
# "calibrated pulse generator + observable but non-driving receptor cascade".
#
# TO RESTORE, once Module 8 is calibrated against the new Module 9/10:
# `DICTY_ACAON=10` (its previous committed value — see the iteration-18 note
# above for why 10 and not 20).  Expect to have to re-balance it against 0.9f,
# since the two arms now sum.
const k_acaOn = parse(Float32, get(ENV, "DICTY_ACAON", "0.0"))   # µM⁻¹s⁻¹  CRAC*-driven ACA activation
# Inactive by default (Module 8 is commented out of CELL_REACTIONS — see the
# DEACTIVATED note above), but kept receptor-model-aware for whenever it is
# re-enabled: in biswas mode there is no single :RC_mem, so the gate reads the
# sum of all five occupied classes (S included — it never desensitizes, but
# it is still "occupied") instead.
r8_11 = Reaction([:CRACa_mem, :ACA_mem], [:ACAa_mem, :CRACa_mem], k_acaOn,
                 BISWAS_RECEPTOR ?
                 @rate(x -> (x[:RC_H_mem] + x[:RC_L_mem] + x[:RC_S_mem] +
                             x[:RpC_H_mem] + x[:RpC_L_mem])^N_ACA /
                            ((KTH_RC * (1f0 + x[:RcptAdapta_cyto] / ADAPT_AREF))^N_ACA +
                             (x[:RC_H_mem] + x[:RC_L_mem] + x[:RC_S_mem] +
                              x[:RpC_H_mem] + x[:RpC_L_mem])^N_ACA)) :
                 @rate(x -> x[:RC_mem]^N_ACA /
                            ((KTH_RC * (1f0 + x[:RcptAdapta_cyto] / ADAPT_AREF))^N_ACA + x[:RC_mem]^N_ACA)))
const K_PKAPKB_BOOST = 30f0
const k_pkaPkb_2nd = Float32(K_PKAPKB_BOOST * min_to_s(0.6f0) / (PKA_pool_v / MOLEC_PER_µM))
r8_13 = Reaction([:PKAa_cyto, :PKBAa_mem], [:PKAa_cyto, :PKBA_cyto], k_pkaPkb_2nd, nothing)

const k_acaOff = 0.0064f0
r8_12 = Reaction([:ACAa_mem], [:ACA_mem], k_acaOff, nothing)

# ============================================================================
# Module 8b — PKBR1 (PKB-Related 1): Fast Membrane-Anchored PKB Isoform
# ----------------------------------------------------------------------------
# PKBR1 (Q55EQ8 / DDB_G0282133) lacks a PH domain and is instead
# constitutively palmitoylated at Cys-2, anchoring it permanently to the
# inner leaflet of the plasma membrane WITHOUT requiring PIP3 (Meili et al.
# 2000 Curr Biol 10:708; Kamimura et al. 2008 Curr Biol 18:1034).  Because
# PKBR1 is pre-positioned at the membrane, TORC2 (which is also constitutively
# membrane-resident) phosphorylates its hydrophobic motif (HM, Ser534 in Dd)
# within 5–10 s of cAMP stimulation — the tightest timing constraint in the
# model (Liu et al. 2017 J Cell Sci 130:1545; crit 8.x #1).
#
# PP2A-type phosphatases rapidly dephosphorylate PKBR1-HM (τ ≈ 2 s), making
# PKBR1 a FAST-ON / FAST-OFF pulse generator: it peaks early, decays within
# ~15 s, and leaves the refractory window open before the slow PKA brake (8.13)
# and the PKBA slow arm can accumulate.  This FAST DECAY is the key property
# that lets a fast cascade avoid population latch:
#   • PKBR1 (τ_off ≈ 2 s) → sharp CRAC pulse → ACA fires → cAMP_i spike
#   • Early cAMP_i → early PKA → inhibits PKBR1a (8b.4) AND PKBa (8.13)
#   • Net: pulse is shorter, not larger → reset window preserved
#
# PKBR1 and PKBA share substrates (CRAC, GefS, GAPA/B, talin, PI4P5K);
# PKA-dependent feedback inhibits both isoforms (Liu 2017).
#
# Cascade-speed tension resolution (TUNING.md iters 2–5): increasing PKBA
# alone (casc_mult) caused population latch because the slow PKBA off-rate
# (τ = 10 s) allowed sustained ACA drive.  PKBR1's 5× faster off-rate (τ = 2 s)
# means PKBR1a is cleared long before PKA activates, capping ACA drive at the
# pulse level rather than sustaining it.
#
# Pool split: PKBA (PKBA_cyto) is reduced from 154 → 77 molec/vox; PKBR1_mem
# is 77 molec/vox → total PKB capacity conserved at 154/vox.  The ~50:50
# split approximates the relative TORC2-mediated phosphorylation observed in
# Kamimura 2008 / Liu 2017 pulse-labelling experiments.
# ============================================================================
const k_pkbr1On  = 30.0f0   # µM⁻¹s⁻¹  TORC2→PKBR1-HM phosphorylation; 3× faster than PKBA
                              # (no cytosol→membrane translocation step; direct HM access)
const k_pkbr1Off = 0.5f0    # s⁻¹       PP2A dephosphorylation of PKBR1-HM; τ ≈ 2 s
                              # (5× faster than PKBA k_pkbOff = 0.1 s⁻¹; Kamimura 2008)
# 8b.1  TORC2*_mem + PKBR1_mem → TORC2*_mem + PKBR1a_mem
#   TORC2 phosphorylates constitutively-membrane PKBR1 at its HM; no recruitment step needed.
#   Rate 3× PKBA (8.7) reflecting the closer proximity / fewer diffusive steps.
r8b_1 = Reaction([:TORC2a_mem, :PKBR1_mem], [:PKBR1a_mem, :TORC2a_mem], k_pkbr1On,    nothing)
# 8b.2  PKBR1a_mem → PKBR1_mem   (PP2A rapid dephosphorylation, τ ≈ 2 s)
#   PP2A is constitutively active at the membrane and rapidly strips the PKBR1-HM
#   phosphate, making PKBR1 a pulse generator rather than a sustained activator.
r8b_2 = Reaction([:PKBR1a_mem], [:PKBR1_mem], k_pkbr1Off, nothing)
# 8b.3  PKBR1a_mem + CRAC_mem → CRACa_mem + PKBR1a_mem
#   PKBR1 activates the MEMBRANE-RECRUITED CRAC pool (same substrate as PKBA 8.9);
#   rate identical to 8.9 (k_cracOn = 10 µM⁻¹s⁻¹) — same phosphorylation chemistry.
#   [2026-07-22] educt CRAC_cyto → CRAC_mem with the three-state restructure: this is
#   the reaction that makes the pkbR1⁻/pkbA⁻ double null abolish ACA relay (Kamimura
#   2008; Cai 2010), because with BOTH kinases gone nothing converts CRAC_mem → CRAC*.
r8b_3 = Reaction([:PKBR1a_mem, :CRAC_mem], [:CRACa_mem, :PKBR1a_mem], k_cracOn, nothing)
# 8b.4  PKA*_cyto + PKBR1a_mem → PKA*_cyto + PKBR1_mem   (PKA refractory brake)
#   PKA inhibits PKBR1 activity by the same mechanism as it inhibits PKBA (8.13):
#   Liu 2017 shows PKA activity is required for the *transient* nature of HM
#   phosphorylation on both isoforms.  Reuses k_pkaPkb_2nd (gain-normalised form).
r8b_4 = Reaction([:PKAa_cyto, :PKBR1a_mem], [:PKAa_cyto, :PKBR1_mem], k_pkaPkb_2nd, nothing)

# ============================================================================
# Module 8.14 — THE OUTPUT EDGE: PIP3-gated relay into the excitable element
# [ADDED 2026-08-09 — calib/CALIBRATION_M8.md]
# ----------------------------------------------------------------------------
# WHAT THIS REPLACES.  8.11 (`CRAC*_mem + ACA_mem → ACA*_mem`) used to be the
# only route from Modules 1–8 into cAMP synthesis, and it is switched off
# (`k_acaOn = 0`) for the reason recorded at its definition: driving ACA
# DIRECTLY makes the receptor arm a SECOND, independent pulse generator running
# in parallel with the calibrated Module 0.9 one.  Two consequences follow from
# that topology and neither is fixable by tuning 8.11's gain:
#   (i)  the induced pulse has its own amplitude and its own shape, set by
#        CRAC*'s time course rather than by the calibrated X spike, so
#        "a relayed pulse ≈ the spontaneous pulse" is a coincidence to be
#        re-fitted every time either arm moves;
#   (ii) the two generators have SEPARATE refractory states — Module 0.9's Y
#        gates only the spontaneous arm, so a receptor-driven pulse can fire
#        into Module 0.9's refractory window and vice versa.  A refractory
#        period that must hold "regardless of source" is then not a property of
#        the cell at all, only of whichever arm happens to fire.
#
# THE FIX is to make the receptor arm a STIMULUS TO THE EXISTING EXCITABLE
# ELEMENT rather than a parallel output stage: CRAC* injects the Module 0.9
# activator X.  Everything downstream — the ERK2*→RegA disarm, the ACA gate
# (0.9f), the cAMP_i spike, the export, the PKA/RegA brake — is then literally
# the same machinery on the same rate constants, so
#   • an induced pulse IS a spontaneous pulse (same amplitude by construction,
#     not by fitting), and
#   • ONE recovery variable (Y, plus the shared Module 10 brake) sets the
#     refractory period for BOTH sources.  That is what makes the 8-minute
#     block source-independent instead of arm-specific.
# This is also the standard reading of the biology: the cAMP relay is an
# excitable system with one excitable core, and receptor occupancy is its
# stimulus, not a second oscillator.
#
# WHY THE GATE IS PIP3.  The requirement is that the cell has already resolved
# the DIRECTION of the signal before it re-amplifies it.  PIP3_mem is the
# output of the Module 3→4 LEGI/STEN stage, i.e. the first species in the
# cascade that is both (a) stimulus-driven and (b) POLARISED; gating the relay
# on it means the pathway from receptor to cAMP re-synthesis is open only when
# a front has actually formed.  Module 4's own refractory arm (4.10–4.14) then
# also contributes to the block: a cell that has just polarised cannot re-form
# a PIP3 domain immediately, so the gate is shut on top of Y being high.
# The Hill form (N_M8X_PIP3, K_M8X_PIP3) is a threshold on the per-voxel PIP3
# count, so a resting cell — PIP3 at its M4_INITIAL rest level — leaves the
# gate essentially closed and the relay silent.
#
# UNITS.  The engine's propensity for the single educt CRAC*_mem is
# rate·n(CRAC*_mem) per voxel; summed over the 13 voxels of the cell that is
# rate·CRAC*_total·gate molecules of X per second, and X is a per-CELL variable
# (see the M_CELL note at Module 0.9).  So the rate is written as
# "molecules of X per second per cell at FULL CRAC* activation" divided by the
# CRAC pool — a number directly comparable to Module 0.9's own k_fhn_k0
# (2.52 molec/s) and to the ~24-molecule excitation threshold.
# ============================================================================
const CRAC_v      = 154                      # CRAC pool molecules/voxel (all three states)
const CRAC_TOT    = CRAC_v * N_VOX_CELL      # 2002 CRAC / cell
# Injection gain: X molecules/s/cell when EVERY CRAC is active and the gate is
# fully open.  The excitation threshold of the Module 0.9 element is ≈24
# molecules of X against a baseline removal of k_fhn_k2 = 1 s⁻¹, so this has to
# clear ~24 molec/s at the CRAC* fraction the cascade actually reaches — and it
# must NOT clear the same threshold against the Y-mediated removal a cell
# carries in its refractory window (0.9y2 adds k_fhn_k5·Y ≈ 0.05·Y s⁻¹).  That
# two-sided constraint is what pins it; see calib/CALIBRATION_M8.md.
# ── WHERE THE GATE HAS TO SIT, AND WHY IT IS STEEP  [MEASURED 2026-08-09] ────
# calib/m8_probe.jl `dose` (six doses, 2 nM – 200 nM, one process each) measures
# PIP3_mem at the peak of the wave and again 150 s after the wave has gone:
#     dose[µM]  PIP3 peak/vox   PIP3 plateau/vox   ratio    CRAC* peak/plateau
#       0.002       17644            6578          2.68        54 % / 44 %
#       0.010       20593            6302          3.27        55 % / 43 %
#       0.060       25672            6168          4.16        56 % / 41 %
#       0.200       26781            6289          4.26        56 % / 42 %
# Two facts follow and both shape this gate:
#
# (1) THE RELAY IS ESSENTIALLY DOSE-INDEPENDENT.  2 nM already drives PIP3 to
#     17644/voxel and fires a full pulse; 100× more dose adds 50 %.  That is the
#     LEGI/STEN stage doing its job (it reports contrast, not level), and it is
#     what makes an all-or-none relay possible at the few-nM a neighbour 100 µm
#     away actually delivers.  So the gate must be OPEN at 17644, not at 25000.
#
# (2) MODULES 3/4 DO NOT RETURN TO REST.  Resting PIP3 is 1670/voxel; 150 s
#     after the stimulus it settles at ~6300, i.e. ~3.8× rest, and stays there.
#     Traced upstream this is Gβγ: r2_5 recapture is SECOND order, so its time
#     constant grows as the pool drains (τ = 132 s at the post-pulse level
#     against 8700 s at rest), and RasG-GTP and PIP3 inherit the tail.  This is
#     a Module 1–4 property, not a Module 8 one, and it is NOT fixed here —
#     re-opening it would re-open CALIBRATION_M12/M3/M4.  But it is the single
#     hardest constraint on this gate: a gate placed anywhere below the plateau
#     leaves a standing relay drive, and the cell then free-runs.  MEASURED at
#     the first-pass setting (K = 9000, n = 4): spikes at 130, 394, 650 s in the
#     refractory run and 70, 386, 684, 991 s in the recovery run — a ~290 s
#     self-sustained oscillation with no stimulus present.
#
# CONSEQUENCE: the gate has to separate 17644 from 6578 — a ratio of only 2.68 —
# by the ~300× that keeps the plateau drive well under Module 0.9's own seed
# rate (k_fhn_k0 = 2.52 molec/s).  A power law is the only thing that can, and
# n = 8 is what 2.68^n ≳ 300 requires.  That steepness is not a claim about one
# cooperative binding event: 8.14 is a LUMPED gate standing in for the whole
# PIP3 → PKB-coincidence → CRAC → cyclase-priming chain, several steps of which
# are themselves thresholded, and the model already uses n = 4 gates of exactly
# this kind (N_ACA, N_COMP, N_ERKGATE).
#
# At K = 13000, n = 8 the gate is 0.92 at the 2 nM peak, 0.998 at 60 nM,
# 4.3e-3 at the plateau and 1e-7 at rest.
const k_m8X_max   = parse(Float32, get(ENV, "DICTY_M8_XGAIN", "200.0"))  # molec X / s / cell
const N_M8X_PIP3  = parse(Int,     get(ENV, "DICTY_M8_NPIP3", "8"))      # gate Hill steepness
# [2026-08-16] 13000 -> 2500. ⚠ THIS CONSTANT WAS CALIBRATED AGAINST A DEFECT.
# CALIBRATION_M8 tuned it to separate a 17644 peak from a 6578 plateau — but
# both of those were inflated by Module 3's non-adapting pedestal (RasG_GTP_mem
# stuck at 5040/voxel instead of its 956 rest, so PIP3_mem floored at ~7150
# instead of ~1500). With the antithetic arm switched back on (M3_ARMS "both",
# 2026-08-16) PIP3_mem peaks at 5007/voxel under a FULL clamped stimulus, so the
# old gate evaluated to (13000/5007)^8 = 2087, i.e. 0.05 % open: measured on
# calib/m8_probe.jl the cell showed X = 0 and cAMP_i = 0 through an entire
# stimulus, and the 250-cell population collapsed from 100 % to 13 % of cells
# firing >=2x because nothing relayed.
#
# 2500 is the geometric mean of the 5007 stimulated peak and the ~1500 rest.
#
# ⚠ [2026-08-23] STALE, UN-REVALIDATED SINCE THE ANTITHETIC ARM WAS REMOVED.
# This whole constant was chosen specifically to work WITH the antithetic
# controller, whose textbook-exact setpoint is what kept RasG_GTP_mem's
# population-context rest near 956/voxel instead of the 5040/voxel pedestal the
# note above measures for CORE+LEGI alone under REPEATED, incompletely-reset
# population stimulation (Gbg_cyto floors at 1699/voxel between relay waves,
# 11x past the LEGI arms' crossover — see the removed `_m3_inhib`/`M3_ARMS`
# block's retained history in this file's git log, "CORE+LEGI has NO
# stimulus-independent setpoint"). That is exactly the class of limit Takeda et
# al. (2012) name for feedforward adaptation — "inherently a tuned ratio, not
# an architectural guarantee" — and exactly the scenario (dense, rapid,
# incompletely-resetting relay, as opposed to the isolated single-pulse/step
# lab protocols Takeda et al. actually tested) where it is most likely to bite.
# Whether K_M8X_PIP3 = 2500 still separates peak from pedestal correctly now
# that the arm holding the pedestal down is gone has NOT been re-measured —
# re-run `calib/m8_probe.jl` and a population run before trusting this gate.
# ── ✓ [2026-09-17] BOTH RE-RUN; SEE THE 2026-09-16 BLOCK BELOW. The answer was
# that the pedestal this warning was about had been closed elsewhere (the M3
# kinetics retune), and the gate had gone stale in the OPPOSITE direction --
# four decades shut at the dose the relay exists to serve. This note is kept
# because it called the right constant; its conclusion is superseded.
# What makes it the right value rather than merely a working one: it returns the
# relay to CALIBRATION_M8's OWN measured numbers, which is not something a free
# re-tune would do —
#     ACA* peak    94 % of pool   (target 94 %)
#     exported     8.28e6 molec   (target ~9e6)
# so the gate constant was the only thing that had gone stale.
#
# Checked for self-triggering, since lowering a threshold is exactly the change
# that caused the documented ~290 s free-running of an earlier first-pass gate.
# calib/m8_probe.jl rest, seed 1, 1800 s, three arms:
#     K=13000 relay live        2 spikes    (no relay: gate shut)
#     K=2500  relay off XGAIN=0 2 spikes    (no relay: output cut)
#     K=2500  relay LIVE        2 spikes    <- the one that matters
# against 2.41 expected from the calibrated 748 s ISI. The relay adds no
# spontaneous firing. Refractory probe also PASSes (blocked at 105 s and 405 s).
# [2026-08-21] 2500 -> 13000, i.e. BACK to its calibrated value. It had been
# dropped to 2500 to reach a PIP3 signal that was ~7x too small; k_pip3 = 6.2
# above restores the level (plateau 7230, peak 18833/voxel) so the gate no longer
# has to chase a collapsed species. At 2500 against the restored level the gate
# would sit BELOW the resting PIP3 and be permanently open — the two constants
# MUST move together.
# ── ⚠ [RE-ANCHORED 2026-09-16] 13000 → 5000.  THE PLATEAU THIS GATE WAS BUILT
# TO CLEAR NO LONGER EXISTS, AND 13000 NOW SITS ABOVE THE SIGNAL. ────────────
# Everything above is a negotiation with ONE defect: "Modules 3/4 do not return
# to rest", PIP3 settling at ~6300/voxel = 3.8× its resting 1670 after a pulse.
# That pedestal is what forced K up to 13000 and n to 8, and it is what the
# note above calls "the single hardest constraint on this gate".
#
# IT IS GONE.  `calib/m8_probe.jl dose` on HEAD, one cell, field clamped, the
# PIP3 PEAK during the stimulus against the PLATEAU 150 s after it has ended —
# the two numbers this gate has to separate, measured the same way as the 2026-
# 08-09 table above:
#     dose[µM]   PIP3 peak/vox   PIP3 plateau/vox   ratio
#       0.002         5026              1072         4.7
#       0.005         7170               930         7.7
#       0.010         5940              1040         5.7
#       0.020         9625              1160         8.3
#       0.060        16813              1007        16.7
#       0.200        21311              1127        18.9
# The plateau is 930-1160 at EVERY dose, i.e. back at the solved resting value
# (993/voxel, `analyse/cell_rest_state.jl`) rather than 3.8× above it.  The
# 2026-09-06 Module 3 kinetics retune (sustained-step plateau 54.4 % → 3.5 % of
# peak) and the STEN brake are what closed it; this gate was never re-examined
# against that.
#
# WHAT IT COSTS TO LEAVE IT AT 13000.  The peak the relay actually has to catch
# is the one a NEIGHBOUR delivers, and CALIBRATION_M8.md's whole point is that
# this is a few nM: at 2 nM the peak is 5026/voxel, so the gate evaluates to
# (13000/5026)^-8 = 5e-4 — FOUR DECADES SHUT at the dose the relay exists to
# serve.  13000 was chosen when the same 2 nM stimulus drove PIP3 to 17644,
# 3.5× higher; the amplitude moved and the gate did not.
#
# WHERE IT GOES.  Same rule the 2500 value used — half-open between the peak the
# relay must catch and the plateau it must ignore — but with today's numbers:
# 5000 sits just under the 2 nM peak and 4.5× above the plateau, giving
#     gate @ 2 nM peak (5026)  = 0.51        gate @ plateau (1100) = 5.5e-6
#     gate @ 60 nM peak (16813) = 0.9999     gate @ rest (993)     = 2.4e-6
# so the relay drive `XGAIN·f(CRAC*)·gate` is 200·0.44·0.51 = 45 molec X/s at
# the weakest dose — against Module 0.9's ≈24-molecule excitation threshold —
# and 4.7e-4 molec/s standing, i.e. 5000× below 0.9a's own seed (2.52/s).  The
# geometric mean (≈2350) would open the gate wider still at 2 nM, and is NOT
# taken: it buys nothing the 2× margin above does not already have, and gives
# away the headroom that protects against a pedestal returning in a POPULATION
# context, where the cells sit in a cAMP background this clamped single-cell
# probe does not reproduce.
#
# MEASURED END-TO-END at this value (`calib/m8_probe.jl single`, with
# `k_cracDeph` re-anchored — the two changes are one repair and neither works
# alone), two doses two decades apart:
#     0.002 µM → FIRED, X_pk 447, cAMP_i 386390, ACA* 94 % of pool, exported 8.9e6
#     0.060 µM → FIRED, X_pk 468, cAMP_i 380992, ACA* 94 % of pool, exported 8.6e6
# ⚠ THE PEDESTAL IS A POPULATION PROPERTY AND THIS CLAMPED PROBE CANNOT SEE IT.
# ── ✓ [2026-09-17] MEASURED IN A POPULATION, AND THE GATE DOES NOT FREE-RUN ──
# 1000 cells / 400^2 / 3000 s, seeds 42/43/44. The two things that would show a
# standing `gate(PIP3_pedestal)` are the period falling below the refractory
# block and the field failing to clear between waves. Neither happens:
#     period          467 +/- 16 s (kymograph) / 480 s (driver), vs a 300-420 s block
#     field trough    whole-line mean p10 = 1.2 / 1.3 / 1.5 nM
# i.e. the trough goes BELOW the 2.5 nM relay threshold in every seed, so the
# population genuinely resets rather than sitting in a bath.
# ⚠ Do NOT read the probe's own "line over threshold in 95 % of records" as
# contradicting that: it counts ANY of 57 line positions being above threshold,
# and with several sources firing per episode something is almost always up
# somewhere. The whole-line MEAN is the statistic that answers "is the
# population bathed"; the per-position count is not.
# If a future population run DOES show intervals well below the block, the
# standing gate is still the first thing to measure and this is still the knob.
# `DICTY_M8_KPIP3=13000` restores the old value.
const K_M8X_PIP3  = parse(Float32, get(ENV, "DICTY_M8_KPIP3", "5000.0")) # PIP3_mem/voxel at half-open gate  [13000 -> 2500 -> 13000 -> 5000]
const KM8P_N      = K_M8X_PIP3^N_M8X_PIP3                                 # baked into the @rate AST
# 8.14  CRAC*_mem → CRAC*_mem + X   (gated on PIP3_mem AND on the refractory clock)
#   The second factor is the SAME gate 0.9a carries (Module 0.9r), which is what
#   makes the 8-minute block source-independent: a cell that has just fired —
#   spontaneously or by relay — cannot be re-fired by either route.
#
#   ⚠ THE PIP3 GATE IS WRITTEN IN RECIPROCAL FORM, AND MUST STAY THAT WAY.
#   The algebraically identical x^n/(K^n + x^n) HANGS THE SIMULATION. Propensity
#   arithmetic is Float32 (core_single/propensity.jl converts every count with
#   `Float32(agent_state[j])`), and at n = 8 the numerator x^8 overflows Float32
#   for x > floatmax(Float32)^(1/8) = 65 536. `PIP3_mem` reaches ~48 800 on the
#   FIRST pulse and, because Modules 3/4 do not return to rest afterwards (it
#   settles ~3.8× rest — see CLAUDE.md), crosses 65 536 on the SECOND. Then
#   x^8 = Inf, the gate is Inf/(K + Inf) = Inf/Inf = NaN, and the SSA loop in
#   gillespie_voxel! never terminates: `a0 <= 0` is false for NaN and `t > dt`
#   is false for NaN, so neither break fires. Measured symptom: the GPU pins at
#   100 %, no error is raised, and the run stops advancing forever at a
#   DETERMINISTIC step (525 969 = t 1051.9 s ≈ one inter-pulse interval, at
#   100 cells/640²/seed 42).
#
#   1/(1 + (K/x)^n) has the same value everywhere and cannot overflow into a
#   NaN: for x ≥ 1 the base K/x ≤ K, and for x = 0 it is Inf, giving Inf^n = Inf
#   and a gate of exactly 0 — which is the correct limit. n = 8 is the only
#   exponent in this network whose Float32 ceiling (65 536) sits inside the
#   range its species actually visit; every other gate here is n ≤ 4, i.e. a
#   ceiling of 4.3e9. Any NEW gate with n ≥ 6 needs this form too.
r8_14 = Reaction([:CRACa_mem], [:CRACa_mem, :X], Float32(k_m8X_max / CRAC_TOT),
                 @rate(x -> (1f0 / (1f0 + (K_M8X_PIP3 / x[:PIP3_mem])^N_M8X_PIP3)) *
                            (KR_N / (KR_N + x[:Refr_cyto]^N_REFR))))

# ── the Module 8 bundle ─────────────────────────────────────────────────────
# WHAT IS IN, AND WHAT IS DELIBERATELY OUT.
#
# OUT — 8.7 / 8.8 / 8.13, the "Module 8 PKB" arm.  These read `PKBA_cyto` /
# `PKBAa_mem`, which are Module 4's PkbA species (4.10–4.13), and Module 4's arm
# is the calibrated one: PIP3 recruitment THEN a kinase-gated hydrophobic-motif
# step — the published PkbA coincidence detection — on a pool of
# PKBA_v = 10^4/voxel.  8.7 is a second, PIP3-INDEPENDENT activation route for
# the same molecule at k = 10 µM⁻¹s⁻¹ on a 77/voxel pool.  Running both means
# (a) two irreconcilable pool sizes for one species and (b) a TORC2-driven
# PKBAa_mem flood into Module 4's r4_14 feedback, which would re-open
# CALIBRATION_M4.md.
#
# ⚠ AND 8.7 IS MIS-TYPED, WHICH IS THE DEEPER REASON IT STAYS OUT: a single
# TORC2-gated step with NO recruitment stage is PkbR1's mechanism, not PkbA's
# (Meili 2000; Kamimura 2008; Kamimura 2010 — see the naming block above 4.10).
# So 8.7 modelled the wrong paralog's kinetics on the right paralog's species.
# Module 8b (8b.1–8b.4) already supplies PkbR1 properly, with its own `PKBR1_v`
# pool, its own 5× faster off-rate and its own PKA brake — so re-enabling 8.7
# would not "restore the missing isoform", it would double-count it.
#
# Division of labour, therefore: **Module 4 owns PkbA, Module 8b owns PkbR1**,
# and both feed the shared substrates.  8.9 reads `PKBAa_mem` and 8b.3 reads
# `PKBR1a_mem` onto the same `CRAC_mem` pool, exactly as Kamimura 2008 requires;
# 4.14 and 4.14b likewise both feed the Charest 2010 Sca1 feedback.  Because
# `PKBAa_mem` is itself PIP3-recruited, the 8.9 edge is a SECOND PIP3 dependence
# in the relay chain on top of 8.9b and the 8.14 gate — while 8b.3 is the
# PIP3-independent one, which is the whole point of having both.
#
# OUT — 8.11, superseded by 8.14 above (it stays defined, and `DICTY_ACAON`
# still re-enables it, so the old direct-to-ACA topology remains reproducible).
const MODULE_8_ON = get(ENV, "DICTY_MODULE8", "1") == "1"
const MODULE_8_REACTIONS = MODULE_8_ON ? [
    r8_1, r8_2, r8_3, r8_4, r8_4b, r8_4c, r8_4d,   # RasGEFA → RasC-GTP + its GAP adaptation
    # Sca1 feedback on RasC (Charest 2010): own species since 2026-10-02 (see 8.2s); DICTY_M8_SCA1FB=0 removes it
    (k_sca1Fb > 0 ? vcat(BRAKE_W_PKBR1 > 0 ? [r8_2p] : Reaction[],
                         (M4_REFRACTORY && BRAKE_W_PKBA > 0) ? [r8_2pA] : Reaction[],
                         [r8_2q, r8_2s]) : Reaction[])...,
    r8_5, r8_6,                                     # RasC → TORC2*
    r8b_1, r8b_2, r8b_3, r8b_4,                     # PKBR1 arm (M8's own kinase) + PKA brake
    r8_9, r8_9b, r8_9c, r8_10, r8_10b,              # the three-state CRAC moiety
    r8_14,                                          # PIP3-gated relay into X
    (k_acaOn > 0 ? [r8_11] : Reaction[])...,        # legacy direct-to-ACA edge, off by default
] : Reaction[]

# ── Module 8 pool sizes, initial state and diffusion ────────────────────────
# Pools are the pre-existing values carried by dicty_minimal.jl/mini_simulation.jl
# (the drivers that last instantiated Module 8), centralised here the same way
# M3_INITIAL/M4_INITIAL centralise theirs.  `PKBA_cyto`/`PKBAa_mem` are absent on
# purpose — they are PkbA and Module 4 owns them (see the bundle note above);
# what lives here is PkbR1's `PKBR1_v` pool, a different protein.
const GEFA_v    = parse(Int, get(ENV, "DICTY_M8_GEFAV", string(M8_GEFA_PAIR ? GEFA_v_pair : GEFA_v_legacy)))  # RasGEFA (Aimless) / voxel
const RasC_v    = 385        # RasC / voxel
const TORC_v    = 154        # TORC2 / voxel
const PKBR1_v   = parse(Int, get(ENV, "DICTY_M8_PKBR1V", "77"))   # PKBR1 / voxel (0 = pkbR1⁻)
const M8_INITIAL = Dict{Symbol, Int}(
    :RasGEFA_cyto  => GEFA_v,   :RasGEFAa_cyto => 0,
    :RasC_GDP_mem  => RasC_v,   :RasC_GTP_mem  => 0,
    :RasCGAP_cyto  => RasCGAP_v, :RasCGAPa_cyto => 0,
    :TORC2_mem     => TORC_v,   :TORC2a_mem    => 0,
    :PKBR1_mem     => PKBR1_v,  :PKBR1a_mem    => 0,
    :CRAC_cyto     => CRAC_v,   :CRAC_mem      => 0,  :CRACa_mem => 0,
    :Sca1_cyto     => Sca1_v,   :Sca1a_mem     => 0,  # Sca1-like PKB substrate (8.2p/q/s), all unphosphorylated at rest
)
const M8_DIFFUSION = Dict{Symbol, Float32}(
    :Sca1_cyto     => D_M3_BRAKE, :Sca1a_mem   => D_M3_BRAKE,  # as the former shared brake species (RasC arm unchanged)
    :RasGEFA_cyto  => 10.0f0,   :RasGEFAa_cyto => 2.0f0,   # as CELL_DIFF's generic cytosolic
    :RasCGAP_cyto  => 10.0f0,   :RasCGAPa_cyto => 10.0f0,
    :RasC_GDP_mem  => 0.1f0,    :RasC_GTP_mem  => 0.1f0,   # prenylated, as RasG (M3_DIFFUSION)
    :TORC2_mem     => 0.1f0,    :TORC2a_mem    => 0.1f0,   # membrane complex
    :PKBR1_mem     => 0.1f0,    :PKBR1a_mem    => 0.1f0,   # palmitoylated anchor
    :CRAC_cyto     => 10.0f0,   :CRAC_mem      => 0.1f0,  :CRACa_mem => 0.1f0,
    :ACA_mem       => 0.1f0,    :ACAa_mem      => 0.1f0,   # membrane cyclase — was MISSING
)

# ============================================================================
# Module 9 — cAMP Synthesis, Export and the PKA brake on ACA
# ----------------------------------------------------------------------------
# [RECALIBRATED 2026-07-30 — calib/CALIBRATION.md Iterations 7, 9]
#
# AMPLITUDE.  The synthesis rate is anchored to the measured relay output, not
# to Laub–Loomis' normalised k9/k11.  Devreotes, Derstine & Steck 1979 (J Cell
# Biol 80:291) measure ~3×10⁶ cAMP molecules released per cell per stimulus in a
# perfusion assay; Grutsch & Robertson 1978 measure 3×10⁷ at a saturating
# stimulus, with amplification complete within 8 s.
#
# [RECALIBRATED 2026-07-31 — Iteration 14] 30 → 20 s⁻¹ total, split 8/12 (same
# 0.4/0.6 ratio as before).  The task specifies a per-ACA turnover "in the
# order of 20 cAMP/s" directly, superseding the Devreotes-anchored 30 s⁻¹; the
# 10⁷-molecule amplitude target is met instead by the ACA pool (r0_9f note)
# with this rate held at the task's value.  20 s⁻¹ is still an ordinary
# adenylyl-cyclase k_cat, and the ~40× discrepancy the old comment noted
# against the pre-recalibration relay-tuning value (k11_eff ≈ 0.9 s⁻¹) stands
# unchanged.
const k9_eff  = 8.0f0    # s⁻¹  9.1 Hunger-scaled synthesis arm
const k11_eff = 12.0f0   # s⁻¹  9.2 constant synthesis arm  (k9+k11 = 20 s⁻¹)
# 9.3 PKA* ⊣ ACA*, LL k2 = 0.9 min⁻¹ in effective-rate (PKA-pool-normalised) form.
const k2_eff  = min_to_s(0.9f0)
const k2_2nd  = Float32(k2_eff / (PKA_pool_v / MOLEC_PER_µM))   # ≈ 0.047 µM⁻¹s⁻¹
# ── 9.2b AbcB3 export ───────────────────────────────────────────────────────
# [RECALIBRATED] 3.0 → 0.5 s⁻¹.  The 3.0 was never a measurement: it was derived
# as k10_eff·k11_eff/k9_eff purely to preserve a secreted:degraded partition
# under the OLD Module 10, an anchor this recalibration removes.  Miranda et al.
# 2015 (Dev Biol 397:203) establish only that AbcB3 efflux is first order in
# cAMP_i.  The binding constraints now are (i) τ_export ≪ pulse width — 2 s
# against a ~40 s pulse, so the transporter still adds no dynamical delay — and
# (ii) Dinauer 1980, cells export most of the cAMP they make, which the RegA
# gate delivers where it should: the exported fraction is 0.05 at rest (RegA
# armed, cAMP_i held near zero) and 0.83 at the pulse peak (RegA disarmed).
# At 3.0 s⁻¹ export outran RegA by 250× and RegA had NO effect on the pulse at
# all — measured: knocking out the whole ERK2→RegA edge changed the exported
# cAMP by a factor of 1.0.
const k_abcB3 = 0.5f0     # s⁻¹  AbcB3 cAMP efflux (first order in cAMP_i)
# 9.1 ACA*_mem → ACA*_mem + cAMP_i  (Hunger scales synthesis up; §0.8)
r9_1 = Reaction([:ACAa_mem], [:ACAa_mem, :cAMP_i], k9_eff,
                @rate(x -> 1f0 + α_H * x[:Hunger] / Float32(H_max_vox)))
# 9.2 ACA*_mem → ACA*_mem + cAMP_i  (constant arm)
r9_2 = Reaction([:ACAa_mem], [:ACAa_mem, :cAMP_i], k11_eff, nothing)
# 9.2b cAMP_i → cAMP_ext  (AbcB3 transporter; Miranda et al. 2015)
r9_2b = Reaction([:cAMP_i], [:cAMP_ext], k_abcB3, nothing)
# 9.3 PKA*_cyto + ACA*_mem → PKA*_cyto + ACA_mem  (negative feedback)
r9_3 = Reaction([:PKAa_cyto, :ACAa_mem], [:PKAa_cyto, :ACA_mem], k2_2nd, nothing)
# 8.12 ACA*_mem → ACA_mem  (intrinsic deactivation)
# [RECALIBRATED] 0.1 → 0.0064 s⁻¹, i.e. t½ = 108 s, which is Cai et al. 2010's
# measured ACA* decline.  At 0.1 s⁻¹ (τ = 10 s) ACA* was slaved to the X spike
# and collapsed with it, so the pulse carried only 8×10⁵ molecules.  At 0.0064
# the brief X spike CHARGES ACA* and the pulse length is set by the
# ACA*/PKA/RegA loop instead of by X.  (The realised decline is faster than
# 108 s — t½ ≈ 37 s — because PKA* strips ACA* via 9.3 on top of this rate;
# Cai's 108 s is measured under a SUSTAINED stimulus, where the PKA arm is at
# steady state rather than pulsing.)

# ============================================================================
# Module 10 — the intracellular brake: RegA (3-state) / ERK2 / PKA
# ----------------------------------------------------------------------------
# Structure follows Module10.md Part III (10.B/10.C) — RegA is a conserved
# three-state moiety with an exit from every state, replacing the old
# unbounded ∅→RegA source + absorbing RegA_inakt pair:
#
#   RegA_cyto  (RegA_u : unphosphorylated at D212, activity 1/20)
#      ⇅  10.B1 / 10.B2   phosphorelay (RdeA-H65 ↔ RegA-D212, Thomason 1998/99)
#   RegAp_cyto (RegA_p : phospho-D212, ACTIVE — the brake)
#      ⇅  10.B3 / 10.B4   ERK2* phosphorylation ↔ dephosphorylation
#   RegAp_i_cyto (RegA_pi : ERK2-phosphorylated, INACTIVE)
#
# Conserved: RegA_cyto + RegAp_cyto + RegAp_i_cyto = RegA_pool_v per voxel.
# 10.B4 is the reaction the old Module 11 never had (Module10.md §I.3) and is
# what terminates each pulse: it is the primary handle on the falling phase.
# ============================================================================
# [RECALIBRATED 2026-07-30] RegA pool 154 → 1540/voxel (0.032 → 0.32 µM).
# Module10.md Part VIII lists RegA's copy number as NEVER MEASURED; the 0.032 µM
# figure comes from Laub–Loomis' normalised units, not from an experiment.  At
# 0.032 µM, RegA needs k_cat/K_M ≈ 3×10⁸ M⁻¹s⁻¹ — at or beyond the diffusion
# limit — before it can compete with AbcB3 for cAMP_i, i.e. before the regA⁻
# phenotype (elevated cAMP, precocious development: Shaulsky 1998, Thomason
# 1998) is reproducible at all.  0.32 µM puts RegA on the same pool scale as
# ERK2 and PKA and brings the required efficiency down to 3×10⁷ M⁻¹s⁻¹.
const RegA_pool_v   = 1540
const k_relayOn     = 0.05f0     # s⁻¹  10.B1 lumped Dhk→RdeA→RegA-D212 transfer [ESTIMATE]
const k_relayOff    = 0.0056f0   # s⁻¹  10.B2 chosen for a 0.9 resting active fraction
# [RECALIBRATED] LL's k8 = 1.3 min⁻¹ is a coefficient in a NORMALISED-ACTIVITY
# ODE, not a molecular rate (Module10.md §I.1), so it carries no unit-bearing
# information.  Used literally it disarms RegA with τ = 144 s while ERK2* is
# only up for ~45 s, and the brake never opens (measured knockdown 1.3×).
# 8.0 µM⁻¹s⁻¹ gives τ_disarm ≈ 0.5 s at the ERK2* peak, so the brake is released
# inside the X spike — the ordering the pulse requires.  Measured knockdown 64×.
const k_regAInh     = 8.0f0      # µM⁻¹s⁻¹ 10.B3 ERK2* inactivates RegA_p
const k_regAReact   = 0.0167f0   # s⁻¹     10.B4 RegA_pi → RegA_p, τ = 60 s (LL k7 = 1.0 min⁻¹)
# [RECALIBRATED] 0.417 → 30 µM⁻¹s⁻¹ = k_cat/K_M = 3×10⁷ M⁻¹s⁻¹, i.e. k_cat ≈
# 150 s⁻¹ at the measured K_M = 5 µM (Thomason 1998) — fast, but comfortably
# sub-diffusion-limited.  0.417 is LL's k10 read as if it were a molecular rate
# and leaves RegA 250× weaker than AbcB3 export, i.e. dynamically absent.
const k_regACat     = 30.0f0     # µM⁻¹s⁻¹ 10.C1 RegA_p hydrolyses cAMP_i
const REGA_U_FRAC   = 0.05f0     # RegA_u activity relative to RegA_p (Thomason 1999: ≥20-fold)
# PKA: activation is cAMP binding to the R subunit, K_d ≈ 0.1 µM for Dd PKA-R,
# so k_on/k_off = 1/K_d with k_off from LL's k4 = 1.5 min⁻¹ (τ = 40 s), which is
# also what sets the minutes-scale refractory period.
const k_pkaOn       = 0.25f0     # µM⁻¹s⁻¹ 10.D1
const k_pkaOff      = 0.025f0    # s⁻¹     10.D2  τ = 40 s
const k_pkaErk      = Float32(min_to_s(0.8f0) / (PKA_pool_v / MOLEC_PER_µM))  # 10.D3, LL k6 = 0.8 min⁻¹

# 10.B1 RegA_cyto → RegAp_cyto            (phosphorelay activation)
r10_B1 = Reaction([:RegA_cyto], [:RegAp_cyto], k_relayOn, nothing)
# 10.B2 RegAp_cyto → RegA_cyto            (phosphorelay reversal)
r10_B2 = Reaction([:RegAp_cyto], [:RegA_cyto], k_relayOff, nothing)
# 10.B3 RegAp_cyto + ERK2a_cyto → RegAp_i_cyto + ERK2a_cyto  (ERK2 ⊣ RegA; Maeda 2004)
r10_B3 = Reaction([:RegAp_cyto, :ERK2a_cyto], [:RegAp_i_cyto, :ERK2a_cyto], k_regAInh, nothing)
# 10.B4 RegAp_i_cyto → RegAp_cyto         (THE reset reaction — Module10.md §I.3)
r10_B4 = Reaction([:RegAp_i_cyto], [:RegAp_cyto], k_regAReact, nothing)
# 10.C1 cAMP_i + RegAp_cyto → RegAp_cyto  (catalytic hydrolysis, cAMP_i ≪ K_M ⇒ linear)
r10_C1 = Reaction([:cAMP_i, :RegAp_cyto], [:RegAp_cyto], k_regACat, nothing)
# 10.C2 cAMP_i + RegA_cyto  → RegA_cyto   (residual 1/20 activity of the unphospho form)
r10_C2 = Reaction([:cAMP_i, :RegA_cyto], [:RegA_cyto], Float32(k_regACat * REGA_U_FRAC), nothing)
# 10.D1 cAMP_i + PKA_cyto → cAMP_i + PKA*_cyto   (cAMP_i catalytic; PKA pool conserved)
r10_D1 = Reaction([:cAMP_i, :PKA_cyto], [:cAMP_i, :PKAa_cyto], k_pkaOn, nothing)
# 10.D2 PKA*_cyto → PKA_cyto
r10_D2 = Reaction([:PKAa_cyto], [:PKA_cyto], k_pkaOff, nothing)
# 10.D3 PKA*_cyto + ERK2a_cyto → PKA*_cyto + ERK2_cyto   (PKA ⊣ ERK2; LL k6)
r10_D3 = Reaction([:PKAa_cyto, :ERK2a_cyto], [:PKAa_cyto, :ERK2_cyto], k_pkaErk, nothing)

# ============================================================================
# Module 11 — extracellular PdsA phosphodiesterase (Module10.md Part III /
# Module11.md, implemented 2026-08-05, fixes defects P1–P9)
# ----------------------------------------------------------------------------
# Full pathway: PdsA_i_cyto (translated) → PdsA_mem (surface-displayed) →
# PdsA_ext (shed, world-shared, diffuses). PdiA_ext is a real 1:1 stoichiometric
# inhibitor (⇌ PdsA_PdiA_ext complex), not a catalytic converter — fixes P4.
# Degradation now covers the complex too — fixes P6/P7 (no more false
# "conserved moiety"; see calib/CALIBRATION_M11.md for the flux-balance check
# that replaces it).
#
# Developmental gate (P_dev, fixes P8): Module11.md §I.1/§I.2 found pdsA/pdiA
# are PULSE-INDEPENDENT starvation genes — induced by density/PSF within ~30
# min of starvation onset (Van Driessche et al. 2002; Faure et al. 1988/1990),
# NOT by cAMP-receptor/pulse signalling. A direct search of Cai et al. (2014)
# GtaC ChIP-seq targets found no pdsA hit, so unlike the old model there is NO
# GtaC or cAMP-pulse-decoder term here — gate on Hunger alone. The window
# OPENS EARLY relative to Hunger's own hours-long ramp to H_max_vox (H_DEV_HALF
# ≪ H_max_vox, matching the fast ~30 min literature onset) and does not close:
# closing is a Weening et al. (2003) hours-to-days effect outside the timescale
# any current driver exercises Module 11 over (documented scope limit).
# ============================================================================
# ── ENV override helper (same convention as Module 4's `_m4`) ───────────────
# Every Module 11 rate constant is overridable from the environment so that
# calib/m11_pop.jl can sweep them without editing this file. Defaults below are
# the committed values; an unset variable therefore changes nothing.
_m11(env, default) = parse(Float32, get(ENV, env, string(default)))

const H_DEV_HALF = parse(Float32, get(ENV, "DICTY_DEV_HALF", "25"))   # Hunger/voxel at half-induction [ESTIMATE]
const N_DEV      = 2                                                    # P_dev Hill steepness
const HD_N       = Float32(H_DEV_HALF)^N_DEV                            # baked into the @rate AST
const _P_dev     = @rate(x -> x[:Hunger]^N_DEV / (HD_N + x[:Hunger]^N_DEV))

# ── 11.D — PdsA abundance, secretion, localisation (per cell) ───────────────
# k_pdsaBas + k_pdsaInd calibrated (Module10.md Part III worked example) to
# reproduce the model's previous [PdsA_ext] ≈ 2 µM at MIN_DIST=6 (≈1 cell/36
# voxels): 173 molecules·cell⁻¹·s⁻¹ ÷ 13 voxels ÷ MOLEC_PER_µM ≈ 2.76e-3 µM/s/
# voxel total at full induction; 10 %/90 % basal/induced split. [ESTIMATE beyond
# that anchor — see calib/CALIBRATION_M11.md for the tuned values.]
# [2026-08-15] COMMITTED G = 16: both constants ×16 over the Module10.md Part III
# anchor (2.8e-4 / 2.5e-3). At G = 1 the population fires as ONE SYNCHRONOUS
# GLOBAL FLASH — measured k_eff = 9.6e-4 s⁻¹ ⇒ λ = sqrt(D_cAMP/k_eff) = 600 µm,
# which at any density used here reaches the whole field, so mean and max
# cAMP_ext rise together, there is no front, and NO cell ever resets below the
# 12 molec/voxel firing threshold (final whole-field mean 781.7/voxel, 65×
# threshold ⇒ 0 % of cells fire ≥2×). At G = 16, k_eff = 0.0154 s⁻¹
# (τ_clear = 65 s), λ = 151 µm ≈ 1.2 cell spacings at MIN_DIST = 50, final mean
# 9.4/voxel (BELOW threshold), 100 % of cells fire ≥2× at a 780–792 s period,
# and the kymograph shows travelling waves at ~6–7 µm/s. The requirement is
# λ ≈ 1–1.5 spacings, which is nearly DENSITY-INVARIANT (available k_eff ∝ ρ,
# required k_eff ∝ 1/spacing² ∝ ρ) ⇒ G ≈ 8–32 over the whole density range.
# See calib/CALIBRATION_M11_POP.md §6. Revert with DICTY_M11_KBAS/KIND.
const k_pdsaBas  = _m11("DICTY_M11_KBAS",  4.48f-3) # µM/s/voxel  constitutive floor (P_dev = 0)
const k_pdsaInd  = _m11("DICTY_M11_KIND",  4.0f-2)  # µM/s/voxel  induced amplitude  (P_dev = 1)
const k_pdsaSec  = _m11("DICTY_M11_KSEC",  1.0f-2)  # s⁻¹  PdsA_i_cyto → PdsA_mem (τ ≈ 100 s) [ESTIMATE, Module10.md 10.D2]
const k_pdsaDeg  = _m11("DICTY_M11_KDEG",  5.0f-4)  # s⁻¹  turnover, all PdsA/PdiA protein states [ESTIMATE, Module10.md 10.D4]
const k_pdsaShed = _m11("DICTY_M11_KSHED", k_pdsaDeg) # s⁻¹  PdsA_mem → PdsA_ext; 50:50 partition start (Palsson 2009 knob)

# 11.D1  ∅ → PdsA_i_cyto   (basal floor + Hunger-gated induction; fixes P1/P2/P8)
r11_D1a = Reaction(Symbol[], [:PdsA_i_cyto], k_pdsaBas, nothing)
r11_D1b = Reaction(Symbol[], [:PdsA_i_cyto], k_pdsaInd, _P_dev)
# 11.D2  PdsA_i_cyto → PdsA_mem   (secretion to the cell surface; fixes P3)
r11_D2  = Reaction([:PdsA_i_cyto], [:PdsA_mem], k_pdsaSec, nothing)
# 11.D3  PdsA_mem → PdsA_ext   (shedding; writes into the cell's own world
#         voxel via world_shared_species, same mechanism cAMP_ext/9.2b use)
r11_D3  = Reaction([:PdsA_mem], [:PdsA_ext], k_pdsaShed, nothing)
# 11.D4  turnover — pure agent-internal states, no world copy needed
r11_D4a = Reaction([:PdsA_mem], Symbol[], k_pdsaDeg, nothing)
r11_D4b = Reaction([:PdsA_i_cyto], Symbol[], k_pdsaDeg, nothing)
# 11.D4 (world-shared states) — dual copy: agent-side (this cell's voxels) +
#        world-side (free voxels, below) so clearance happens everywhere,
#        fixing P6 (old code degraded PdsA_ext but never the inactive complex)
r11_D4c = Reaction([:PdsA_ext], Symbol[], k_pdsaDeg, nothing)
w11_D4c = Reaction([:PdsA_ext], Symbol[], k_pdsaDeg, nothing)
r11_D4d = Reaction([:PdsA_PdiA_ext], Symbol[], k_pdsaDeg, nothing)
w11_D4d = Reaction([:PdsA_PdiA_ext], Symbol[], k_pdsaDeg, nothing)

# ── 11.E — PdiA inhibition (reversible, stoichiometric; fixes P4/P5) ────────
# 11.F — cAMP_slow: a per-voxel low-pass filter of local cAMP_ext, so PdiA
# repression tracks the pulse ENVELOPE rather than individual pulses
# (Module10.md 10.E1). k_campTrack is grounded in real kinetics: the PdiA
# transcript is repressed/derepressed by cAMP with <30 min response time in
# BOTH directions (PDI-cloning paper, PubMed 2169446; Module11.md §II.3) —
# unity-gain filter, τ ≈ 1/k_campTrack ≈ 42 min, same order as that anchor.
const k_campTrack = _m11("DICTY_M11_KTRACK", 4.0f-4)   # s⁻¹  cAMP_slow tracking/relaxation rate
r11_F1 = Reaction([:cAMP_ext], [:cAMP_ext, :cAMP_slow], k_campTrack, nothing)   # catalytic: cAMP_ext not consumed
r11_F2 = Reaction([:cAMP_slow], Symbol[], k_campTrack, nothing)

# PdiA is co-induced by the SAME starvation gate as PdsA (Module11.md §II.4:
# PSF induces PdsA and PdiA together), repressed by cAMP_slow (Yeh et al. 1978).
const k_pdiSyn = _m11("DICTY_M11_KPDISYN", 2.5f-3)  # µM/s/voxel  induced amplitude, same order as k_pdsaInd (co-induced)
const K_pdi    = _m11("DICTY_M11_KPDI",    5.0f-2)  # µM  cAMP_slow half-repression [ESTIMATE]
const N_PDI    = 2
const KPDI_N   = K_pdi^N_PDI
# 11.E1  ∅ → PdiA_ext   (agent-only: synthesis is cell-local)
r11_E1 = Reaction(Symbol[], [:PdiA_ext], k_pdiSyn,
                  @rate(x -> (x[:Hunger]^N_DEV / (HD_N + x[:Hunger]^N_DEV)) *
                             (1f0 / (1f0 + x[:cAMP_slow]^N_PDI / KPDI_N))))
# 11.E2  PdiA_ext → ∅   (dual copy — same slow protein turnover as PdsA)
r11_E2 = Reaction([:PdiA_ext], Symbol[], k_pdsaDeg, nothing)
w11_E2 = Reaction([:PdiA_ext], Symbol[], k_pdsaDeg, nothing)
# 11.E3  PdsA_ext + PdiA_ext ⇌ PdsA_PdiA_ext   (dual copy — real 1:1 binding,
#         titratable, replacing the old catalytic "PdsA + PDI → PdsA_inakt + PDI")
const k_pdiOn  = _m11("DICTY_M11_KPDION",  1.0f0)   # µM⁻¹s⁻¹  [ESTIMATE, K_d ≈ 1 µM target]
const k_pdiOff = _m11("DICTY_M11_KPDIOFF", 1.0f0)   # s⁻¹  fast, near-equilibrium binding (protein-protein
                              # complex formation, not an enzymatic step) — must be
                              # fast relative to k_pdsa's few-second clearance timescale
                              # or the inhibition never has time to act before a cAMP
                              # pulse is already cleared (calib/CALIBRATION_M11.md,
                              # calib/m11_probe.jl TEST titration)
r11_E3f = Reaction([:PdsA_ext, :PdiA_ext], [:PdsA_PdiA_ext], k_pdiOn, nothing)
w11_E3f = Reaction([:PdsA_ext, :PdiA_ext], [:PdsA_PdiA_ext], k_pdiOn, nothing)
r11_E3r = Reaction([:PdsA_PdiA_ext], [:PdsA_ext, :PdiA_ext], k_pdiOff, nothing)
w11_E3r = Reaction([:PdsA_PdiA_ext], [:PdsA_ext, :PdiA_ext], k_pdiOff, nothing)
# 11.E4  cAMP_ext + PdsA_PdiA_ext → PdsA_PdiA_ext   (dual copy — residual
#         catalysis: PdiA raises K_M ~1000×, doesn't abolish activity)
const k_pdsa = _m11("DICTY_M11_KPDSA", min_to_s(4.9f0) * 1.5f0)   # µM⁻¹s⁻¹ effective cAMP_ext hydrolysis
const k_pdsaResidual = k_pdsa / 1000f0
r11_E4 = Reaction([:cAMP_ext, :PdsA_PdiA_ext], [:PdsA_PdiA_ext], k_pdsaResidual, nothing)
w11_E4 = Reaction([:cAMP_ext, :PdsA_PdiA_ext], [:PdsA_PdiA_ext], k_pdsaResidual, nothing)
# 11.2   cAMP_ext + PdsA_ext → PdsA_ext   (dual copy — primary clearance, kept
#         from the pre-refactor code)
r11_2 = Reaction([:cAMP_ext, :PdsA_ext], [:PdsA_ext], k_pdsa, nothing)
w11_2 = Reaction([:cAMP_ext, :PdsA_ext], [:PdsA_ext], k_pdsa, nothing)

# ── Density scaling of the extracellular PdsA field [2026-08-10] ────────────
# `PdsA_ext` is produced per CELL and lost per VOXEL, so its steady-state level
# — and with it the ENTIRE clearance strength of the model — is proportional to
# the cell density ρ [cells/voxel]. Two consequences that a population driver
# must not get wrong:
#
#  1. There is no such thing as "the" PdsA_ext warm-start level. The 1 µM that
#     drivers used to hardcode is the value at Module10.md Part III's worked
#     density (1 cell per 36 voxels). At a population run's typical 10⁻³ – 10⁻²
#     cells/voxel the cells' own synthesis sustains 10–30× less, so a 1 µM seed
#     is an over-clearing transient that decays away on a slow 1/k_pdsaDeg ≈
#     2000 s timescale — i.e. the run is non-stationary for longer than it lasts,
#     with the signal suppressed for the part that matters. Seed the DENSITY-
#     MATCHED value (`pdsa_steady_µM`) and the field is stationary from t = 0.
#
#  2. k_eff = k_pdsa·[PdsA_ext] ∝ ρ, hence the signal range
#     λ = sqrt(D_cAMP/k_eff) ∝ ρ^(−1/2) — exactly the scaling of the mean cell
#     spacing dx·ρ^(−1/2). λ expressed IN CELL SPACINGS is therefore DENSITY-
#     INVARIANT, and is the one number that says whether a pulse reaches the
#     next cell. It is set by the rate constants alone, not by the geometry, so
#     it is a genuine Module 11 tuning target (see calib/CALIBRATION_M11.md).
const _f_sec  = k_pdsaSec  / (k_pdsaSec  + k_pdsaDeg)   # fraction of PdsA_i reaching the membrane
const _f_shed = k_pdsaShed / (k_pdsaShed + k_pdsaDeg)   # fraction of PdsA_mem shed to the world
P_dev(hunger) = hunger^N_DEV / (HD_N + hunger^N_DEV)
"PdsA_ext molecules shed into the world per second by one fully-induced cell."
pdsa_source_per_cell(; hunger = H_max_vox, n_vox = 13) =
    n_vox * (k_pdsaBas + k_pdsaInd * P_dev(Float32(hunger))) * Float32(MOLEC_PER_µM) *
    _f_sec * _f_shed
"Steady-state [PdsA_ext] in µM at a cell density ρ [cells per voxel]."
pdsa_steady_µM(ρ; hunger = H_max_vox, n_vox = 13) =
    ρ * pdsa_source_per_cell(; hunger, n_vox) / (k_pdsaDeg * Float32(MOLEC_PER_µM))

"""
    m11_world_steady(ρ; hunger = H_max_vox, n_vox = 13)

Steady state of the three EXTRACELLULAR Module 11 states, in molecules per
voxel, at cell density `ρ` [cells/voxel]: `(PdsA_ext, PdiA_ext, PdsA_PdiA_ext)`.

Both proteins are shed per cell and lost per voxel at the same `k_pdsaDeg`, so
their TOTAL pools are `ρ·source/k_pdsaDeg`; the 11.E3 binding equilibrium then
partitions each total between free and complexed. Solving the equilibrium here
rather than seeding free protein only matters because `k_pdiOn`/`k_pdiOff`
equilibrate in ~1 s while the pools themselves take ~2000 s — seed the pools
wrong and the binding follows instantly, seed the binding wrong and it corrects
instantly, but seed only free `PdsA_ext` (as drivers did) and the run starts
with more ACTIVE enzyme than the network sustains.

PdiA repression by `cAMP_slow` is evaluated at `cAMP_slow = 0` (a rested
population, which is the state every driver starts from).
"""
function m11_world_steady(ρ; hunger = H_max_vox, n_vox = 13)
    P_tot = ρ * pdsa_source_per_cell(; hunger, n_vox) / k_pdsaDeg
    I_tot = ρ * n_vox * k_pdiSyn * P_dev(Float32(hunger)) * Float32(MOLEC_PER_µM) / k_pdsaDeg
    Kd    = Float64(k_pdiOff / k_pdiOn) * Float64(MOLEC_PER_µM)   # molecules/voxel
    b     = Float64(P_tot) + Float64(I_tot) + Kd
    C     = (b - sqrt(max(b^2 - 4*Float64(P_tot)*Float64(I_tot), 0.0))) / 2
    return (Float64(P_tot) - C, Float64(I_tot) - C, C)
end

"""
    m11_world_profile(gw, gh, lo_i, hi_i, lo_j, hi_j, S_tot; D, k, dx, bc) -> Matrix

The SUSTAINED extracellular enzyme field, in molecules per voxel, on a
`gw × gh` lattice whose edge is ABSORBING (default; `bc = :reflecting` or
`:periodic` solves the matching non-leaking world boundary instead), fed by a total source `S_tot`
[molec/s] spread uniformly over the voxel box the cells occupy.

⚠ WHY `m11_world_steady` IS NOT ENOUGH, AND WHY A SINGLE NUMBER CANNOT BE.
`m11_world_steady` divides the source by `k_pdsaDeg` alone, i.e. it solves a
CLOSED system. The engine's world edge is absorbing (`core_single/diffuse_kernel.jl`
gates every jump on `nExists` and lets a molecule leave), and `PdsA_ext` has
D = 65 µm²/s against `k_pdsaDeg` = 5e-4 s⁻¹, so its decay length is
    sqrt(D/k) = 361 µm
which is COMPARABLE TO THE WHOLE DOMAIN at the grid sizes this repo runs. The
edge, not degradation, is then the dominant sink, and the closed-system answer
is too high by a factor that depends on the geometry:

    geometry                seeded   sustained   over-seed   λ/spacing
    400²  margin 20  1000c   15700       3133       5.0×     1.33 → 2.97
    1000² margin 55  4000c   10048       6502       1.5×     1.34 → 1.67
    1400² margin 100 10⁴ c   12816      11471       1.1×     1.39 → 1.47

⚠ SO "λ IN CELL SPACINGS IS DENSITY-INVARIANT, SET BY THE RATE CONSTANTS ALONE
AND NOT BY THE GEOMETRY" IS FALSE AS IMPLEMENTED — it is 2.97 / 1.67 / 1.47
across those three. It is true of the RATES; it is not true once the domain is
finite and open. Quote λ per geometry.

And the sustained field is NOT UNIFORM: centre/edge-of-cells is 4.4–4.7×. A
uniform seed of any value therefore leaves a transient, and the one it leaves is
slow — building the profile means transporting enzyme to the edges, which takes
L²/D ≈ 10⁴ s, LONGER THAN THE RUN. Hence a profile, not a number.

Verified three ways at 400²/1000 cells: SOR relaxation on the lattice 2596,
this spectral solution 2618, and the grid mean the engine itself ends at 2616.

The solution is the sine expansion, exact for a rectangular Dirichlet box:
    n = Σ_mn  s_mn / (D((mπ/L_x)² + (nπ/L_y)²) + k) · sin(mπx/L_x) sin(nπy/L_y)
evaluated as U'·A·U so the cost is O(N²G), milliseconds at N = 161 modes.
"""
function m11_world_profile(gw::Int, gh::Int, lo_i::Int, hi_i::Int, lo_j::Int, hi_j::Int,
                           S_tot::Real; D::Real = 65.0, k::Real = Float64(k_pdsaDeg),
                           dx::Real = 2.0, nmode::Int = 161, bc::Symbol = :absorbing)
    Lx, Ly = gw*dx, gh*dx
    ax, bx = (lo_i-1)*dx, hi_i*dx
    ay, by = (lo_j-1)*dx, hi_j*dx
    nsrc   = (hi_i-lo_i+1)*(hi_j-lo_j+1)
    nsrc <= 0 && error("m11_world_profile: empty source box")
    s = S_tot/nsrc                                   # molec/voxel/s
    bc === :absorbing || return _m11_profile_noflux(gw, gh, ax, bx, ay, by, s, D, k, dx, nmode, bc,
                                                    (lo_i, hi_i, lo_j, hi_j))
    ms = collect(1:2:nmode)
    Ix = [(Lx/(m*pi))*(cos(m*pi*ax/Lx) - cos(m*pi*bx/Lx)) for m in ms]
    Iy = [(Ly/(n*pi))*(cos(n*pi*ay/Ly) - cos(n*pi*by/Ly)) for n in ms]
    A  = [ ((4/(Lx*Ly))*s*Ix[p]*Iy[q]) /
           (D*((ms[p]*pi/Lx)^2 + (ms[q]*pi/Ly)^2) + k) for p in eachindex(ms), q in eachindex(ms) ]
    Ux = [sin(m*pi*((i-0.5)*dx)/Lx) for m in ms, i in 1:gw]
    Uy = [sin(n*pi*((j-0.5)*dx)/Ly) for n in ms, j in 1:gh]
    max.(Ux' * A * Uy, 0.0)        # the truncated series can ring slightly negative
end

# [2026-09-23] The same steady state under a NON-LEAKING world boundary
# (`Environment.world_boundary = :reflecting | :periodic`, see
# calib/CALIBRATION_BOUNDARY_2026-09-23.md). A reflecting lattice edge is a
# Neumann wall at the voxel face, so the eigenfunctions are cosines, m = 0
# included — the m = 0 term is the closed-system mean S_tot/(k·gw·gh), which is
# why `m11_world_steady` becomes EXACT again when the cells fill the grid. A
# periodic grid with a box centred in it is even about the centre, hence even
# about every face, hence the same cosine series; an off-centre box on a torus
# is a different problem and is refused rather than silently mis-solved.
function _m11_profile_noflux(gw, gh, ax, bx, ay, by, s, D, k, dx, nmode, bc, box)
    bc in (:reflecting, :periodic) || error("m11_world_profile: unknown bc = :$bc")
    if bc === :periodic
        lo_i, hi_i, lo_j, hi_j = box
        (lo_i - 1 == gw - hi_i && lo_j - 1 == gh - hi_j) ||
            error("m11_world_profile(bc = :periodic): the cell box must be centred " *
                  "(margins $(lo_i-1)/$(gw-hi_i) × $(lo_j-1)/$(gh-hi_j))")
    end
    Lx, Ly = gw*dx, gh*dx
    ms = collect(0:nmode)
    Ix = [m == 0 ? (bx - ax) : (Lx/(m*pi))*(sin(m*pi*bx/Lx) - sin(m*pi*ax/Lx)) for m in ms]
    Iy = [n == 0 ? (by - ay) : (Ly/(n*pi))*(sin(n*pi*by/Ly) - sin(n*pi*ay/Ly)) for n in ms]
    εx = [m == 0 ? 1/Lx : 2/Lx for m in ms]
    εy = [n == 0 ? 1/Ly : 2/Ly for n in ms]
    A  = [ (εx[p]*εy[q]*s*Ix[p]*Iy[q]) /
           (D*((ms[p]*pi/Lx)^2 + (ms[q]*pi/Ly)^2) + k) for p in eachindex(ms), q in eachindex(ms) ]
    Ux = [cos(m*pi*((i-0.5)*dx)/Lx) for m in ms, i in 1:gw]
    Uy = [cos(n*pi*((j-0.5)*dx)/Ly) for n in ms, j in 1:gh]
    max.(Ux' * A * Uy, 0.0)
end

"""
    m11_world_seed(gw, gh, lo_i, hi_i, lo_j, hi_j, ncell; hunger, n_vox, bc)

Per-voxel seed for the three extracellular Module 11 states, as
`(PdsA_ext, PdiA_ext, PdsA_PdiA_ext)` matrices in molecules/voxel: the sustained
PdsA and PdiA profiles from `m11_world_profile`, partitioned voxel by voxel
through the same 11.E3 binding equilibrium `m11_world_steady` uses. The binding
equilibrates in ~1 s while the pools take ~2000 s, so seeding the pools right
and letting the binding follow is what matters — but the partition is level
dependent, so it has to be done per voxel rather than once at the mean.
"""
function m11_world_seed(gw::Int, gh::Int, lo_i::Int, hi_i::Int, lo_j::Int, hi_j::Int,
                        ncell::Int; hunger = H_max_vox, n_vox = 13, bc::Symbol = :absorbing)
    P_src = ncell * pdsa_source_per_cell(; hunger, n_vox)
    I_src = ncell * n_vox * k_pdiSyn * P_dev(Float32(hunger)) * Float32(MOLEC_PER_µM)
    P = m11_world_profile(gw, gh, lo_i, hi_i, lo_j, hi_j, P_src; bc)
    I = m11_world_profile(gw, gh, lo_i, hi_i, lo_j, hi_j, I_src; bc)
    Kd = Float64(k_pdiOff / k_pdiOn) * Float64(MOLEC_PER_µM)
    C  = similar(P)
    @inbounds for idx in eachindex(P)
        b = P[idx] + I[idx] + Kd
        C[idx] = (b - sqrt(max(b^2 - 4*P[idx]*I[idx], 0.0)))/2
    end
    (P .- C, I .- C, C)
end

# ── Module 11 agent-internal rest state, DERIVED (not hand-set) ─────────────
# `CELL_INITIAL` carried no Module 11 entry at all, so every driver started its
# cells with PdsA_i_cyto = PdsA_mem = 0 and had to spend ~1/k_pdsaDeg ≈ 2000 s
# — longer than most population runs — filling the secretion chain before the
# cells contributed any enzyme of their own. Same defect, and same fix, as the
# rest-state derivations in CALIBRATION_M34.md/M3/M4: solve the chain instead of
# guessing a number.
#   PdsA_i_cyto : source/(k_pdsaSec + k_pdsaDeg)
#   PdsA_mem    : k_pdsaSec·PdsA_i_cyto/(k_pdsaShed + k_pdsaDeg)
#   cAMP_slow   : unity-gain low-pass of cAMP_ext ⇒ 0 in a rested population
const _pdsa_i_ss = (k_pdsaBas + k_pdsaInd * P_dev(Float32(H_max_vox))) *
                   Float32(MOLEC_PER_µM) / (k_pdsaSec + k_pdsaDeg)
const M11_INITIAL = Dict{Symbol, Int}(
    :PdsA_i_cyto => round(Int, _pdsa_i_ss),
    :PdsA_mem    => round(Int, k_pdsaSec * _pdsa_i_ss / (k_pdsaShed + k_pdsaDeg)),
    :cAMP_slow   => 0,
)

# ============================================================================
# Reaction bundles consumed by dicty_simulation.jl
# ----------------------------------------------------------------------------

# Module 0.8 Hunger + 0.9 excitable element and its output stage
const MODULE_0_REACTIONS = [
    r0_8a, r0_8b,
    r0_9a, r0_9b, r0_9c, r0_9d, r0_9y1, r0_9y2, r0_9y3,
    r0_9r1, r0_9r2,                     # [2026-08-09] the shared refractory clock
    r0_9e, r0_9e2, r0_9f,
]
# Module 9 cAMP synthesis/export + Module 10 RegA/ERK2/PKA brake
const MODULE_910_REACTIONS = [
    r8_12,
    r9_1, r9_2, r9_2b, r9_3,
    r10_B1, r10_B2, r10_B3, r10_B4, r10_C1, r10_C2,
    r10_D1, r10_D2, r10_D3,
]
# Module 1 (receptor) and Module 2 (G-protein cycle), named separately so the
# calibration harnesses in calib/ can build a receiver cell out of exactly these
# two modules without dragging in the rest of the cascade.
const MODULE_1_REACTIONS = BISWAS_RECEPTOR ? MODULE_1_REACTIONS_BISWAS :
    [r1_1, r1_2, r1_3, r1_4, r1_5, r1_5b, r1_6, r1_7, r1_8a, r1_8b]
# [2026-08-07] `r2_1` used to appear TWICE in this list (commit 39aa40a,
# "tuning") — the engine sums propensities channel by channel, so a duplicated
# entry is nothing but a silent ×2 on `k_gef`, applied where no one reading the
# rate constants would see it and where `calib/consistency.jl` cannot check it.
# The factor is now folded into `k_gef` itself (M2_GEF_GAIN above), so the
# effective rate is visible at its definition and ENV-overridable like every
# other constant.
const MODULE_2_REACTIONS = BISWAS_RECEPTOR ? MODULE_2_REACTIONS_BISWAS :
    [r2_0, r2_1, r2_2, r2_3, r2_5]
# [2026-08-06] restored — dropped by the "new receptors" commit (1496f2e) that
# added Module 1b; calib/m12_probe.jl and calib/m3_probe.jl both build their
# receiver cell out of this constant and were left broken (UndefVarError)
# until this fix.
const MODULE_12_REACTIONS = vcat(MODULE_1_REACTIONS, MODULE_2_REACTIONS)
# ── Module 3: the complete network — incoherent feedforward, no other arm ────
# [2026-08-23] Module 3 is now exactly TWO superimposed mechanisms, both always
# on, matching the incoherent-feedforward (IFFL) architecture in Takeda K et al.
# (2012) "Incoherent feedforward control governs adaptation of activated Ras in
# a eukaryotic chemotactic signaling pathway." Sci Signal 5(205):ra2 — the
# direct experimental/modelling precedent for RasG-GTP adaptation in
# Dictyostelium — and in the spatial extension Shi C, Huang CH, Devreotes PN,
# Iglesias PA (2013) PLoS Comput Biol 9(7):e1003122 / the Gα2-Ric8
# direction-sensing model (PMC4859573). See DICTY_CELL.md's own Module 3
# section, which specifies this same architecture and cites Takeda 2012:
#
#   CORE  3.1/3.2/3.3 + 3.5/3.7 — the GEF arm (Gβγ → RasGEFRa → RasG-GTP,
#         FAST) and the GAP arm's own deactivation/hydrolysis. On its own this
#         is a monotone follower with no adaptation of any kind.
#   LEGI  3.4/3.4e/3.4b — the GLOBAL inhibitor (RasGAPa_cyto, D = 20 µm²/s,
#         SLOWER than the GEF arm by construction — see 3.4/3.5's rate note).
#         This is the entire adaptation mechanism: a fast local excitor and a
#         slower global inhibitor racing on the same substrate is what makes
#         RasG-GTP a transient differentiator instead of a level detector, and
#         it is precisely Takeda et al.'s finding, not an approximation of it.
#
# A THIRD arm — a two-species antithetic integral feedback (AIF) controller,
# reactions numbered 3.8-3.14 — was carried here 2026-07-31–2026-08-22 and has
# been DELETED, not merely disabled; see the citation block above the reaction
# definitions (search "Module 3 is IFFL-ONLY") for the literature comparison
# that motivated removing it and what it cost while present. `DICTY_M3_ARMS`
# is gone with it — there is only one network now, so nothing to switch.
#
# `DICTY_M3_INHIB` (rc | gbg | both) remains: it selects which upstream signal
# arms the LEGI inhibitor, a question orthogonal to whether an AIF arm exists.
# The Module 4-coupled reactions (3.3b-3.4d) stay out of this list — they need
# `PIP3_mem`/`PI34P2_mem`, which the Module 3 calibration harnesses do not
# instantiate, and they are Module 3/4 coupling rather than Module 3 itself.
# ── [2026-08-26] 3.3b + ITS BRAKE RESTORED.  This is the STEN's amplifier. ──
# They were commented out of this list in e991744 and have not run since, while
# CLAUDE.md and this file both went on describing 3.3b as "the spatial knob".
#
# WHY IT HAS TO BE HERE AND NOT IN MODULE 4.  Module 4's PIP3/PTEN mutual
# inhibition CANNOT be the amplifier at any parameter setting: its loop gain is
# arm1*arm2 with arm1 = PTEN's share of PIP3 removal and arm2 =
# n*H*(1-H)*k_ptenDispC/off, and those two trade off against each other — pushing
# PTEN_v up wins arm1 but saturates PTEN and kills arm2.  Scanned over
# PTEN_v x{4..30} × k_ship/{1..10} × N_DISP{2,3,4,6} × K_DISP the gain TOPS OUT
# AT 0.918 (calib/analysis_2026-08-26/amplifier_audit.jl §4/§5).  ⚠ Raising
# N_DISP makes it WORSE, not better — a steeper gate drives H away from 0.5
# where H(1-H) peaks — so this repo's earlier "n>=3 required" note does not
# apply to this gate.
#
# The literature puts the gain at Ras for the same reason: Kataria 2013 PNAS
# ("symmetry breaking occurs between heterotrimeric G protein signaling and Ras
# activation"), Fukushima 2019 (excitable Ras dynamics trigger self-organised
# PIP3), and Xu 2005 (G-protein activation is graded and shallow across the
# WHOLE cell surface while PHcrac/PIP3 is spatially amplified into a crescent).
# LEGI alone provides no amplification — 3.4/3.4e set the SPATIAL comparison,
# 3.3b provides the gain, 3.3c-e and the PKB/TORC2 arm terminate it.
#
# ⚠ THE BRAKE IS NOT OPTIONAL.  3.3b alone is a positive feedback with only Ras
# substrate depletion to stop it; 3.3c-e is what makes the result EXCITABLE
# (fires, amplifies, RESETS) rather than LATCHED, which is the failure mode this
# repo has already measured once (calib/CALIBRATION_M4.md, "test that it turns
# off").  They are restored together and DICTY_M3_STEN=0 removes both.
const M3_STEN = get(ENV, "DICTY_M3_STEN", "1") == "1"
const MODULE_3_REACTIONS = vcat(
    [r3_1, r3_2, r3_3, r3_5, M3_ZOU ? r3_7z : r3_7],                 # CORE (3.7z: NF1 zero-order, opt-in)
    M3_INHIB in ("rc", "both")  ? [r3_4]  : Reaction[],               # LEGI inhibitor, RC_mem-driven
    M3_INHIB in ("gbg", "both") ? [r3_4e] : Reaction[],               # LEGI inhibitor, Gβγ-driven (default)
    [r3_4b],                                                          # basal arming, either way
    # PIP3 amplifier: 3.3b PIP3 -> RasG-GTP (DEFAULT since 2026-10-02) or 3.3g PIP3 -> Gβγ release (DICTY_M3_PIP3PATH=gbg)
    M3_STEN ? (M3_PIP3PATH == "ras" ? [r3_3b] : [r3_3g]) : Reaction[],
    (M3_STEN && M3_PIP3GAP) ? [r3_4d] : Reaction[],     # 3.4d PIP3 -| RasGAP* (2nd amplifier, OFF)
    # 3.3c/d/e the STEN brake — the delayed negative feedback 3.3b has been
    # running without.  ON by default since 2026-09-06 (DICTY_M3_BRAKE=0 removes
    # it); only meaningful with the amplifier present, so gated on M3_STEN too.
    # ⚠ THIS IS NOT A REINSTATEMENT OF THE ANTITHETIC (AIF) ARM, WHICH STAYS
    # DELETED.  The AIF arm was an INTEGRAL controller: a global, memory-carrying
    # species pair with a setpoint (k_aifRef/k_aifSense) that integrated error to
    # zero.  Takeda et al. 2012 rejected exactly that for Dictyostelium Ras — it
    # produces spurious oscillations and stimulus-dependent kinetics — in favour
    # of the incoherent feedforward loop 3.1–3.7.  The brake has NO setpoint and
    # NO integrator: it is a local, membrane-bound, first-order negative feedback
    # (activation ∝ PIP3, decay τ = 10 s) acting on the STEN amplifier 3.3b, of
    # the same class as Module 0.9's validated FitzHugh–Nagumo inhibitor.
    # Adaptation still comes from the IFFL and only from the IFFL; the brake's
    # job is to keep 3.3b's autocatalysis from latching.  Checked for the failure
    # mode Takeda names: `SINGLE` (no 2nd excursion) PASSES at 29.8 %, and
    # m3_probe reports 0 sign changes in the smoothed post-peak decay — no
    # oscillation.
    # [2026-10-02] reader: RasG-GTP (3.3c_R, DEFAULT), PKB (3.3c′, =pkb) or PIP3 (3.3c, =pip3).
    # With the RasG-GTP reader the activation is ∝ RasG-GTP rather than PIP3; the rest of
    # this comment (written for the PIP3 reader) still holds: first order, no setpoint.
    (M3_STEN && M3_BRAKE) ? vcat(
        M3_BRAKE_READ == "ras" ? [r3_3cG] :
        M3_BRAKE_READ == "pkb" ?
            vcat((MODULE_8_ON && BRAKE_W_PKBR1 > 0) ? [r3_3cR] : Reaction[],
                 (M4_REFRACTORY && BRAKE_W_PKBA > 0) ? [r3_3cA] : Reaction[]) :
            [r3_3c],
        [r3_3d, r3_3e]) : Reaction[],
)
# Module 11 PdsA/PdiA — the agent-side half (synthesis/secretion/shedding are
# cell-only; the six dual-copy reactions also run on free voxels via
# WORLD_REACTIONS below).
const MODULE_11_REACTIONS = [
    r11_D1a, r11_D1b, r11_D2, r11_D3, r11_D4a, r11_D4b, r11_D4c, r11_D4d,
    r11_F1, r11_F2,
    r11_E1, r11_E2, r11_E3f, r11_E3r, r11_E4,
    r11_2,
]

# Modules 1–8: receptor → G-protein → Ras → PIP3 → Rac → Ca → RasC/TORC2/PKB/CRAC
# ============================================================================
# MODULE MEM — a PIP3-written directional MEMORY (opt-in, DICTY_MEMORY=1), 2026-09-27
# ----------------------------------------------------------------------------
# WHAT IT IS.  A membrane "memory mark" that is WRITTEN where PIP3 is high and
# ERASED slowly, and nothing else: no adaptation, no feedback into Ras, PIP3 or
# any other module — the rest of the network cannot see it.  The motility reads
# it (DICTY_MOVE_MARKER=mem), so the cell steers on WHERE PIP3 HAS BEEN over
# the last minutes instead of on the instantaneous PIP3 snapshot.
#
#   M.1  Mem_cyto → Mema_mem      k_memOn · H(PIP3_mem),  H = P^n/(K^n+P^n)
#   M.2  Mema_mem → Mem_cyto      k_memOff
#
# Mem_cyto is cytosolic (D = 10, mixed over the cell in ~2.5 s), Mema_mem is
# membrane-bound and nearly immobile (D = 0.01), so the per-voxel Mema_mem
# pattern is the time integral of the local H(PIP3), leaking with τ = 1/k_memOff.
# Its centroid therefore averages the PIP3 direction over the writing window and
# HOLDS it after PIP3 falls (a uniform first-order decay leaves the centroid
# unchanged) until the mark fades.
#
# THE BIOLOGY — what is established, and what is not (calib/PERIODIC_MEMORY_
# 2026-09-26.md §5, all sources re-read):
#   * a memory exists: Dictyostelium keeps its direction ~2 min in the back of
#     6–10 min waves and ≥ 5 min after a gradient reversal (Skoge et al. 2014
#     PNAS 111:14448); polarity "can be maintained in the absence of marked Ras
#     activation" (Nakajima et al. 2014 Nat Commun 5:5367);
#   * its carrier is UNKNOWN ("the molecular basis of cellular memory remains to
#     be determined", Skoge 2014); it is not the Ras-GTP level in the back of the
#     wave and it does not need F-actin (latrunculin, Skoge 2014).
#   So Mem/Mema is a PHENOMENOLOGICAL species standing for that unknown element.
#   It is written by PIP3 here (the user's specification) — note that no verified
#   candidate is both PIP3-driven and long-lived (PKB substrates return to
#   baseline in ~2 min, Kamimura 2008 / Tang 2011; PI(3,4)P2 is anti-correlated
#   with PIP3, Li 2018), and that Skoge's M additionally feeds back onto Ras,
#   which this module deliberately does not.
#
# CONSTANTS.
#   MEM_V      5000 molec/voxel      [ESTIMATE]  pool; large so the mark's centroid
#                                                 is not shot-noise limited
#   K_MEM      5000 molec/voxel      [CALIBRATED] PIP3 half-write level: above the
#                                                 resting PIP3 (1 200–1 600) and
#                                                 inside the rising-edge range
#                                                 (3 700–8 300), wave_probe 2026-09-27
#   N_MEM      3                     [ESTIMATE]  H(rest) ≈ 0.03, H(7 500) ≈ 0.77
#   k_memOn    0.004 s⁻¹             [ESTIMATE]  integrator regime: < ~40 % of the
#                                                 pool is written per wave, so the
#                                                 mark does not saturate (a saturated
#                                                 mark is uniform and has no direction)
#   k_memOff   1/180 s⁻¹ (τ = 3 min) [SOURCED, order] Skoge 2014 ~2 min (waves) to
#                                                 ≥ 5 min (reversal); Shi et al.
#                                                 2013 polarity persistence ~2 min
#   D(Mema_mem) 0.01 µm²/s           [ESTIMATE]  spreads √(4Dτ) ≈ 2.7 µm in 3 min
const MEMORY_ON = get(ENV, "DICTY_MEMORY", "0") == "1"
_mem(env, d) = parse(Float64, get(ENV, env, string(d)))
const MEM_V     = round(Int, _mem("DICTY_MEM_V", 5000))
const K_MEM     = _mem("DICTY_MEM_K", 5000.0)
const N_MEM     = _mem("DICTY_MEM_N", 3.0)
const K_MEM_N   = Float32(K_MEM^N_MEM)
const N_MEM_F   = Float32(N_MEM)
const k_memOn   = Float32(_mem("DICTY_MEM_KON", 0.004))
const k_memOff  = Float32(1.0 / _mem("DICTY_MEM_TAU", 180.0))
const D_MEMA    = Float32(_mem("DICTY_MEM_D", 0.01))
rM_1 = Reaction([:Mem_cyto], [:Mema_mem], k_memOn,
                @rate(x -> x[:PIP3_mem]^N_MEM_F / (K_MEM_N + x[:PIP3_mem]^N_MEM_F)))
rM_2 = Reaction([:Mema_mem], [:Mem_cyto], k_memOff, nothing)
const MODULE_MEM_REACTIONS = MEMORY_ON ? [rM_1, rM_2] : Reaction[]
const MEM_DIFFUSION = MEMORY_ON ?
    Dict{Symbol, Float32}(:Mem_cyto => 10.0f0, :Mema_mem => D_MEMA) : Dict{Symbol, Float32}()
const MEM_INITIAL   = MEMORY_ON ?
    Dict{Symbol, Int}(:Mem_cyto => MEM_V, :Mema_mem => 0) : Dict{Symbol, Int}()

const CELL_REACTIONS = [
    MODULE_MEM_REACTIONS...,
    MODULE_0_REACTIONS...,
    MODULE_1_REACTIONS...,
    MODULE_2_REACTIONS...,
    MODULE_3_REACTIONS...,
    MODULE_4_REACTIONS...,   # re-enabled 2026-08-02 — see M4_INITIAL/M4_DIFFUSION note above
    #r5_1, r5_2, r5_3, r5_4, r5_5, r5_6, r5_7, r5_8, r5b_1, r5b_2,
    #r6_1, r6_2, r6_3, r6_4, r6_5, r6_6, r6_7,
    #r7_1, r7_2, r7_3, r7_4, r7_6, r7_7,   # r7_5 omitted: same channel as r6_2
    # [2026-08-09] Module 8 RE-ENABLED, re-topologised and re-calibrated —
    # see the "Module 8.14" block above and calib/CALIBRATION_M8.md.  It is the
    # receptor arm's route into cAMP re-synthesis: PIP3-gated CRAC* → X, i.e. a
    # stimulus to the SAME excitable element the spontaneous pulse uses, so the
    # relayed pulse and the noise-driven pulse are one pulse with one refractory
    # period.  DICTY_MODULE8=0 removes it (restores the pre-2026-08-09 network).
    MODULE_8_REACTIONS...,
    MODULE_910_REACTIONS...,
    MODULE_11_REACTIONS...,
]

# Reactions running in the world (agent-free voxels): Module 11 PdsA system —
# the free-voxel copies of the six dual-copy reactions in MODULE_11_REACTIONS
# (synthesis/secretion/shedding are cell-only and have no world counterpart).
const WORLD_REACTIONS = [
    w11_D4c, w11_D4d, w11_E2, w11_E3f, w11_E3r, w11_E4, w11_2,
]

# ============================================================================
# Diffusion constants — SINGLE SOURCE OF TRUTH for every diffusion coefficient
# used anywhere in the project.
#
# Collected 2026-08-01 from the driver that carried the fullest, most-audited
# set (dicty_simulation.jl's `agent_diff`/`world_diff`, already merging
# M12_DIFFUSION/M3_DIFFUSION and the Module 0.9 X/Y pair) — several other
# scripts (dicty_minimal.jl, mini_simulation.jl) carried older, hand-copied
# subsets of the same table that had silently drifted (missing :X/:Y, so
# Module 0.9's own activator/recovery species diffused at D = 0 there).
# CELL_DIFF/WORLD_DIFF below are that fullest table; every driver now sources
# its diffusion constants from here (filtering to the species it actually
# instantiates where it intentionally runs a reduced network) instead of
# maintaining its own copy of the numbers. See Diffusion_constants.md for the
# literature derivation of every value.
#
# CFL: D × dt / dx² ≤ 0.25; for dt=0.002 s, dx=2 µm: D_max = 500 µm²/s.
# ============================================================================
const CELL_DIFF = merge(
    Dict{Symbol, Float32}(
        # Cytosolic — free diffusion. Generic cytosolic proteins raised
        # 5→10 µm²/s: 5 sat below the eukaryotic in-cell range (10-30 µm²/s);
        # no Dicty-specific cytosolic measurement exists, so 10 is the
        # literature-centred choice (Diffusion_constants.md §4.5/§6). Ca_i/
        # cAMP_i/AA/_mem values below were each separately audited and are
        # already at their justified values.
        :Gbg_cyto        => 10.0f0,
        :Ca_i            => 20.0f0,   # buffered effective D_Ca; no Ca-buffer species in network -> keep (§4.4)
        :cAMP_i          => 30.0f0,   # r11_1 is catalytic (cAMP_i not consumed) -> buffering not explicit -> keep 30-60 range (§4.3)
        :RasGEFR_cyto    => 10.0f0,
        :RasGEFRa_cyto   => .0f0,
        :RcptAdapt_cyto  => 10.0f0,   # FCD adapter (Rev 4 #1)
        :RcptAdapta_cyto => 10.0f0,
        :RasGAP_cyto     => 10.0f0,
        :RasGAPa_cyto    => 10.0f0,
        :PI3K_cyto       => 10.0f0,
        :PTEN_cyto       => 10.0f0,
        :RacGEF_cyto     => 10.0f0,
        :RacGEFa_cyto    => 10.0f0,
        :RacB_GDP        => 10.0f0,
        :RacB_GTP_mem    => 0.1f0,    # membrane-anchored (prenylated CAAX); same as PIP3/PTEN_mem
        :PAKc_cyto       => 10.0f0,
        :PAKca_mem       => 0.1f0,    # membrane-recruited by RacB_GTP_mem
        :PLA2            => 10.0f0,
        :PLA2a           => 10.0f0,
        :AA              => 10.0f0,
        :RasGEFA_cyto    => 10.0f0,
        :RasGEFAa_cyto   => 2.0f0,
        :RasCGAP_cyto    => 10.0f0,
        :RasCGAPa_cyto   => 10.0f0,
        :PKBA_cyto       => 10.0f0,
        :PKBR1_mem       => 0.1f0,   # palmitoylated membrane anchor — slow lateral diffusion
        :PKBR1a_mem      => 0.1f0,
        :CRAC_cyto       => 10.0f0,
        :CRAC_mem        => 0.1f0,   # membrane-associated → slow lateral diffusion, as PKBR1_mem
        :RegA_cyto       => 10.0f0,
        :RegAp_cyto      => 10.0f0,
        :RegAp_i_cyto    => 10.0f0,
        :ERK2_cyto       => 10.0f0,
        :ERK2a_cyto      => 10.0f0,
        :PKA_cyto        => 10.0f0,
        :PKAa_cyto       => 10.0f0,
        :X               => 10.0f0,   # Module 0.9 FHN activator (cytosolic-like)
        :Y               => 10.0f0,   # Module 0.9 FHN recovery variable (cytosolic-like)
        # Fast, so the cell-wide refractory state is well mixed over the 13
        # voxels in ~2.5 s — the gate is evaluated per voxel, so a slowly
        # diffusing Refr would let one voxel escape the block on its own.
        :Refr_cyto       => 10.0f0,   # Module 0.9r relay refractory clock
        :PIP2_mem        => 0.1f0,
        :PIP3_mem        => 0.1f0,
        :PI34P2_mem      => 0.1f0,
        # Module 11 (2026-08-05): PdsA_i_cyto is generic cytosolic (pre-
        # secretion); PdsA_mem is surface-displayed like the other _mem
        # species; cAMP_slow is a pure per-voxel low-pass-filter state (does
        # not itself spread — it exists to read the local cAMP_ext history).
        :PdsA_i_cyto     => 10.0f0,
        :PdsA_mem        => 0.1f0,
        :cAMP_slow       => 0.0f0,
    ),
    # Modules 1/2 own their diffusion constants (M12_DIFFUSION above). This
    # OVERRIDES the generic 10 µm²/s that :Gbg_cyto and :RcptAdapt*_cyto get
    # above, and supplies the receptor/Gα2 values. See calib/CALIBRATION_M12.md §4.
    M12_DIFFUSION,
    # Module 3 owns RasG_GDP_mem/RasG_GTP_mem's diffusion (0.1 µm²/s,
    # prenylated membrane anchor — M3_DIFFUSION above).
    M3_DIFFUSION,
    # Module 4 fills in PI3K_mem/PTEN_mem, missing from the generic dict above
    # until the 2026-08-02 re-enable (M4_DIFFUSION note above r4_1).
    M4_DIFFUSION,
    # Module 8 fills in RasC/TORC2/CRAC*/ACA — all of them missing from the
    # generic dict above (a species absent here silently gets D = 0 in
    # prepare.jl, which for ACA_mem/ACAa_mem was a real gap: the cyclase is
    # membrane-resident and should diffuse like PIP3/PTEN_mem, not not at all).
    M8_DIFFUSION,
    MEM_DIFFUSION,   # empty unless DICTY_MEMORY=1
)

const WORLD_DIFF = Dict{Symbol, Float32}(
    # cAMP_ext: THE PHYSICAL VALUE, 350 µm²/s (Diffusion_constants.md §4.1/§6 Set A =
    # 400 free-solution × 0.87 tortuosity). Cells are NOT excluded volumes for these
    # species (world_shared species below), so the tortuosity haircut must live in D → Set A.
    # dt=0.002 gives a CFL ceiling of 500 µm²/s, so 350 fits (0.175 ≤ 0.25).
    :cAMP_ext       => 350.0f0,
    # [2026-08-05] Module 11 rebuild: PdsA_ext now models the SECRETED form
    # only (the surface-bound population is the explicit, non-diffusing
    # PdsA_mem agent species above) — this resolves the "total pool vs.
    # secreted-only" ambiguity Diffusion_constants.md §4.2 Q2 left open, in
    # favour of secreted-only, so the value moves from the old effective 20
    # (a 25%-free/75%-surface-bound mix) to the secreted-only Set A figure
    # that document already derived: 65 µm²/s. PdsA_PdiA_ext is the same
    # secreted enzyme in complex with the inhibitor, so it gets the same D.
    :PdsA_ext       => 65.0f0,
    :PdsA_PdiA_ext  => 65.0f0,
    :PdiA_ext       => 60.0f0,    # soluble only, no surface-bound form (§4.2/§5) — unchanged
)

"D_CAMP: the extracellular cAMP diffusion constant [µm²/s], read out of WORLD_DIFF."
const D_CAMP = WORLD_DIFF[:cAMP_ext]

# ============================================================================
# Cell initial conditions — SINGLE SOURCE OF TRUTH for every agent-species
# starting count used anywhere in the project.
#
# Collected 2026-08-02 from two_cells.jl (calib/engine_check.jl's own IC:
# pools inactive, RegA at its resting 0.9-phospho set point, Hunger at max —
# a fully competent, starving cell), merged with M12_INITIAL/M3_INITIAL the
# same way CELL_DIFF merges M12_DIFFUSION/M3_DIFFUSION above. Every driver
# now sources its initial state from here (filtering to the species it
# actually instantiates) instead of maintaining its own copy of the numbers.
#
# ⚠ THIS IS THE *LEGACY* IC AND IT IS NOT A REST STATE.  Every entry below is a
# closed-form balance solved for ONE module in isolation, and the whole-cell
# mean field (cell_steady_state.jl) puts its fastest species six orders of
# magnitude off balance at t = 0: `PTEN_cyto`/`PTEN_mem` alone carry a residual
# of 3.0e6 molecules/voxel/s against a pool of 10 000, i.e. they are at the
# WRONG END of an equilibrium that re-establishes itself in ~3 ms, and
# `PIP3_mem` (0.94 /s) and `RasG_GTP_mem` (1.6 /s) are not far behind.  It is
# kept because the reduced-network probes in calib/ (m3_probe, m4_probe,
# m34_step_probe, …) instantiate only a SUBSET of the modules and need the
# matching per-module balances — a whole-cell rest state is the wrong IC for a
# network that is not the whole cell.  `CELL_INITIAL` below selects between this
# and the solved rest state; see `DICTY_IC`.
# ============================================================================
const CELL_INITIAL_LEGACY = merge(
    Dict{Symbol, Int}(
        :Hunger       => H_max_vox,
        :X            => 0,
        :Y            => 0,
        :Refr_cyto    => 0,             # rested cell: the relay block is lifted
        :ERK2_cyto    => ERK2_pool_v,               :ERK2a_cyto   => 0,
        # [FIXED 2026-08-09] 154 → 865/voxel.  865 is the CALIBRATED pool
        # (calib/CALIBRATION.md Iteration 14, and the value calib/engine_check.jl
        # — the script every single-cell number in that document was measured on —
        # has used since); CELL_INITIAL was left at the pre-iteration-14 154, so
        # every driver sourcing its IC from here (two_cells.jl, dicty_simulation.jl)
        # was running the cyclase at 18 % of the calibrated pool and could not
        # reproduce the calibrated export pulse.  See the r0_9f pool note.
        :ACA_mem      => ACA_v,                     :ACAa_mem     => 0,
        :cAMP_i       => 0,
        :RegA_cyto    => round(Int, 0.1 * RegA_pool_v),
        :RegAp_cyto   => round(Int, 0.9 * RegA_pool_v),
        :RegAp_i_cyto => 0,
        :PKA_cyto     => PKA_pool_v,                :PKAa_cyto    => 0,
    ),
    M12_INITIAL,
    M3_INITIAL,
    M4_INITIAL,
    M8_INITIAL,
    M11_INITIAL,   # [2026-08-10] was MISSING entirely — cells used to start with
                   # an empty PdsA secretion chain and took ~2000 s to fill it
    MEM_INITIAL,   # empty unless DICTY_MEMORY=1
)

# ============================================================================
# The HUNGER-DEPENDENT REST STATE, and which IC the drivers get
# ----------------------------------------------------------------------------
# `cell_steady_state.jl` solves the whole cell at once off `CELL_REACTIONS`
# instead of module by module — read its header for why that is the only way to
# get an IC once 3.3b makes Modules 3 and 4 a loop, and for the three things
# (engine units, Hunger as a PARAMETER, the excitable element held quiescent)
# that make the answer agree with the engine.  `analyse/cell_rest_state.jl`
# prints it and checks it; `analyse/m34_steady_state.jl` cross-checks it against
# an independent Catalyst/ODE integration; `calib/m34_rest.jl` measures the same
# fixed point in the SSA and is the acceptance test.
#
# WHAT HUNGER ACTUALLY CHANGES.  Hunger is a slow developmental clock, not a
# signalling variable: 0.8a/0.8b hold it at `H_max_vox` in the engine, and the
# model is entitled to move it underneath the network.  It must therefore be a
# PARAMETER of the rest state, never something a steady-state solve is allowed
# to relax — solve it and the answer is "Hunger = H_max_vox" by construction,
# which is not a statement about anything.  Held fixed and swept, it turns out
# to enter in exactly three places:
#   * `P_dev(H)` — Module 11 PdsA/PdiA synthesis.  This is the ONLY part of the
#     rest state that moves with Hunger (`PdsA_i_cyto`/`PdsA_mem` fall ~9x from
#     H = 100 to H = 0, and the extracellular background with them).
#   * `comp(H)` — the competence gate on 0.9a, i.e. the IGNITION RATE.  This is
#     what Hunger is FOR, and it is invisible to any steady-state calculation:
#     ignition is a rare noise event in the SSA, so the deterministic rest state
#     is the same at H = 0 and H = 100 while the spontaneous pulse interval is
#     not.  Measure that on the engine (`calib/engine_isi.sh`), not here.
#   * `α_H` on 9.1 — scales cAMP synthesis by `ACAa_mem`, which is 0 at rest, so
#     it is silent in the rest state and only shows up inside a pulse.
#
# `DICTY_IC`:
#   auto    (default) the solved rest state if the network HAS one, else the
#           legacy per-module ICs — decided by the POOL-INVERSION check in
#           `cell_rest_state` (does the resting network empty the state the IC
#           declares as rest?), so a network that spends a whole substrate pool
#           at rest can never get an IC written from the inside of that
#           collapse.  Warns, naming the pools, when it falls back.
#   steady  force the solved rest state even if a pool is inverted.
#   legacy  force the pre-2026-08-27 per-module ICs (bit-identical).
# `DICTY_IC_HUNGER` / `DICTY_IC_DENSITY` / `DICTY_IC_AMBIENT` set the operating
# point the rest state is solved AT (default: full competence, the repo's
# 200-cells-on-100² density, zero ambient cAMP).
# ============================================================================
include(joinpath(@__DIR__, "cell_steady_state.jl"))
using .CellRest

"""
    cell_world_background(; hunger, density, ambient_nM)

The extracellular state a resting cell sits in, per voxel — what
`cell_rest_state` holds fixed while it relaxes the cell.  `cAMP_ext` is a free
choice (the ambient the cell is bathed in); the Module 11 enzymes are NOT, they
are the density-matched field `m11_world_steady` already derives, and they scale
with cell density because `PdsA_ext` is produced per CELL and lost per VOXEL.
"""
function cell_world_background(; hunger = H_max_vox,
                                 density = 200 / (100 * 100),
                                 ambient_nM = 0.0)
    Pfree, Ifree, C = m11_world_steady(density; hunger = hunger)
    Dict{Symbol,Float64}(
        :cAMP_ext      => ambient_nM * 1e-3 * Float64(MOLEC_PER_µM),
        :PdsA_ext      => Pfree,
        :PdiA_ext      => Ifree,
        :PdsA_PdiA_ext => C,
    )
end

# ── THE COMMITTED M3/M4 REST STATE (reference values, molecules/voxel) ───────
# `DICTY_IC=auto` SOLVES this at include time from the committed rate constants
# rather than reading a literal table, so it cannot go stale when a rate moves —
# which is the defect this repo has hit repeatedly ("the rest state and the rate
# constants are ONE coupled set"). The numbers below are recorded only so the
# expected state is auditable at a glance; they are NOT the source of truth.
# Measured at k_pip3Ras = 55, K_DISP = 0.22, k_gapBasal = 2.5e-3, brake OFF,
# Hunger = H_max_vox, zero ambient cAMP; solver residual ~5e-14 /s:
#     RasGEFR_cyto  138593    RasGEFRa_cyto      7
#     RasGAP_cyto    42944    RasGAPa_cyto    3256
#     RasG_GDP_mem   38479    RasG_GTP_mem      21   (Ras 0.05 % GTP-loaded)
#     PI3K_cyto     199966    PI3K_mem          34
#     PIP2_mem      769102    PIP3_mem         129   (PIP3 0.02 % of the lipid pool)
#     PTEN_cyto        401    PTEN_mem        9599   (PTEN 96 % on the membrane)
#     PKBA_cyto      10000    PKBAa_mem          0
# Sanity: no conserved pool is inverted, so `analyse/cell_rest_state.jl` installs
# rather than refuses. Before 2026-08-27 the resting cell held RasG-GTP at
# 99.88 % and the whole lipid pool as PIP3, i.e. there was no rest state at all.
const IC_MODE    = get(ENV, "DICTY_IC", "auto")
IC_MODE in ("auto", "steady", "legacy") ||
    error("DICTY_IC must be auto|steady|legacy, got $(repr(IC_MODE))")
const IC_HUNGER  = parse(Float64, get(ENV, "DICTY_IC_HUNGER",  string(H_max_vox)))
const IC_DENSITY = parse(Float64, get(ENV, "DICTY_IC_DENSITY", string(200 / (100 * 100))))
const IC_AMBIENT = parse(Float64, get(ENV, "DICTY_IC_AMBIENT", "0.0"))

const CELL_REST = IC_MODE == "legacy" ? nothing :
    cell_rest_state(CELL_REACTIONS;
                    hunger = IC_HUNGER,
                    world  = cell_world_background(; hunger = IC_HUNGER,
                                                     density = IC_DENSITY,
                                                     ambient_nM = IC_AMBIENT),
                    guess  = CELL_INITIAL_LEGACY)

const CELL_INITIAL = begin
    if CELL_REST === nothing
        CELL_INITIAL_LEGACY
    elseif IC_MODE == "auto" && !isempty(CELL_REST.inverted)
        @warn """DICTY_IC=auto: keeping the LEGACY initial conditions — the network has no rest state to install.
        The resting cell INVERTS $(length(CELL_REST.inverted)) conserved pool(s): $(join(["$(x.species) $(round(100*x.ref/x.total, digits=1))% -> $(round(100*x.rest/x.total, sigdigits=2))% (now $(x.winner))" for x in CELL_REST.inverted], "; ")).
        With no stimulus at all these substrates are spent / these switches sit fully ON, so the solved numbers are a collapse rather than an initial condition.
        Diagnose with `julia analyse/cell_rest_state.jl`. DICTY_IC=steady installs them anyway; DICTY_IC=legacy silences this."""
        CELL_INITIAL_LEGACY
    else
        merge(CELL_INITIAL_LEGACY, CELL_REST.state)
    end
end
