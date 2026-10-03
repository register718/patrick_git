# Candidate papers for the remaining model assumptions

Save each PDF as `sources/all/<key>.pdf` (key = suggested BibTeX key). Abstracts were read on PubMed; "supports" means what the abstract says, to be confirmed in the full text.
Priority A = fixes a claim the paper currently makes without a source, B = strengthens, C = optional.

| prio | key | paper | DOI | use |
|---|---|---|---|---|
| A | `dinauer1980iv` | Dinauer, Steck, Devreotes 1980, JCB 86:545 — Recovery of the cAMP signaling response after adaptation | 10.1083/jcb.86.2.545 | de-adaptation first order, t½ = 3–4 min → duration of refractory clock (M0) |
| A | `lee2005` | Lee et al. 2005, Mol Biol Cell 16:4572 — TORC2 integrates chemotaxis and signal relay | 10.1091/mbc.e05-04-0342 | TORC2 needed for PKB/PKBR1 phosphorylation **and** ACA activation; Ras regulates TORC2 (M8) |
| A | `lim2001` | Lim, Spiegelman, Weeks 2001, EMBO J 20:4490 — RasC required for ACA and Akt/PKB activation | 10.1093/emboj/20.16.4490 | RasC → ACA and RasC → PKB (M8) |
| A | `chen1997` | Chen, Long, Devreotes 1997, Genes Dev 11:3218 — Pianissimo | 10.1101/gad.11.23.3218 | Pia (TORC2 subunit) **and** CRAC are both essential for ACA: parallel requirement, not "PKB → CRAC" (M8) |
| A | `lilly1995` | Lilly & Devreotes 1995, JCB 129:1659 — CRAC translocation | 10.1083/jcb.129.6.1659 | CRAC is recruited to membranes by Gβγ-dependent sites, independent of PKB (M8, 8.9c) |
| A | `adhikari2020` | Adhikari, Kuburich, Hadwiger 2020, Microbiology 166:129 — MAPK regulation of RegA | 10.1099/mic.0.000868 | Erk2 downregulates RegA, binds RegA (M10, reaction 10.B3) |
| A | `kuburich2019` | Kuburich, Adhikari, Hadwiger 2019, Cell Signal 57:65 — RegA phosphorylation sites | 10.1016/j.cellsig.2019.02.005 | T676 (Erk2 site) phosphorylation **reduces** RegA function (M10) |
| B | `dinauer1980v` | Dinauer, Steck, Devreotes 1980, JCB 86:554 — Adaptation during cAMP stimulation | 10.1083/jcb.86.2.554 | adaptation time course (M1, M0) |
| B | `vanhaastert1983` | Van Haastert & Van der Heijden 1983, JCB 96:347 — excitation/adaptation/deadaptation | 10.1083/jcb.96.2.347 | recovery first order, t½ = 1–2 min (cGMP response) |
| B | `wu1995` | Wu, Valkema, Van Haastert, Devreotes 1995, JCB 129:1667 — Gβ essential | 10.1083/jcb.129.6.1667 | Gβγ links receptor to ACA/GC (M2, M8) |
| B | `bolourani2006` | Bolourani, Spiegelman, Weeks 2006, Mol Biol Cell 17:4543 — RasG vs RasC | 10.1091/mbc.e05-11-1019 | RasC more important for ACA, RasG for chemotaxis (M3/M8) |
| B | `insall1996` | Insall, Borleis, Devreotes 1996, Curr Biol 6:719 — Aimless RasGEF | 10.1016/s0960-9822(09)00453-9 | RasGEF needed for ACA activation and chemotaxis (M8, 8.1) |
| B | `nicholson1998` | Nicholson & Syková 1998, Trends Neurosci 21:207 — extracellular space diffusion | 10.1016/s0166-2236(98)01261-2 | tortuosity reduces effective D (M11, D of cAMP); brain data, not Dictyostelium |
| B | `malchow1972` | Malchow et al. 1972, Eur J Biochem 28:136 — membrane-bound cAMP PDE | 10.1111/j.1432-1033.1972.tb01894.x | surface-bound PDE of responding cells (M11, PdsA_m) |
| B | `barra1980` | Barra, Barrand, Yeh, Coukell 1980 — *pdsA* gene, active PDE production during starvation | 10.1007/bf00272671 | starvation-induced PdsA (M11 gate) |
| C | `rathi1991` | Rathi et al. 1991, Dev Genet 12:82 — prestarvation factor induces gene expression | 10.1002/dvg.1020120115 | density/starvation induction (check that PDE is among the genes) |
| C | `lee1999` | Lee, Parent, Insall, Firtel 1999, Mol Biol Cell 10:2829 — RIP3 | 10.1091/mbc.10.9.2829 | RIP3 (TORC2 component) needed for ACA (M8) |
| C | `iglesias2012b` | Iglesias 2012, Sci Signal 5:pe8 — Adaptation and amplification | 10.1126/scisignal.2002897 | LEGI and amplification review (M3) |
| C | `roos1977` | Roos, Nanjundiah, Malchow, Gerisch 1977 — adenylyl cyclase and differentiation | 10.1016/0045-6039(77)90018-5 | measured ACA activities (M9) |
| B | `xu2017` | Xu et al. 2017, PNAS 114:E10092 — GPCR-controlled membrane recruitment of C2GAP1 | 10.1073/pnas.1703208114 | primary evidence that C2GAP1 activation requires Ras (M3 brake 3.3c$_R$; now cited only via `xu2022`); check whether C2GAP1 acts on RasG |

## Model assumptions for which no supporting paper was found

* **PKB phosphorylates CRAC** (reactions 8.9, 8b.3): no paper shows this. Known: TORC2/PKB and CRAC are *both* needed for ACA (`chen1997`, `lee2005`, `kamimura2008`). Best option: cite them and call the edge an assumption, or restructure as "TORC2$^*$ AND CRAC$^*$".
* **Pool sizes** (RasG, GEF, GAP, PI3K, PIP2 10⁷/cell, PTEN, PKB, RegA, PKA, ERK2, ACA, CRAC, TORC2) and **G:R = 60**: no quantitative copy numbers found; state as design choices (already done).
* **PdsA–PdiA K_d = 1 µM** (measured about 10⁻¹⁰ M, `franke1981`) and **PdsA turnover/shedding rates**: no data; design choices.
* **PIP3 gate exponent n = 8**, **Hunger competence gate**: phenomenological.

## Discrepancies found in papers already in `all/` (paper text already corrected)

* Export/secretion: measured 0.34–0.94 min⁻¹ and intracellular PDE ≈ 1.7 min⁻¹ (`dinauer1980`) vs. model 30 min⁻¹ and much faster RegA.
* PKBR1 hydrophobic-motif phosphorylation peaks at 30–60 s, half-life ≈ 40 s (`cai2010`); model 5–10 s.
* Receptors: ≈ 7 × 10⁴ sites/cell after 4 h (`johnson1991`); model 4 × 10⁴.
