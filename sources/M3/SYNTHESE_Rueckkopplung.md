# Modul 3 – PIP3-Rückkopplung und verzögerte Bremse: Synthese und Modellvorschlag

Stand 2026-10-02. Die PDFs liegen in `sources/all/<key>.pdf`.
**Umgesetzt** am 2026-10-02 in `sections/m3.tex` (mit Tabellen), `sections/m8.tex` und `dicty_reactions.jl`. Dabei wurden nur Belege aus Repo-PDFs verwendet und markiert; [P]-Belege stehen nicht im TeX. Abweichung zum Vorschlag unten: τ_I = 10 s (Takeda), die Bremse bleibt global (D = 20). Details in `tools/marks/verify_m3.md`. „S.“ ist die Zeitschriftenseite.
Kennzeichnung je Beleg: **[R]** im Repo-PDF geprüft, **[P]** nur als PMC-Volltext gelesen (PDF fehlt noch, daher Abschnitt statt Seite), **[A]** nur Abstract.

## 1 Kurzfassung

1. **Die PIP3-Rückkopplung gehört an Ras, nicht ans G-Protein.** 3.3g (PIP3 → Gα2βγ-Dissoziation) steht im Widerspruch zu G1–G7; vorgeschlagen ist, es zu streichen und 3.3b (PIP3 → RasG-GDP→GTP) wieder zum Standard zu machen. Das ist die Form des veröffentlichten Ras/PIP3-Modells (Fukushima 2019, Gl. 11).
2. **Der cAMP-induzierte Ras-Puls braucht weder PIP3 noch PKBR1 noch F-Aktin** (R1–R3). Die IFFL 3.1–3.7z bleibt deshalb unverändert der Kern für Aktivierung und Adaptation.
3. **Auch die verzögerte Bremse, die Erregbarkeit und Refraktärzeit erzeugt, braucht weder PIP3 noch TorC2** (E1, E2). Ein Leser PKBR1* ist damit nicht haltbar. Vorschlag: Leser ist RasG-GTP, entsprechend einer Ras-abhängig rekrutierten GAP (C2GAP1-artig, K1–K5). Die Zeitkonstante τ ≈ 10 s ist jetzt durch die gemessene Refraktärzeit verankert (E4, E5).
4. **Die Bremse PKB → Sca1 bleibt erhalten, aber nur für RasC** (Modul 8, 8.2s, K10) und mit eigener Spezies.
5. **Sieben Kalibriertests mit Quellen (Abschnitt 5) ersetzen die bisherigen reinen Designziele.**

## 2 Befunde

### G – G-Protein-Ebene

| ID | Befund | Quelle | Zitat | |
|---|---|---|---|---|
| G1 | Bei Dauerreiz bleibt die G-Protein-Aktivierung bestehen; abgeschaltet wird erst weiter unten. | janetopoulos2001, S. 2408 | „Even though physiological responses subsided, the activation did not decline.“ | R |
| G2 | Bei zwei Stufen steigt die Dissoziation stufenweise und bleibt; Ras antwortet dagegen transient. | xu2022, S. 3 | „The kinetics of G-protein dissociation showed a pattern of two step–like persistent increases“ | R |
| G3 | Auch im Gradienten bleibt das G-Protein über die ganze Oberfläche aktiv. | xu2005, S. 676 | „G-proteins were persistently activated throughout the entire cell surface“ | R |
| G4 | Elzie hat PI3K nicht gemessen. LatA ändert die G-Protein-Aktivierung nicht. | elzie2009, S. 2598 | „…PI3K activity or actin polymerization would have subsided“; „no influence of the actin cytoskeleton on G-protein activation“ | R |
| G5 | Mobilität und Domänenbildung von Gβ hängen nicht von PI3K ab (gemessen ist Mobilität, nicht Dissoziation). | vanhemert2010, S. 2927 | „the F-actin-dependent domain formation was independent of PI3K activity“ | R |
| G6 | Unter LY bleibt die G-Protein-abhängige cGMP-Antwort erhalten. | loovers2006, S. 1503 | „LY294002 did not reduce cAMP-mediated cGMP production“ | R |
| G7 | PI3K-abhängige Ras-Aktivität gibt es auch ohne Gβ. | Sasaki 2007, JCB 178:185, Abschn. „Gβ-independent, PI3K-dependent Ras activation…“ | „…a receptor-mediated, Gβ-dependent pathway and a Gβ-independent, PI3K-dependent pathway“ | P |
| G8 | C2GAP1 bindet Gα2 und dämpft die G-Protein-Aktivierung (c2gapA⁻ zeigt stärkeren FRET-Verlust). Im Vorschlag nicht übernommen, siehe Abschnitt 6. | Xu 2026, Cells 15:819, Abschn. 3.7 | „c2gapA− cells exhibited a notably greater loss of FRET“ | P |

### R – Ras-Puls ohne Folgewege

| ID | Befund | Quelle | Zitat | |
|---|---|---|---|---|
| R1 | Für die Ras-Aktivierung sind PI3K, TorC2, PLA2 und sGC nicht nötig. | kortholt2011, S. 1273 | „The signalling enzymes PI3K, TorC2, PLA2 and sGC are not required for Ras activation“ | R |
| R2 | Sind alle vier Wege gehemmt, bleibt die Ras-Antwort auf gleichmäßiges cAMP gleich, auch mit LatA und auch in *sgc/pla2/pkbR1*-null + LY. | kortholt2011, S. 1273–1274 | „essentially the same response as wild-type cells“; mit LatA „indistinguishable from that of wild-type cells“; „Similar results… sgc/pla2/pkbR1null + 90 µM LY“ | R |
| R3 | LatB ändert weder Dosis-Wirkungs-Kurve noch Antwortzeit. | takeda2012, S. 2 | „latrunculin B treatment did not have a major effect on either the dose-response curve… or the characteristic response time“ | R |
| R4 | RasG wird innerhalb von 5 s aktiviert. | kae2004, S. 604 | „An increase of RasG–GTP was also detected within 5 s“ | R |
| R5 | Im Wildtyp ist Ras-GTP nach 20 s wieder auf Basalniveau (Pan-Ras). | Lee 2010, MBoC 21:1810, Results | „wild-type cells rapidly decrease to basal levels by the 20-s time point“ | P |

### P – PIP3 → Ras

| ID | Befund | Quelle | Zitat | |
|---|---|---|---|---|
| P1 | LY senkt die Ras-Aktivierung, verhindert sie aber nicht. | sasaki2004, S. 511 | „…the Ras activation level was reduced in the presence of LY294002 or LatA“ | R |
| P2 | Der LY-Effekt ist PI3K-spezifisch; PIP3 stabilisiert die Erregbarkeit. | fukushima2019, S. 3, 4 | „…those in pi3k1-5-null cells were relatively unaffected“; „feedback from PIP3 to Ras enhances the excitability“ | R |
| P3 | Im veröffentlichten Modell wirkt PIP3 im GEF-Term auf den Ras-Austausch, also wie 3.3b. | fukushima2019, S. 10 (Gl. 11) | „One term defines the basal activity of Ras and the other defines feedback from PIP3.“ | R |
| P4 | Der Ras-Anteil, der von PIP3 abhängt, ist klein. Weniger PIP3 erhöht PI(3,4)P2. | li2018, S. E9128, E9130 | „largely dependent on decreased PI(3,4)P2. Nevertheless, there is a contribution of feedback from PI(3,4,5)P3.“; bei RG2/RG3 nur „a small increase, <15%“ | R |
| P5 | Die PIP3 → Ras-Rückkopplung wirkt auch ohne Zytoskelett. | Huang 2013, NCB 15:1307, Results Fig. 2i–j | „a positive feedback loop from PIP3 to Ras activity that is independent of cytoskeletal activities“ | P |
| P6 | Es gibt eine Ras-eigene, PIP3-unabhängige positive Rückkopplung über RasGEFX/B (spontane Erregung). | iwamoto2025, S. 10 | „further suggesting positive feedback between RasGEFX/B and Ras-GTP“ | R |

### E – Erregbarkeit und verzögerte Hemmung

| ID | Befund | Quelle | Zitat | |
|---|---|---|---|---|
| E1 | Ras-Wellen gibt es auch, wenn PIP3-, TorC2-, PLA2- und sGC-Weg alle gehemmt sind. | fukushima2019, S. 2 | „Even in the cells with all major four pathways (the PIP3, TorC2, PLA2 and sGC pathways) inhibited, Ras waves were still observed“ | R |
| E2 | Ras ist auch ohne Folgewege und ohne Aktin erregbar. | fukushima2019, S. 7 | „Ras-GTP constitutes an excitable network in the absence of downstream parallel PI3K, TorC2, GC and PLA2 pathways and of a functional actin cytoskeleton“ | R |
| E3 | Erregbarkeit setzt eine verzögerte negative Rückkopplung voraus. | fukushima2019, S. 1 | „delayed negative-feedback for the temporal decline of the response followed by the refractory period“ | R |
| E4 | Refraktärzeit des RBD-Signals (LatA, cAMP-Doppelpuls): absolut 8,8 s, Erholungs-Halbwertszeit 7,4 s. | Huang 2013, Fig. 2h | „absolute refractory period (8.8 ± 1.1 seconds) and recovery half-life (7.38 ± 1.74 seconds)“ | P |
| E5 | Mit LimE (F-Aktin) gemessen: absolut <12 s, Halbwertszeit ~7 s, nach 45 s voll erholt. | artemenko2016, S. E7503 | „recovered with a halftime of ∼7 s, as observed for chemoattractant-induced stimulation“ | R |
| E6 | Das Molekül der verzögerten Hemmung ist unbekannt. | li2018, S. E9134; Huang 2013, Discussion | „hypothetical F and R states“; „Although the molecular basis is not known“ | R/P |

### K – Kandidaten für die Bremse

| ID | Befund | Quelle | Zitat | |
|---|---|---|---|---|
| K1 | C2GAP1 wird nur rekrutiert, wenn Ras an der Membran ist; das entspricht einer negativen Rückkopplung (NFBLB). | xu2022, S. 2 | „membrane translocation and activation of C2GAP1 requires Ras proteins on the membrane… indicating the involvement of an NFBLB mechanism“ | R |
| K2 | Eine Ras-GTP-aktivierte GAP (NFBLB) passt genauso gut zu den Ras-Daten wie die IFFL; der Parameterraum wurde aber nicht abgesucht. | xu2022, S. 4, 13 | „we identified a set of parameters for each model“; „we did not determine the space of possible parameter values“ | R |
| K3 | Takeda hat nur die *integrale* Variante verworfen (Abbau nullter Ordnung). Der Puffer muss langsam sein, sonst schwingt das System. | takeda2012_supp, S. 3, 9; takeda2012, S. 3 | `dGAP/dt = kGAP·Ras-GTP − k−GAP`; „needs to be slow to avoid oscillations“ | R |
| K4 | Ohne C2GAP1 ist die basale Ras-Aktivität erhöht. | xu2021c2gap, S. 4 | „c2gapA− cells displayed an enhanced basal Ras activity“ | R |
| K5 | C2GAP1 kommt nach Ras: Membran bei ~15 s, Rückzug bei ~30 s, zweite Phase bei ~1 min. Ras hat seinen Peak ≤5 s (R4). | Xu 2026, Abschn. 3.7 | „initial plasma membrane translocation of C2GAP1 at ~15 s, followed by withdrawal at ~30 s, a second translocation at ~1 min“ | P |
| K6 | NfaA ist eine RasG-spezifische GAP: Ohne sie ist RasG verzögert und verlängert aktiv, RasC/RasD/Rap1 bleiben unverändert. | zhang2008, S. 1587 | „…RasG… delayed and extended considerably in nfaA− cells… RasD, Rap1, and RasC… is unaffected“ | R |
| K7 | PKBs sind nur Kandidaten (RBD). PKBR1 hängt kaum von PIP3 ab. | miao2017, S. 331; kamimura2008, S. 1035 | „could act as a delayed negative feedback regulator“; „PKBR1… relatively independent of PIP3“ | R |
| K8 | PKA hemmt RasG und Rap1, ist aber langsam (über ACA/cAMP). | scavello2017, S. 1550 | „RasG and Rap1 display elevated basal activity levels, and elevated and extended cAMP-induced activation, in pkaC null cells“ | R |
| K9 | Ohne Myosin II ist Ras verlängert aktiv. Da der Ras-Puls ohne F-Aktin normal bleibt (R2, R3), ist Myosin nicht der Kern. | Lee 2010, Results | „There was an insignificant extension of RasC activation, but… we cannot determine if… RasG“ | P |
| K10 | Die Sca1-Rückkopplung betrifft nur RasC. | charest2010, S. 740 | „the Sca1 complex regulates RasC activity while RasG activation is unaffected“ | R |
| K11 | Im veröffentlichten Modell wird die RasGAP durch Gβγ aktiviert, also wie 3.4e. | meierschellersheim2006, PDF-S. 5 | „‘RasGAP’ … deactivates Ras after activation by Gbc“ | R |

### L – Ruhezustand und Einrasten (Latch)

| ID | Befund | Quelle | Zitat | |
|---|---|---|---|---|
| L1 | Ein dauerhaft aktiver Zustand entsteht nur durch Doppelstörungen und ist reversibel. | Edwards 2018, PNAS 115:E3722, Introduction | „despite the assessment of numerous single perturbations, cells were not seen to reach a persistently activated state“ | P |
| L2 | In pten⁻ endet das Ras-Signal nach ~20 s, während PIP3 >1 min bleibt. | Sasaki 2007, Results | „PHcrac was retained at the plasma membrane for >1 min, whereas RBD returned to the cytosol in ∼20 s“ | P |
| L3 | Nach Ende des Reizes ist das Netz in <1 min wieder im Ruhezustand. | xu2022, S. 3 | „returned to prestimulus states in less than 1 min“ | R |
| L4 | Bistabil ist PIP3/PTEN, nicht Ras. | matsuoka2018, S. 4; fukushima2019, S. 6 | „suggesting the characteristics of bistability“; „the PIP3–PTEN bistable state is subordinate to the Ras-GTP excitable state“ | R |
| L5 | Modellbefund: Ein PKB-Leser ist in Ruhe ≈0, deshalb hält nichts gegen die Rückkopplung. | `dicty_reactions.jl` L1285–1287 | „PKB (≈0 at rest) nothing opposes 3.3b and the rest-state solve lands on the LATCHED branch“ | – |

### S – Was misst die RBD-Sonde?
- kae2004, S. 602 [R]: Im Pull-down bindet Raf1-RBD an RasG, nicht an RasC.
- iwamoto2025, S. 2 [R]: In lebenden Zellen hängt das RBD-Signal vor allem von RasG und Rap1 ab, mit kleinem RasC-Anteil („these three small GTPase are the primary factors“).
- scavello2017, S. 1550 [R]: RasG gilt im Aggregationsstadium als das vorherrschende Ras.

Daraus folgt: RBD-Daten dürfen für RasG verwendet werden, wenn man den Rap1-Anteil erwähnt.

## 3 Folgerungen

- **F1 – 3.3g streichen.** Nach G1–G7 folgt das G-Protein dem Reiz und nicht PIP3. Die PIP3-abhängige Ras-Aktivität bleibt auch ohne Gβ erhalten.
- **F2 – PIP3 wirkt auf Ras (3.3b).** Diese Form hat das veröffentlichte Modell (P3) und sie ist mit P1, P2 und P5 vereinbar. Welches Molekül angegriffen wird (GEF oder GAP, vgl. P4), ist nicht gemessen; im Text als offen benennen. Für den Netto-Austausch sind beide Formen dynamisch gleichwertig.
- **F3 – Der Gain der PIP3-Schleife muss klein sein.** Der cAMP-Puls ändert sich ohne PIP3 kaum (R2) oder nur abgeschwächt (P1). Der PIP3-Anteil ist klein (P4). → Test T2.
- **F4 – Die Bremse liest RasG-GTP, nicht PKBR1\* oder PIP3.** Erregbarkeit gibt es ohne PIP3, TorC2 und F-Aktin (E1, E2); der cAMP-Puls ist auch in *pkbR1*-null + LY normal (R2). Belegt ist eine Ras-abhängige GAP-Rekrutierung (K1, K4, K5). Die Topologie ist NFBLB (K2). Takedas verworfene integrale Form (K3) liegt nicht vor, solange B* mit erster Ordnung zerfällt und die Adaptation weiter aus der IFFL kommt. Dynamisches Argument (Modell, L5): Ein RasG-GTP-Leser bemerkt den eingerasteten Zustand direkt, ein PKB-Leser nicht.
- **F5 – τ ist jetzt gemessen verankert.** Annahme: Die Erholung ist durch den Zerfall von B* bestimmt. Dann gilt τ = t½/ln 2 = 7,38 s / 0,693 = 10,6 s (E4); E5 gibt ~10 s. Also bleibt `k_B,off` = 0,1 s⁻¹.
- **F6 – Sca1 nur für RasC.** K10 und R2 sprechen dafür; 8.2s bekommt eine eigene Spezies Sca1*, die von PKBR1* gelesen wird.
- **F7 – Kein Latch bei Einzelstörungen.** Nach L1–L3 darf das Modell bei *pten*⁻ und LY nicht einrasten. → Test T3.
- **F8 – Reichweite der Bremse.** C2GAP1 sammelt sich an der Front (Xu 2026 [P], Abschn. 4: „preferential C2GAP1 accumulation at the front“). B* sollte also membranständig und lokal sein (`DICTY_M3_DBRAKE=0.4` statt 20). Das widerspricht der derzeitigen globalen Einstellung und muss getestet werden.

## 4 Vorgeschlagene Reaktionen

| ID | Reaktion | Rate | Status | Belege |
|---|---|---|---|---|
| 3.1/3.2/3.3 | Gβγ aktiviert GEF_R, GEF_R* macht RasG-GTP | unverändert | belegt | kae2007, kortholt2013, takeda2012 |
| 3.4e/3.4b/3.5/3.7z | Gβγ- und basal aktivierte GAP, zero-order NF1 | unverändert | belegt bzw. Modellform | takeda2012, meierschellersheim2006 (K11), zhang2008 (K6), xu2021c2gap (K4), nakajima2014 |
| ~~3.3g~~ | ~~Gαβγ + PIP3 → Gα2 + Gβγ + PIP3~~ | – | **streichen** | G1–G7 |
| 3.3b | RasG^GDP + PIP3 → RasG^GTP + PIP3 | `k_PIP3,R [PIP3][RasG^GDP]`, optional sättigend `V·[PIP3]/(K+[PIP3])·[RasG^GDP]` | Modellform aus fukushima2019 Gl. 11; Angriffspunkt offen | P1–P5 |
| 3.4d | PIP3 inaktiviert GAP* | – | bleibt aus (nur <15 % Effekt, P4) | li2018 |
| 3.3c_R (neu) | RasG^GTP + B → B* + RasG^GTP | `k_B,on [RasG^GTP][B]` | Kandidatenmechanismus (C2GAP1-artig); Molekül nicht gesichert | K1, K2, K4, K5, E1, E2 |
| 3.3d | B* → B | `k_B,off` = 0,1 s⁻¹ (τ ≈ 10 s) | τ aus Messung (F5) | E4, E5 |
| 3.3e | B* + RasG^GTP → RasG^GDP + B* | `k_B,hyd [B*][RasG^GTP]` | wie bisher; Stärke über T4/T7 kalibrieren | K1 |
| S.1/S.2 (M8, neu) | PKBR1* + Sca1 → Sca1* + PKBR1*; Sca1* → Sca1 | wie bisher 3.3c′/3.3d | belegt für RasC | charest2010, K10 |
| 8.2s | Sca1* + RasGEFA* → RasGEFA + Sca1* | `k_sca1Fb` | belegt für RasC | charest2010 |
| optional O1 | Ras-eigene positive Rückkopplung (Ras-GTP → GEF) | – | nur nötig, wenn spontane Ras-Wellen ohne PIP3 nachgebildet werden sollen | P6, E1 |

Umsetzung im Code:
- `DICTY_M3_PIP3PATH=ras` gibt es schon.
- Neu nötig sind der Leser `DICTY_M3_BRAKEREAD=ras` (Reaktion `RasG_GTP_mem + RasBrake_cyto → RasBrakea_mem + RasG_GTP_mem`) und eine eigene Sca1-Spezies für 8.2s, damit RasBrakea_mem nicht mehr beide Zweige bedient.
- `DICTY_M3_DBRAKE=0.4` gibt es schon (lokale Bremse).

## 5 Kalibriertests (ersetzen Designziele)

| Test | Protokoll im Modell | Sollwert | Quelle |
|---|---|---|---|
| T1 | gleichmäßiger cAMP-Schritt | RasG-Peak ≤5 s; Rückkehr zum Basalwert bis ~20 s; Plateau nahezu adaptiert, kleiner Rest erlaubt | R4, R5; takeda2012 S. 2; xu2022 S. 3 |
| T2 | „LY“: PI3K bzw. 3.3b aus | cAMP-Puls bleibt erhalten, höchstens abgeschwächt | R2 (gleich) vs. P1 (reduziert) → als Spanne |
| T3 | „pten⁻“: PTEN_v → 0 | kein Einrasten; RasG endet nach ~20 s, PIP3 bleibt länger | L2, L1 |
| T4 | Doppelpuls 2 s, Abstand Δt | zweite Antwort fehlt bei Δt < ~9 s; Erholung mit t½ ≈ 7 s; voll erholt bei 45 s | E4, E5 |
| T5 | Reiz entfernen | Rückkehr in den Ruhezustand in <1 min | L3 |
| T6 | G-Protein-Dissoziation bei zwei Stufen | stufenförmig und anhaltend, unabhängig von PIP3 | G1, G2 |
| T7 | Einzelschritt, Abklingphase | keine Oszillation (Puffer langsam genug) | K3 |

## 6 Was offen und nicht belegt bleibt

- **Molekül der RasG-Bremse:** C2GAP1 ist der am besten belegte Kandidat. Ob es an RasG bindet, steht in xu2017 – das PDF fehlt.
- **Angriffspunkt der PIP3-Rückkopplung:** GEF oder GAP (P3 gegen P4) ist nicht gemessen.
- **Widersprüche:**
  - P1 (Sasaki 2004: LY senkt Ras) gegen R2 (Kortholt 2011: gleich).
  - Miao 2017 (LY+PP242 erhöht RBD; Suppl. Fig. 6 nicht gesehen) gegen P5 (Huang 2013).
  - Cheng & Othmer 2016 („no known direct feedback“) gegen P5.
  - Takeda 2012 (nahezu perfekte Adaptation) gegen xu2022 (unvollständige Adaptation).
- **G8 (C2GAP1–Gα2) nicht übernommen:** Die Bindung besteht schon vor dem Reiz, ihre Dynamik ist nicht quantifiziert, und die Quelle liegt nur als PMC-Text vor.

## 7 Neue Claims für `tools/marks/m3.json` (Vorschlag)

| Id | Aussage | Paper, Seite |
|---|---|---|
| ras_no_downstream | Der cAMP-Ras-Puls ist ohne PI3K, TorC2/PKBR1, PLA2, sGC und F-Aktin unverändert. | kortholt2011 S. 1273–1274 |
| ras_excit_intrinsic | Ras ist ohne die vier Folgewege und ohne Aktin erregbar. | fukushima2019 S. 2, 7 |
| pip3_gef_form | Im veröffentlichten Modell geht PIP3 in den Ras-GEF-Term ein. | fukushima2019 S. 10 |
| pip3_small | Der PIP3-Anteil an der Ras-Aktivität ist klein. | li2018 S. E9128 |
| c2gap_ras_dep | Die C2GAP1-Rekrutierung braucht Ras (NFBLB). | xu2022 S. 2 |
| nfblb_fits | Ein NFBLB-Modell beschreibt die Ras-Daten ebenfalls. | xu2022 S. 4 |
| integral_def | Verworfen wurde die integrale Form (Abbau nullter Ordnung). | takeda2012_supp S. 3 |
| refractory | Refraktärzeit: Erholungs-Halbwertszeit ~7 s. | artemenko2016 S. E7503 |
| g_steps | Die G-Protein-Dissoziation ist stufenförmig und anhaltend. | xu2022 S. 3 |
| reset_1min | Nach Ende des Reizes ist das Netz in <1 min im Ruhezustand. | xu2022 S. 3 |
| rbd_rasg_rap1 | RBD meldet in Zellen RasG und Rap1. | iwamoto2025 S. 2 |

## 8 Noch fehlende PDFs (für Seitenzahlen bei [P])

Sortiert nach Priorität:
- `xu2017` (PNAS 114:E10092, C2GAP1)
- `huang2013` (NCB 15:1307, Refraktärzeit)
- `sasaki2007` (JCB 178:185)
- `edwards2018` (PNAS 115:E3722)
- `lee2010` (MBoC 21:1810)
- `xu2026` (Cells 15:819)
- `xu2007` (JCB 178:141)
- optional: `banerjee2025` (JCS 138:jcs263634), `cheng2016` (PLoS CB 12:e1004900)
