# Modul 1: Abweichungen zwischen Modell, Text und Literatur

Geprüft: [sections/m1.tex](sections/m1.tex), die Tabellen `tab_m1_*.tex` und [dicty_reactions.jl](dicty_reactions.jl) (Zeilen 522–810) gegen die Quellen in `sources/M1/`.
Stand: 2026-10-02.

**Konsistenz Code ↔ Tabellen:** Alle Werte stimmen überein (R_tot 3077/Voxel = 40 001/Zelle, k_on 7,5, k_off 0,45, ε_p 1/8, k_ph 0,02, k_ph0 5e-4, k_dp 5e-3, k_b 0,00578, k_A 0,07, k_A,off 0,02).
Offen ist nur, ob der Text die Werte rechtfertigt. Die Befunde stehen unten nach Gewicht sortiert.

Stationäre und zeitabhängige Werte wurden für die Vier-Zustands-Kette R, RC, R_p, R_pC mit den Konstanten der Tabelle nachgerechnet (Rest = 0 nM cAMP).

---

## A. Abweichungen von der Literatur (Modellverhalten)

### A1. Ruhezustand: zu wenig niedrigaffine Rezeptoren
- **Text:** Die niedrigaffine Klasse wird als phosphorylierter Rezeptor R_p dargestellt (Text, Absatz „Role“).
- **Literatur:** Van Haastert 1984: etwa 40 % hochaffin (K_d 60 nM) und etwa 60 % niedrigaffin (450 nM).
- **Modell:** Ruhepool R_p = k_ph0/(k_ph0 + k_b) = 5e-4/(5e-4 + 0,00578) ≈ **8 %**. Es sind also 92 % hochaffin und nur 8 % niedrigaffin.
- Der Code-Kommentar nennt „~5 %“ (k_ph0/k_b = 0,087 ergibt aber 8 %, nicht 5 %).
- Im Text steht nichts dazu.

### A2. „80 %“ wird nicht erreicht, und die 10 % hochaffin schon gar nicht
- **Text/Tabelle:** k_ph/(k_ph + k_dp) = 80 %, passend zu „Verlust von mehr als 80 % der Bindungsplätze“ (Johnson 1991) und „hochaffiner Anteil fällt von 40 % auf 10 %“ (Van Haastert 1984).
- **Modell (stationär):**

  | cAMP | phosphoryliert | hochaffin (R + RC) |
  |---|---|---|
  | 0 nM | 8 % | 92 % |
  | 100 nM | 69 % | 31 % |
  | 1 µM | **78 %** | **22 %** |

- Die Formel ignoriert Basalpfade und Bindungsgleichgewichte. Der hochaffine Anteil fällt von 92 % auf 22 %, die Literatur sagt 40 % auf 10 %. Nur das Verhältnis (≈4) passt, die Absolutwerte nicht.
- **Johnson 1991:** Wildtyp-Zellen verlieren mindestens 75 %. Die „über 80 %“ stehen für cAR1-überexprimierende Zellen. Der Text nennt 80 % ohne diese Einschränkung.
- **Vaughan 1988** nennt 80–90 % der Rezeptoren im phosphorylierten D-Zustand unter Sättigung und passt als Quelle für die 80 % besser.

### A3. Phosphorylierung etwas zu schnell
- **Literatur:** Vaughan 1988: t½ = 45 s bei sättigendem cAMP. Bei 5 nM: t½ = 1,5 min.
- **Modell:** k_ph = 0,02 s⁻¹ entspricht t½ = 35 s. Die gekoppelte Rechnung ergibt ca. 27 s bei 1 µM und ca. 35 s bei 100 nM, also bis zu 1,7× schneller als gemessen.
- Der Text spricht von „of the order of“, das ist vertretbar. „Literaturverankert“ (Text zu 1.5) wäre zu großzügig.

### A4. Dosis-Wirkung der Phosphorylierung nach rechts verschoben
- **Literatur:** Vaughan 1988: halbmaximal bei **5 nM**, gesättigt bei 100 nM.
- **Modell:**
  - 5 nM → 25 % phosphoryliert
  - 20 nM → 48 %
  - 100 nM → 69 %
  - Sättigung erst bei etwa 1 µM
- Halbmaximal (zwischen 8 % und 78 %) liegt bei etwa 15 nM, um den Faktor 3 verschoben. Keine grobe Abweichung, aber nicht erwähnt.
- Zum Vergleich: Die halbmaximale H→L-Konversion in Van Haastert liegt bei 12,5 nM, passt also.

### A5. Dephosphorylierung nur eine Exponentialfunktion
- **Literatur:** Vaughan 1988: t½ = 2 min bei 22 °C (R-Form-Wiederkehr 2,5 min), aber zweiphasig. Schnelle Phase in 0–3 min, langsamere danach, Reste noch nach 30–45 min.
- **Modell:** Einphasig, t½ = 120 s (k_b). Das stimmt für die Zahl, nicht für die Form.

---

## B. Aussagen im Text, die nicht ganz stimmen

### B1. k_off,p = k_off („gemeinsam in beiden schnellen Klassen“)
- **Text (Tabelle):** „the fast dissociation half-life of ≈ 1 s is common to both fast classes“ (Van Haastert 1984).
- **Quelle:** H: t½ ≈ 1,5 s (k = 0,45 s⁻¹). L: t½ ≈ 0,7 s (k₋L = 1,0 s⁻¹). Biswas Tabelle 1: k₋L = 1,0, k_L = 2,2 µM⁻¹s⁻¹ (K_d = 454 nM).
- **Modell:** Der ganze Affinitätsverlust sitzt auf der On-Rate (k_on,p = 0,938, k_off,p = 0,45). R_p dissoziiert damit 2,2× langsamer als die gemessene L-Klasse.

### B2. Verhältnis „8 = Verhältnis der beiden schnellen Klassen“
- 450 nM / 60 nM = 7,5, nicht 8. Der Text schreibt „equals the ratio“ (Marke M1-ratio75), korrekt wäre „≈“.
- ε_p = 1/8 wurde laut Code auf 480 nM abgestimmt, passend zur gemessenen Welle (Peak 204 nM beim Sender). Das ist ein Modellziel, keine Literaturgröße. Der Text nennt beides, der Fokus liegt aber auf der Literatur.

### B3. Basal-Phosphorylierung 40× langsamer (k_ph0)
- **Tabelle:** Basis „E“ (Schätzung), Text und Tabelle: „40× slower; a basal component exists (Vaughan, Hereld)“.
- **Literatur:** Kein Beleg für 40:1.
  - Vaughan 1988: Einbaurate im D-Zustand ≈ ⅕ der induzierten Rate und ≈ 2× der Basalrate. Daraus folgt induziert : basal ≈ **10 : 1**. Das ist eine Ableitung aus zwei Zahlen und misst den ³²P-Einbau pro Rezeptor, nicht eine Ratenkonstante.
  - Vaughan 1988: Phosphorylierung steigt bei Stimulation um mindestens das 4-Fache.
  - Hereld 1994: Einbau steigt um etwa das 5-Fache, etwa 4 zusätzliche Mol Phosphat pro Mol Rezeptor.
- **Folge:** Bei 10× wäre k_ph0 ≈ 0,002 s⁻¹ und der Ruhepool R_p ≈ 26 % statt 8 % (näher an A1). Das ändert die Kalibrierung und wurde nicht gerechnet.

### B4. k_dp = 5e-3 s⁻¹ (gebundener Rezeptor)
- Wurde nur gewählt, um 80 % zu erreichen (Basis „C“). Gemessen wurde die Dephosphorylierung nur nach cAMP-Entzug (freier Rezeptor, siehe k_b).
- Der Wert hat deshalb keine eigene Literaturgrundlage. Für die 80 % gilt zudem A2.

---

## C. Fehler im Code-Kommentar (nicht im Paper-Text)

### C1. Biswas-Behauptung zu „on-rate only“
- **Kommentar (Zeilen 3–4 im Block 1.3/1.4):** „That is the form Biswas measured … k_PHC = 0.04 vs k_H = 7.5 … 187× slower on-rate“.
- **Biswas Tabelle 1:** Bei PH:C sinkt auch die Off-Rate um den gleichen Faktor (k₋PHC = 2,42·10⁻³ gegenüber 0,45 s⁻¹). K_d(PH:C) = 2,42e-3/0,04 = **60 nM**, die Affinität bleibt gleich.
- Bei L: k_PLC = 5, k₋PLC = 1,5, also K_d = 300 nM statt 454 nM (Affinität steigt sogar).
- Die „3–5-fach“-Abnahme steht im Biswas-Text (Seite 5), aber nicht in den Tabellenwerten.
- Biswas modelliert also Verlangsamung, nicht Schwächung. Das Modell hier schwächt nur die On-Rate (ε_p = 1/8), das ist nicht Biswas.
- Die „t½ = 198 s“ der phosphorylierten Zustände bei Biswas entsteht aus den langsamen Phosphat-Schritten und entspricht im Modell k_dp/k_b, nicht der Bindungsrate.

### C2. Veraltete Kommentare
- „30× less sensitive“ (Block 1.3/1.4) und „(30× weaker on-rate)“ (Kommentar zu r1_3): Seit der Änderung 1/30 → 1/8 am 2026-08-30 stimmt das nicht mehr (K_d,p = 480 nM statt 1800 nM).
- Kommentar zu `k_phos`: Der Block „0.02 → 1.0“ ist durch „1.0 → 0.02“ überholt, steht aber noch ausführlich im Code.

### C3. Rezeptorzahl
- Der Code-Kommentar nennt Saxe 1996 und Pupillo 1992 für 40 000 Rezeptoren. Beide sind weder im Text noch in `sources/M1/`. Der Paper-Text bezeichnet den Wert korrekt als Design-Wahl.

---

## D. Bekannte und im Text benannte Abweichungen (nur zur Vollständigkeit)

- Die Slow-Klasse (4 %) fehlt (benannt).
- Affinitätskonversion H→L hat t½ ≈ 9 s, das Modell t½ ≈ 35 s (benannt in „Limitations“).
- ε_p = 8 liegt über den 3–5 der Literatur (benannt).
- R_tot = 40 001 statt 7·10⁴ (benannt).
- Der Adapter A hat keine Wirkung im Standardnetz (benannt).

---

## E. Gut belegt

- K_d 60 und 450 nM, 4 % langsame Klasse, 40/60 % Verteilung, 9 s Konversion, 10 % Minimum der hochaffinen Klasse (Van Haastert 1984).
- k_on = 7,5 und k_off = 0,45 (Biswas Tabelle 1, Zeile 1/2, aus Van Haastert 1984).
- 7·10⁴ Plätze in 4 h entwickelten Zellen (Johnson 1991) und 7·10⁴ ± 5·10³ (Biswas).
- D_R = 0,024 ± 0,004 µm²/s (Takebayashi, n = 27) und 0,027 µm²/s (Ueda 2001: 2,7·10⁻¹⁰ cm²/s).
- Dephosphorylierung t½ = 2 min nach Entzug (Vaughan 1988).
- Phosphorylierung an Serin-Clustern der C-terminalen Domäne (Vaughan 1988, Hereld 1994), Affinitätsverlust statt Entfernen von der Oberfläche (Caterina 1995).
