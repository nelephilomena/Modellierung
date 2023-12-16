# Modellierung
Skript und Funktionen zur Modellierung eines Schwarmverhaltens mit der Berücksichtigung einer Orientierungspräferenz der simulierten Agenten sowie einer externen Störung.

## Funktionen
CircStat2012a
- circ_ang2rad           -> Winkel von Winkelmaß in Bogenmaß
- circ_rad2ang           -> Winkel von Bogenmaß in Winkelmaß
- circ_mean              -> Zirkulärer Mittelwert in Bogenmaß
- circ_r                 -> Durchschnittliche Resultantlänge der Winkel in Bogenmaß
- circ_rtest             -> Rayleigh Test der Winkel in Bogenmaß
- circ_dist              -> Paarweise Distanz von zwei Winkeln in Bogenmaß
- circ_std               -> Zirkuläre Standardabweichung in Bogenmaß
- circ_vmrnd             -> Simuliert n Zufallswinkel im Bogenmaß aus einer von Mises-Verteilung, mit bevorzugter 
                            Richtung und Konzentrationsparameter </p>
P. Berens, CircStat: A Matlab Toolbox for Circular Statistics, Journal of Statistical Software, Volume 31, Issue 10, 2009
http://www.jstatsoft.org/v31/i10

Eigene Funktionen
- read_tab_group         -> Einlesen der Tabellen mit den Positionen (x,y) der Fischgruppen
- read_tab_single_phi    -> Einlesen der Tabelle mit allen Positionswinkel von allen einzelnen Fischen
- read_tab_single_xy     -> Einlesen der Tabellen mit den Positionen (x,y) der einzelnen Fischen
- modell_schwarm         -> Modellierung von Trajektorien eines Fisches bzw. einer Fischgruppe für gegebende Parameter: lambda, eta, N, time_sim, psi 

## Skript
- main_modellierung        -> Modelliert für variierende Parameter eta und lambda die Trajektorien von einzelnen Fischen und 
                            einer Fischgruppe sowie der Berechnung von 
                            statistischen Größen
- Randbedingung            -> Modellierung von Trajektorien eines einzelnen Fisches (N = 1) mit unterschiedlichen Startpositionen und Anfangsschwimmrichtungen, um den Einfluss des Beckenrandes zu simulieren (beschränkte Randbedingungen)

## Daten
### Dodos
- 'Dodos_angles_all_Gabi.xlsx'
  - Tabelle mit den Winkeln in den Spalten pro Beobachtungszeitpunkt 
  - Für jeden Versuch und für alle einzeln geteste Fische der Art 
- '4.D.G2.3_47B.3D11D12D13D14.3.xlsx'
  - 215 einzelne Tabellen mit den x- und y-Werten der Positionen der Fische 
  für jeden Beobachtungszeitpunkt
  - Jeweils die ersten x- und y-Werte geben das Zentrum des Beckens an
  - Nicht vollständig (nicht für jeden Versuch von einem Fisch vorhanden)

### Fischgruppen
- 49 Tabellen mit jeweils den x-und y-Werten der Positionen der Fischgruppe
  - Es lassen sich nicht die Datenpunkte exakt einem Fisch zuordnen 
     aufgrund der Zeit zwischen den Beobachtungszeitpunkten
  - Jeweils die ersten x- und y-Werte geben das Zentrum des Beckens an



