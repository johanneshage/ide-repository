"""
Diskretisierte Version der Instanz aus dem Beweis zu Theorem 4.1 in https://arxiv.org/abs/2007.07808v2
fineness gibt an wie fein die Diskretisierung ist
timehorizon gibt die Laenge der Netzwerkeinflusszeit an
"""

fineness = 2
timehorizon = 20

R = []
# Fuer 10 Zeiteinheiten
for i in range(0,timehorizon*fineness):
    # Einfluss von zwei Paketen pro Zeiteinheit
    for k in 2*fineness*[i]:
        R.append(k)

# Alle Spieler reisen von s1 nach t1
ST = 2*timehorizon*fineness*fineness*[(1,1)]
alpha = 2*timehorizon*fineness*fineness*[1/2]
variante = 'A'

# Keine virtuellen Kanten
kanten_queue = []  # virtuelle Warteschlange
start_queue = []  # Startzeitpunkt des virtuellen Einflusses
ende_queue = []  # Endzeitpunkt des virtuellen Einflusses
y0_queue = []  # Menge des virtuellen Einflusses

# graph, falls kein Graph spezifiziert, so wird Graph aus Datei "myGraph.gexf" eingelesen
graph = {'s1': {'v': (2*fineness,fineness), 'w': (2*fineness,fineness)},
         'v': {'t1': (fineness,fineness)},
         't1': {},
         'w': {'t1': (fineness,2*fineness)}
         }