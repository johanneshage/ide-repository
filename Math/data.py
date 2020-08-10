"""
Enthält die für die Erzeugung eines Objekts der Klasse "Main" notwendigen Parameter:
    R: Liste aller Startzeitpunkte, indiziert in der Reihenfolge der Spieler
    ST: Liste von Tupeln, wobei k-tes Tupel (i,j) beschreibt, dass Spieler k +1 Quelle s_i und Senke t_j besitzt
    alpha: aus Intervall [0,1]. Legt Gewichtung des Einflusses der Reise- und Wartezeiten auf die Kosten
         fest, beispielsweise:
            alpha = 0: nur die Reisedauer bestimmt Kosten
            alpha = 1: nur die Wartezeit bestimmt Kosten
            alpha = 1/2 (Standardwert): Reisedauer und Wartezeit nehmen gleichen Einfluss auf Kosten
    variante: gibt Variante zur Kostenberechnung an; Möglichkeiten: 'A', 'B', 'C', 'D'
    kanten_queue: Liste, die alle Kanten mit virtueller Warteschlange, als Tupel der Form ('v','w'), enthält
    start_queue: Liste die zu den Einträgen in 'kanten_queue' die entsprechenden Startzeitpunkte des virtuellen Ein-
    flusses enthält (i-ter Eintrag in 'start_queue' bezieht sich auf i-ten Eintrag in 'kanten_queue'
    ende_queue: Analog zu 'start_queue' ist dies eine Liste, die die Endzeitpunkte des virtuellen Einflusses enthält
    y0_queue: Liste mit den Einflussgrößen des virtuellen Flusses, indiziert wie 'kanten_queue', 'start_queue',
                                                                                 'ende_queue'
    graph: Gerichteter Graph als Dictionary, falls kein Graph spezifiziert, so wird Graph aus Datei "myGraph.gexf"
    eingelesen
"""


R = 3*[0]+2*[1] + 2*[2] +2*[3] + [4]
# i-tes Tupel enthält Nummer des Start- und des Endknotens von Spieler i
ST = [(1, 1), (1, 2)] + 3*[(2, 1)] + [(2, 2), (1, 1), (1, 2), (1, 1), (2, 1)]
alpha = 10*[1/2]
variante = 'A'
kanten_queue = [('w', 't1'), ('w', 'v')]  # virtuelle Warteschlange
start_queue = [2, 2]  # Startzeitpunkt des virtuellen Einflusses
ende_queue = [3, 2]  # Endzeitpunkt des virtuellen Einflusses
y0_queue = [4, 4]  # Menge des virtuellen Einflusses

# graph, falls kein Graph spezifiziert, so wird Graph aus Datei "myGraph.gexf" eingelesen
graph = {'s1': {'v': (1,2), 'w': (1,1)},
         's2': {'u': (3,5)},
         'v': {'t1': (2,7)},
         'w': {'v': (2,1), 't1': (2,6)},
         'u': {'v': (1,1), 'w': (1,1), 'x': (2,1)},
         'x': {'t1': (1,3), 't2': (1,4)},
         't1': {'t2': (2,1)},
         't2': {}
         }