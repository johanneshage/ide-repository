"""
Enthält die für die Erzeugung eines Objekts der Klasse "ContMain" notwendigen Parameter:
    u: Liste aller Einflussraten. Enthält für jeden Knoten eine Liste von 2-Tupeln der Form
       (<Startzeitpunkt>, <Einflussmenge>). <Startzeitpunkt> gibt dabei den Beginn des Einflusses von <Einflussmenge>
       Flusseinheiten in den Knoten mit entsprechendem Index an. Dieser Einfluss bleibt bestehen bis zum
       <Startzeitpunkt> des nächsten 2-Tupels in der selben Liste.
    graph: Gerichteter Graph als Dictionary, falls kein Graph spezifiziert, so wird Graph aus Datei "myGraph.gexf"
           eingelesen
    table_output: Hat Wert 'True' oder 'False' und legt fest, ob eine 'OutputTable' erzeugt wird oder nicht. Falls
                  nicht, so werden die darin enthaltenen Daten stattdessen als Listen ausgegeben (siehe dafür auch
                  "ContMain"). Falls nicht spezifiziert, so wird der Standardwert 'True' verwendet.
"""

"""graph = {
    's1': {'v': (1,1), 'w': (1,2)},
    's2': {'x': (2,3), 'w': (1,1)},
    'v': {'t1': (1,5), 'a': (1,1)},
    'w': {'a': (1,2), 'b': (2,2)},
    'x': {'a': (1,1), 'b': (2,1)},
    'a': {'t1': (2,4), 'b': (1,1)},
    'b': {'a': (1,1), 'v': (1,1)},
    't1': {}
}
"""
graph = {
    's1': {'v': (1,1), 'w': (1,2), 'x': (1,0)},
    's2': {'x': (2,3), 'w': (1,1)},
    'v': {'t1': (1,5), 'a': (1,1)},
    'w': {'a': (1,2)},
    'x': {'a': (1,1)},
    'a': {'t1': (2,4)},
    't1': {}
}

u = [[(0.1, 4.5), (4, 1.5), (5,0)]] + [[(0, 4.2), (1.1, 4.4), (2, 0)]] + 5 * [[]]

"""
graph = {
    's1': {'v': (1.5,2.5), 'W': (1,3), 'u': (3,10), 't1': (1,20)},
    's2': {'W': (5,3.25), 'w': (2,3)},
    'v': {'w': (1,1.25), 'W': (1.75,2.5)},
    'w': {'u': (3.5,3.5), 't1': (2.25,7), 'x': (1,2)},
    'W': {'u': (1,3), 'x': (2,3)},
    'x': {'u': (2,1), 's1': (1,1)},
    'u': {'t1': (2.5,1)},
    't1': {}
}

u = [[(0,11.5), (1.9,2.1), (3.22,1.4), (5.2,0)]] + [[(0.9,12.1), (3.3,1), (6.4,0)]] + 6 * [[]]

"""

"""graph = {
    's1': {'v': (1.1,2.5), 'W': (0.9,3), 'u': (3,10), 't1': (1,25)},
    's2': {'W': (3,3.25), 'w': (2,2.8)},
    's3': {'s1': (3,1), 'W': (2, 1.5), 'x': (1.2, 3), 'a': (2,7)},
    'v': {'w': (1,1.25), 'x': (1.5, 4.5)},
    'w': {'u': (3.5,3.5), 't1': (2.25,9), 'x': (1,2), 'c': (2, 8)},
    'W': {'u': (1,3), 'x': (2,3)},
    'x': {'u': (2,1), 's1': (1,1), 'w': (2, 3.4), 'a': (1.1, 2.1)},
    'u': {'t1': (2.5,2)},
    'a': {'u': (1, 10), 'v': (2, 3.5), 'c':(1.5, 5)},
    'b': {'a': (2,2), 's3': (3, 1.5), 'W': (1,0.5)},
    'c': {'t1': (0.4, 0.9), 'w': (4, 1.5)},
    't1': {}
}

u = [[(0, 16), (0.6, 8), (1.2, 9), (4, 8.1), (7,0)]] + [[(0.2, 8), (1, 5), (1.3, 8), (5.5, 0)]] + \
    [[(1.5, 16), (2, 14), (3, 7), (7,0), (1024,1), (1025,0)]] + 9 * [[]]"""


table_output = False
