import numpy as np
import networkx as nx
import Math.data as data
from Math.Application import Application
from shutil import copyfile


class Main:
    def __init__(self, graph, R, ST, alpha, variante, kanten_queue, start_queue, ende_queue, y0_queue):
        """
        liest Parameter aus "Math/data.py" ein und startet das Programm (also erzeugt eine "Application" und ruft deren
         "runner" Methode auf).
        :param graph: Gerichteter Graph als Dictionary, falls kein Graph spezifiziert, so wird Graph aus Datei
         "GraphenGEXF/myGraph.gexf" eingelesen
        :param R: Liste aller Startzeitpunkte, indiziert in der Reihenfolge der Spieler
        :param ST: Liste von Tupeln, wobei k-tes Tupel (i,j) beschreibt, dass Spieler k +1 Quelle s_i und Senke t_j
         besitzt
        :param alpha: aus Intervall [0,1]. Legt Gewichtung des Einflusses der Reise- und Wartezeiten auf die Kosten
         fest, beispielsweise:
            alpha = 0: nur die Reisedauer bestimmt Kosten
            alpha = 1: nur die Wartezeit bestimmt Kosten
            alpha = 1/2 (Standardwert): Reisedauer und Wartezeit nehmen gleichen Einfluss auf Kosten
        :param variante: gibt Variante zur Kostenberechnung an; Möglichkeiten: 'A' (standard), 'B', 'C', 'D'
        :param kanten_queue: Liste, die alle Kanten mit virtueller Warteschlange, als Tupel der Form ('v','w'), enthält
        :param start_queue: Liste die zu den Einträgen in 'kanten_queue' die entsprechenden Startzeitpunkte des
         virtuellen Einflusses enthält (i-ter Eintrag in 'start_queue' bezieht sich auf i-ten Eintrag in 'kanten_queue'
        :param ende_queue: Analog zu 'start_queue' ist dies eine Liste, die die Endzeitpunkte des virtuellen Einflusses
         enthält
        :param y0_queue: Liste mit den Einflussgrößen des virtuellen Flusses, indiziert wie 'kanten_queue',
         'start_queue', 'ende_queue'
        """
        if graph is None:
            graph, posit = self.import_gexf_graph()
        else:
            posit = None
        app = Application(graph, R, ST, alpha, posit, variante, kanten_queue, start_queue, ende_queue, y0_queue)
        # starte 'app.runner()', diese Funktion sorgt für den wiederholten Aufruf von 'app.run()'
        app.button_win.after(1000, app.runner)
        app.button_win.mainloop()

    @staticmethod
    def import_gexf_graph():
        """
        Funktion zum Einlesen eines Graphen aus einer .gexf Datei
        :return: graph: gerichteter Graph als Dictionary,
                 posit: Knotenpositionen als passende Koordinaten
        """

        # Funktion aus 'shutil', erzeugt Kopie von 'myGraph.gexf')
        copyfile("GraphenGEXF/myGraph.gexf", "GraphenGEXF/myGraph.gexf" + "~")

        # füge Ausdruck 'defaultedgetype="directed"' zu 'myGraph.gexf' hinzu (falls noch nicht vorhanden), da sonst
        # Kantenorientierung nicht beachtet wird
        source = open( "GraphenGEXF/myGraph.gexf" + "~", "r")
        gexfstring = source.read()
        if "defaultedgetype=\"directed\"" not in gexfstring[:500]:
            destination = open( "GraphenGEXF/myGraph.gexf", "w" )
            print("string:", gexfstring)
            gexfstring1, gexfstring2 = gexfstring.split("graph mode=\"static\"")
            destination.write(gexfstring1 + "graph mode=\"static\" defaultedgetype=\"directed\"" + gexfstring2)
            destination.close()
        source.close()

        G = nx.read_gexf("GraphenGEXF/myGraph.gexf")
        nodes = G.nodes()
        edges = G.edges()
        graph = {}
        posit = {}
        # folgende Werte werden verwendet um Koordinaten passend zu transformieren
        x_min, x_max, y_min, y_max = 1, 1, 1, 1
        for v in nodes:
            graph[nodes[v]['label']] = {}
            x = nodes[v]['viz']['position']['x']  # x- und y-Position von Knoten v
            y = nodes[v]['viz']['position']['y']
            if abs(x) > x_max:  # merke betraglich größte Koordinaten
                x_max = abs(x)
            if abs(y) > y_max:
                y_max = abs(y)
            posit[nodes[v]['label']] = np.array([x, y])
            for m in edges:
                if m[0] == v:
                    cap_trav = edges[m]['label'].split("/")
                    # speichern der Kapazität und der Reisezeit von Kante m, entnommen aus der .gexf-Datei
                    graph[nodes[v]['label']][nodes[m[1]]['label']] = (int(cap_trav[0]), int(cap_trav[1]))
        # x_max, bzw. y_max beschreiben nun die betraglich größte x- bzw. y-Koordinate
        # mit Hilfe dieser werden nun die Koordinaten aller Knotenpositionen auf das Intervall [-1,1]^2 transformiert
        for v in nodes:
            posit[nodes[v]['label']][0] = posit[nodes[v]['label']][0]/x_max
            posit[nodes[v]['label']][1] = posit[nodes[v]['label']][1]/y_max
        return graph, posit


def main():
    try:
        graph = data.graph
    except AttributeError:
        graph = None
    try:
        alpha = data.alpha
    except AttributeError:
        # falls 'alpha' nicht spezifiziert, verwende Standardwert 1/2
        alpha = 1/2
    try:
        variante = data.variante
    except AttributeError:
        # falls 'variante' nicht spezifiziert, verwende Standardwert 'A'
        variante = 'A'
    try:
        kanten_queue = data.kanten_queue
        start_queue = data.start_queue
        ende_queue = data.ende_queue
        y0_queue = data.y0_queue
    except AttributeError:
        kanten_queue = []
        start_queue = []
        ende_queue = []
        y0_queue = []

    Main(graph,data.R,data.ST,alpha,variante,kanten_queue,start_queue,ende_queue,y0_queue)
    return 0


main()
