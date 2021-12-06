import networkx as nx
import Math.cont_data as data
from Math.ContApp import ContApp
from shutil import copyfile


class ContMain:
    def __init__(self, graph, u, output):

        # liest Graph aus "GraphenGEXF/myGraph.gexf", falls keiner in "Math/cont_data.py" angegeben ist
        if graph is None:
            graph = self.import_gexf_graph()

        self.app = ContApp(graph, u, output)

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
        for v in nodes:
            graph[nodes[v]['label']] = {}
            for m in edges:
                if m[0] == v:
                    cap_trav = edges[m]['label'].split("/")
                    # speichern der Kapazität und der Reisezeit von Kante m, entnommen aus der .gexf-Datei
                    graph[nodes[v]['label']][nodes[m[1]]['label']] = (float(cap_trav[0]), float(cap_trav[1]))
        return graph


def main():
    """
    Je nach Wert der Variable "data.table_output" wird einer der folgenden beiden Outputs erzeugt:
    Fall 1: 'data.table_output' = True: Es wird eine Tabelle ausgegeben, mit der alle Daten eingesehen werden können
    Fall 2: 'data.table_output' = False: Es werden alle benötigten Daten in 'cont_main.app.output' gespeichert. Diese
            enthält eine Liste mit folgenden 8 Datenmengen:
                cont_main.app.output[0]: Terminationszeitpunkt als Float
                cont_main.app.output[1]: Einflussraten pro Kante als Liste von 2-Tupeln (Zeit, Einflussrate)
                cont_main.app.output[2]: Ausflussraten pro Kante als Liste von 2-Tupeln (Zeit, Ausflussrate)
                cont_main.app.output[3]: Liste der Kanten des Graphen, selbe Indizierung wie Einflussraten/Ausflussraten
                cont_main.app.output[4]: Liste, welche für jeden Zeitpunkt eine Liste aller Warteschlangen enthält
                                         (Reihenfolge der Warteschlangen wie in cont_main.app.output[3], für Zeitpunkte
                                         siehe  cont_main.app.output[5])
                cont_main.app.output[5]: Liste aller Zeitpunkte/Phasen
                cont_main.app.output[6]: Knotenlabels zu jedem Zeitpunkt für alle Knoten
                cont_main.app.output[7]: Liste der Knoten
    """
    try:
        u = data.u
    except AttributeError:
        raise AttributeError('Daten unvollständig!')
    try:
        output = data.table_output
    except AttributeError:
        output = True
    try:
        graph = data.graph
    except AttributeError:
        graph = None

    cont_main = ContMain(graph, u, output)
    if not output:
        print("Terminationszeitpunkt")
        print(cont_main.app.output[0])
        print("Einflussraten pro Kante")
        print(cont_main.app.output[1])
        print("Ausflussraten pro Kante")
        print(cont_main.app.output[2])
        print("Reihenfolge Kanten")
        print(cont_main.app.output[3])
        print("Warteschlangen pro Zeitpunkt")
        print(cont_main.app.output[4])
        print("Zeitpunkte")
        print(cont_main.app.output[5])
        print("Knotenlabels pro Zeit für jeden Knoten")
        print(cont_main.app.output[6])
        print("Knotenreihenfolge")
        print(cont_main.app.output[7])
    return 0


main()
