import Math.cont_data as data
from Math.ContApp import ContApp


class ContMain:
    def __init__(self, graph, u, output):

        self.app = ContApp(graph, u, output)


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
        graph = data.graph
        u = data.u
    except AttributeError:
        raise AttributeError('Daten unvollständig!')
    try:
        output = data.table_output
    except AttributeError:
        output = True

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
