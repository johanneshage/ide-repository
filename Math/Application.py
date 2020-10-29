from Graphics.ButtonWindowFrame import ButtonWindowFrame
from Graphics.AbfrageVirtuelleSpieler import AbfrageVirtuelleSpieler
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from collections import deque
import sys


class Application:
    """
    Objekt der Klasse "Application" stellt eine Anwendung des Hauptprogramms dar
    """

    def __init__(self, G, R, ST, alpha, posit=None, variante='A', kanten_queue=[], start_queue=[],
                 ende_queue=[], y0_queue=[]):
        """
        Initialisiert Graph und alle Variablen, sowie Kantenkosten und erzeugt Plot zum Zeitpunkt 0 (vor eventuellem
        Hinzufügen weiterer virtueller Spieler durch "AbfrageVirtuelleSpieler.add_parameter()" ). Alle Eingabeparameter
        außer "posit" und "variante" werden in "data.py" spezifiziert und haben folgende Funktion:
        :param G: Gerichteter Graph als Dictionary
        :param R: Liste aller Startzeitpunkte, indiziert in der Reihenfolge der Spieler
        :param ST: Liste von Tupeln, wobei k-tes Tupel (i,j) beschreibt, dass Spieler k +1 Quelle s_i und Senke t_j
         besitzt
        :param alpha: aus Intervall [0,1]. Legt Gewichtung des Einflusses der Reise- und Wartezeiten auf die Kosten
         fest, beispielsweise:
            alpha = 0: nur die Reisedauer bestimmt Kosten
            alpha = 1: nur die Wartezeit bestimmt Kosten
            alpha = 1/2 (Standardwert): Reisedauer und Wartezeit nehmen gleichen Einfluss auf Kosten
        :param posit: Knotenpositionen, werden in "Main.py" zusammen mit dem Graphen aus Datei "myGraph.gexf"
         eingelesen, falls in "data.py" kein Graph angegeben, Standardwert None
        :param variante: gibt Variante zur Kostenberechnung an; Möglichkeiten: 'A', 'B', 'C', 'D'. Standardwert 'A'.
        :param kanten_queue: Liste, die alle Kanten mit virtueller Warteschlange, als Tupel der Form ('v','w'), enthält
        :param start_queue: Liste die zu den Einträgen in "kanten_queue" die entsprechenden Startzeitpunkte des
         virtuellen Einflusses enthält (i-ter Eintrag in "start_queue" bezieht sich auf i-ten Eintrag in "kanten_queue"
        :param ende_queue: Analog zu 'start_queue' ist dies eine Liste, die die Endzeitpunkte des virtuellen Einflusses
         enthält
        :param y0_queue: Liste mit den Einflussgrößen des virtuellen Flusses, indiziert wie "kanten_queue",
         "start_queue", "ende_queue"
        """
        self.kanten_queue = kanten_queue  # Liste der Kanten mit virtuellen Spielern
        self.start_queue = start_queue  # Liste der Startzeitpunkte der virtuellen Zuflüsse
        self.ende_queue = ende_queue  # Liste der Endzeitpunkte der virtuellen Zuflüsse
        self.y0_queue = y0_queue  # Liste der virtuellen Zuflussmengen
        self.button_win = ButtonWindowFrame(self)  # erzeuge Abfragefenster
        self.maxiter = 150
        self.E = []  # Kanten
        self.nu = []  # Kapazitaeten
        self.r = []  # Reisezeiten
        self.c = []  # Kosten
        self.label = []  # Knotenlabels
        self.num = len(R)  # Anz. Spieler
        self.I = range(self.num)  # Liste der Spieler
        self.z = np.negative(np.ones(self.num))  # Liste der Ankunftszeiten
        self.V = list(G.keys())  # Liste der Knoten
        self.n = len(self.V)  # Anzahl Knoten
        self.fp = []  # f^+
        self.fm = []  # f^-
        self.ST = ST  # Zuordnung Start- und Zielknoten zu Spieler
        self.alpha = alpha
        self.R = R  # Startzeitpunkte
        init_m = []
        init_p = []
        for i in self.I:
            init_p.append(None)
            init_m.append(None)
        self.fm.append(init_m)
        self.fp.append(init_p)
        self.variante = variante
        self.G = G  # Graph
        self.ankunft = 0  # zählt angekommene Spieler
        self.deques = []  # Liste aller deques (= Warteschlangen)
        # Liste "leaveQueue" enthält für jeden Zeitpunkt für jede Kante eine Liste der Spieler, die zum gegebenen
        # Zeitpunkt die Warteschlange dieser Kante verlassen
        self.leaveQueue = [[]]
        self.items = G.items()
        self.keys = G.keys()

        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        fig_manager = plt.get_current_fig_manager()
        # fig_manager.full_screen_toggle()  # setze Plot auf Vollbild
        fig_manager.window.state('zoomed')  # maximiere Plot
        self.fig.show()  # zeige Plot an
        plt.ylim(-1.5, 1.5)  # setzte Höhe des Plots, damit Rechtecke in jedem Plot gleich groß
        plt.xlim(-1.25, 1.25)
        if posit is None:
            self.posit = nx.shell_layout(G)  # Layout für Knotenpositionen
        else:
            self.posit = posit  # falls graph aus gexf-datei importiert wurde, verwende auch dessen Layout

        # beendet gesamtes Programm bei Klicken des 'x' - Buttons
        fig_manager.window.protocol('WM_DELETE_WINDOW', sys.exit)

        colorMap = []
        # Farben werden für unterschiedliche Zielknoten verwendet (gibt es mehr Zielknoten als "len(self.colors)", so
        # werden die Farben mehrfach verwendet)
        self.colors = ["limegreen", "orange", "magenta", "turquoise", "gold", "darkviolet"]
        graph = nx.DiGraph()  # erzeuge networkx-Graph
        for v in self.V:
            graph.add_node(v)
            name = str(v)
            if (name.startswith("s") or name.startswith("t")) and len(name) > 1:
                try:
                    # bestimme Nummer des Start-/Zielknotens; ValueError, falls keine gültige Nummer vorhanden
                    farbe = int(name[1:])
                    if name.startswith("s"):
                        # weise Quellknoten die jeweilige Farbe zu
                        colorMap.append("deepskyblue")
                        continue
                    else:
                        # weise Zielknoten die jeweilige Farbe zu
                        colorMap.append(self.colors[(farbe-1) % len(self.colors)])
                        continue
                except ValueError:
                    # Knoten ist kein Start- oder Zielknoten, sondern ein anderer Knoten dessen Name mit "s" oder "t"
                    # beginnt
                    colorMap.append("paleturquoise")  # Standardfarbe für alle anderen Knoten
                    continue
            colorMap.append("paleturquoise")  # Standardfarbe für alle anderen Knoten

        #newpath = "GraphenGEXF"
        #if not os.path.exists(newpath):
        #    os.makedirs(newpath)
        #nx.write_gexf(graph, "GraphenGEXF\g1.gexf")
        #data = nx.readwrite.json_graph.node_link_data(graph)
        #with open('g1.json', 'w') as json_file:
        #    json.dump(data, json_file, indent=4, sort_keys=True)
        '''def read_json_file(filename):
            caps = []
            trav = []
            with open(filename) as f:
                json_string = json.load(f)
            for e in json_string["edgevalues"]:
                caps.append(e["capacity"])
                trav.append(e["traveltime"])
            graph = json_graph.node_link_graph(json_string)
            return graph, list(graph.nodes), list(graph.edges), caps, trav
        print(read_json_file('person.json'))'''

        for delta in self.items:
            for w in list(delta[1].keys()):
                self.E.append((delta[0], w))  # speichere alle Kanten in E
                self.r.append(list(self.items)[list(self.keys).index(delta[0])][1][w][1])  # speichere Reisezeiten in r
                self.nu.append(list(self.items)[list(self.keys).index(delta[0])][1][w][0])  # Kapazitäten in nu
                deq = deque()
                # "deques" enthält für jede Kante die entsprechende Warteschlange !!RECHTS IST VORNE, LINKS IST HINTEN!!
                self.deques.append(deq)
                graph.add_edge(delta[0], w)
        self.m = len(self.E)  # Anz. Kanten

        self.spieler = []  # Rechtecke für Spieler
        self.numbers = []  # Nummerierung der Spieler
        self.spielerV = []  # Rechtecke für virtuelle Spieler
        self.rec_height = 0.05
        self.rec_width = 0.04
        for i in self.I:
            # Rechtecke repräsentieren Spieler
            self.spieler.append(patches.Rectangle(self.posit['s{}'.format(self.ST[i][0])], self.rec_width,
                                                  self.rec_height, linewidth=1, edgecolor='black',
                                                  facecolor= self.colors[self.ST[i][1]-1 % len(self.colors)]))
            # hinzufügen Spieler zu subplot
            self.ax.add_patch(self.spieler[-1])
            # hinzufügen Nummern (für Spieler) zu subplot
            self.numbers.append(self.ax.text(self.posit['s{}'.format(self.ST[i][0])][0],
                                             self.posit['s{}'.format(self.ST[i][0])][1], str(i+1), fontsize=6))
        for e in self.E:
            self.spielerV.append(deque())

        nx.draw(graph, self.posit, node_color=colorMap, with_labels=True, alpha=0.8)  # zeichnen des Graphen
        self.fig.canvas.draw()

        self.graphReversed = self.reverse_graph(G)

        init = []
        for i in self.I:
            init.append(None)
        self.fm.append(init)  # Initialisierung "self.fm" für alle Spieler
        self.currentPos = self.num * [None]  # Initialisierung "self.currentPos" für alle Spieler

        self.c.append(np.zeros((self.num, self.m)))
        # gibt es eine virtuelle Queue beginnend bei theta = 0, so muss diese bei der Initialisierung von "self.c"
        # beachtet werden
        if 0 in self.start_queue:
            anfang = self.start_queue.index(0)
            while True:
                try:
                    # "index" durchläuft alle Indizes von in "start_queue" vorkommenden 0ern
                    index = self.start_queue[anfang:].index(0) + anfang
                except ValueError:  # alle 0er durchlaufen
                    break
                for i in self.I:
                    # setzen der Warteschlangenlängen durch virtuelle Spieler
                    self.c[0][i][self.E.index(self.kanten_queue[index])] += self.y0_queue[index]
                anfang += 1

        for i in self.I:  # Initialisierung "self.c" (Kosten)
            if variante == 'A':
                for e in self.E:
                    # zusätzliche Kosten durch virtuelle Spieler: (self.alpha[i]/(1 - self.alpha[i])) *
                    # ("Anzahl virtueller Spieler auf Kante e" / "Kantenkapazität von e")
                    self.c[0][i][self.E.index(e)] = (self.alpha[i]/(1 - self.alpha[i])) * \
                                                    (self.c[0][i][self.E.index(e)]/float(self.nu[self.E.index(e)])) + \
                                                    self.r[self.E.index(e)]
            elif variante == 'B':
                startknoten = list(set([f[0] for f in self.ST]))
                # zähle Spieler in den Quellknoten zum Zeitpunkt 0 für Variante 'B'. Diese Anzahlen werden in "j"
                # gespeichert und wirken sich gleichverteilt auf die Kosten der vom jeweiligen Startknoten ausgehenden
                # Kanten aus
                j = np.zeros(len(startknoten))
                for s in startknoten:
                    for i_prime in self.I:
                        # bestimme Anzahl Spieler mit Quellknoten "s" und Startzeitpunkt 0
                        if self.R[i_prime] == 0 and self.pos(i_prime, 0) == 's{}'.format(s):
                            j[startknoten.index(s)] += 1
                    # <Anzahl Spieler in "s"> durch <Anzahl von "s" ausgehende Kanten>
                    j[startknoten.index(s)] /= float(len(graph['s{}'.format(s)].keys()))
                for e in self.E:
                    name = str(e[0])
                    # Kosten der von Startknoten ausgehenden Kanten, diese berücksichtigen j-Werte der Spieler, die zum
                    # Zeitpunkt 0 im jeweiligen Quellknoten starten
                    if name.startswith("s") and int(name[1:]) in startknoten:
                        self.c[0][i][self.E.index(e)] = self.r[self.E.index(e)] + \
                                                        (self.alpha[i]/(1 - self.alpha[i])) * \
                                                        ((self.c[0][i][self.E.index(e)] +
                                                          j[startknoten.index(int(name[1:]))]) /
                                                         float(self.nu[self.E.index(e)]))
                    # alle anderen j-Werte sind 0 (da sich Spieler momentan nur in Quellknoten befinden)
                    else:
                        self.c[0][i][self.E.index(e)] = (self.alpha[i]/(1 - self.alpha[i])) * \
                                                        (self.c[0][i][self.E.index(e)]/float(self.nu[self.E.index(e)]))\
                                                        + self.r[self.E.index(e)]

        # erste Liste in "first" enthält alle vorkommenden Positionen, zweite Liste den Spieler mit kleinstem Index in
        # dieser Position
        first = [[], []]
        for i in self.I:  # Plot vor Start des Programms
            self.currentPos[i] = self.pos(i, 0)  # berechne Kosten Spieler "i" zum Zeitpunkt 0
            if self.currentPos[i] not in first[0]:  # erster Spieler in Position "self.currentPos[i]"
                first[0].append(self.currentPos[i])
                first[1].append(i)
                self.spieler[i].set_xy(self.posit[self.currentPos[i]]) # neue Position Spieler
                self.numbers[i].set_x(self.posit[self.currentPos[i]][0]+self.rec_width/3) # neue Position Nummer
                self.numbers[i].set_y(self.posit[self.currentPos[i]][1]+self.rec_height/4)
            else:  # alle weiteren Spieler in Position "self.currentPos[i]"
                count = self.currentPos[:i].count(self.currentPos[i])
                x, y = self.spieler[first[1][first[0].index(self.currentPos[i])]].get_xy()
                # setze Rechtecke der Spieler in gleicher Position übereinander
                self.spieler[i].set_xy((x,y+count*self.rec_height))
                self.numbers[i].set_x(x + self.rec_width/3)
                self.numbers[i].set_y(y+count*self.rec_height + self.rec_height/4)
        plt.title('Theta = 0')  # setze Titel des Plots
        self.fig.canvas.draw()

        AbfrageVirtuelleSpieler(self)

    def runner(self):
        """
        ist das Programm nicht pausiert, so wird jede Sekunde "run" aufgerufen, und "runner" erneut aufgerufen
        :return: 0, falls alle Spieler ihr Ziel erreicht haben,
                 1, falls Programm nicht pausiert. In diesem Fall wird "runner" erneut aufgerufen,
                 -1, falls Programm pausiert. In diesem Fall wird "runner" nicht erneut aufgerufen.
        """
        while self.button_win.get_zeit() < self.maxiter:
            if self.ankunft >= self.num:   # Abbruchkriterium
                self.button_win.nex.config(state="disabled")
                return 0
            if self.button_win.get_unpaused():
                if self.button_win.prev['state'] == 'disabled':
                    self.button_win.prev.config(state="normal")
                # Vergrößere Liste "self.leaveQueue", falls nötig (also falls "self.run(self.button_win.get_zeit() + 1)"
                # vorher noch nicht aufgerufen wurde)
                if len(self.leaveQueue) < self.button_win.get_zeit() + 2:
                    init_p = []
                    init_m = []
                    for i in self.I:
                        init_p.append(None)
                        init_m.append(None)
                    self.fm.append(init_m)
                    self.fp.append(init_p)
                    self.leaveQueue.append([])
                self.button_win.set_zeit(self.button_win.get_zeit() + 1)
                self.run(self.button_win.get_zeit())
                self.button_win.after(1000, self.runner)
                return 1
            else:
                self.button_win.after(1000, self.runner)
                #print("pausiert")
                return -1

    def g(self, i, theta, e):
        """
        bestimmt, ob sich Spieler "i" zum Zeitpunkt "theta" auf Kante "e" befindet
        :param i: Spielerindex
        :param theta: Zeitpunkt
        :param e: Kante
        :return: 1, falls sich Spieler "i" zum Zeitpunkt "theta" auf Kante "e" befindet
                 0, sonst.
        """
        summe = 0
        for t in range(theta + 1):
            if self.fp[t][i] == e:
                summe += 1
            elif self.fm[t][i] == e:
                summe -= 1
        return summe

    def pos(self, i, theta):
        """
        bestimmt Position von Spieler "i" zum Zeitpunkt "theta".
        :param i: Spielerindex
        :param theta: Zeitpunkt
        :return: Position von Spieler "i" zum Zeitpunkt "theta"
        """
        # Fall 4 der Definition von Pos_i(theta) aus dem paper
        if theta <= self.R[i]:
            return 's{}'.format(self.ST[i][0])
        # Fall 5
        if self.z[i] >= 0 and self.z[i] < theta:
            return 't{}'.format(self.ST[i][1])
        # hat diese Funktion hier noch nicht terminiert, so sind wir in einem der Fälle mit "g(i,theta - 1, e)" = 1
        # für genau eine Kante e (siehe Beweis von Lemma 1.4)
        # bestimme also nun diese Kante e
        for edge in self.E:
            if self.g(i, theta - 1, edge):
                e = edge
                break
        # Fall 1
        if self.fm[theta][i] == e:
            return e[1]
        insert = theta - 1
        while True:
            if self.fp[insert][i] == e:
                break
            else:
                insert -= 1
        for t in reversed(range(insert, theta)):
            if self.spieler[i] in self.leaveQueue[t][self.E.index(e)]:
                return (e, theta - t)  # Fall 2
        return (e,0)  # Fall 3

    def reverse_graph(self, graph):
        """
        ersetzt alle Kanten in "graph" (gegeben als Dictionary) durch die jeweilige Rückwärtskante, Kosten bleiben
        gleich
        :param graph: Eingabegraph, dessen Kantenorientierungen vertauscht werden sollen
        :return: "graphReversed", entspricht "graph" mit vertauschten Kantenorientierungen als Dictionary
        """
        graphReversed = {}
        for v in graph.keys():
            dictionary = {}
            for key, delta in graph.items():
                for node in delta.keys():
                    if node == v:
                        dictionary[key] = delta[node]
            graphReversed[v] = dictionary
        return graphReversed

    def varianten(self, alpha, G, var='A', anz=None):
        """
        # berechnet Kosten abhängig von der übergebenen Variante
        :param alpha: aus Intervall [0,1]. Legt Gewichtung des Einflusses der Reise- und Wartezeiten auf die Kosten
         fest, beispielsweise:
            alpha = 0: nur die Reisedauer bestimmt Kosten
            alpha = 1: nur die Wartezeit bestimmt Kosten
            alpha = 1/2: Reisedauer und Wartezeit nehmen gleichen Einfluss auf Kosten
        :param G: Graph als Dictionary
        :param var: Variante zur Berechnung der Kantenkosten, Standardwert: 'A'
        :param anz: nur für Variante 'B'. Liste über alle Knoten (selbe Indizierung), enthält Anzahl Spieler, die sich
         momentan in den jeweiligen Knoten befinden.
        :return: "kosten": Liste über alle Spieler, enthält für jeden Spieler eine Liste mit allen Kantenkosten (gleiche
         Indizierung wie "self.E")
        """
        wartezeit = np.zeros(self.m)
        kosten = np.zeros((self.num, self.m))

        for e in self.E:
            if var == 'A':
                wartezeit[self.E.index(e)] = len(self.deques[self.E.index(e)])/float(self.nu[self.E.index(e)])
            elif var == 'B':
                wartezeit[self.E.index(e)] = (len(self.deques[self.E.index(e)]) +
                                              anz[self.V.index(e[0])]/ float(len(G[e[0]].keys()))) /\
                                             float(self.nu[self.E.index(e)])
            elif var == 'C':
                raise TypeError('C fehlt noch')
            # wie Variante 'A', mit dem Unterschied, dass Wartezeiten aufgerundet werden, da diskretes Modell
            elif var == 'D':
                wartezeit[self.E.index(e)] = math.ceil(len(self.deques[self.E.index(e)]) /
                                                       float(self.nu[self.E.index(e)]))
            else:
                raise TypeError('Ungültige Variante')

        for i in self.I:
            for e in self.E:
                # berechne Kosten aus Wartezeiten
                kosten[i][self.E.index(e)] = 1/(1 - alpha[i]) * ((1 - alpha[i]) * self.r[self.E.index(e)] +
                                                                 alpha[i]*wartezeit[self.E.index(e)])
        return kosten

    def remove_artist(self, ind):
        """
        entfernt virtuellen Spieler aus Plot
        :param ind: Index des zu entfernenden virtuellen Spielers
        :return: kein Rückgabewert
        """
        artist = self.spielerV[ind].popleft()
        artist.remove()
        return

    # Quelle: http://www.gilles-bertrand.com/2014/03/dijkstra-algorithm-python-example-source-code-shortest-path.html
    def dijkstra(self, graph, src, dest, i, theta, visited=[], distances={}, predecessors={}):
        """
        Berechnet rekursiv kürzesten Weg und dessen Kosten von "src" zu "dest" im Graph "graph" für die Kosten vom
        Spieler mit Index "i" zum Zeitpunkt "theta"
        :param graph: Graph als Dictionary
        :param src: Startknoten
        :param dest: Zielknoten
        :param i: Spielerindex
        :param theta: Zeitpunkt
        :param visited: Liste der bereits besuchten Knoten, anfangs leer, muss nicht beachtet werden, da nur intern für
         die Funktion benötigt
        :param distances: Liste der bereits berechneten Distanzen aller Knoten zu "src", anfangs leer, muss nicht
         beachtet werden, da nur intern für die Funktion benötigt
        :param predecessors: Liste der bereits ermittelten Vorfahren aller Knoten, anfangs leer, muss nicht beachtet
         werden, da nur intern für die Funktion benötigt
        :return: "output": Liste über alle Knoten mit deren Distanzen zur ersten(!) "src", also zu dem Knoten, der beim
                           ersten Aufruf von "dijkstra" als "src" übergeben wurde
                 "self.dijkstra(graph, x, dest, i, theta, visited, distances, predecessors)": sorgt für rekursives
                           Aufrufen dieser Funktion, endet mit return von "output"
        """
        """ calculates a shortest path tree routed in src
        """
        # a few sanity checks
        if src not in graph:
            raise TypeError('The root of the shortest path tree cannot be found')
        if dest not in graph:
            raise TypeError('The target of the shortest path cannot be found')
            # ending condition
        if src == dest:
            output = []
            if len(distances) == 0:
                for v in self.V:
                    if v == dest:
                        output.append(0)
                    else:
                        output.append(float('inf'))
                return output
            elif dest not in distances:
                raise TypeError('The target is not reachable from the root')

            for v in self.V:
                if v in distances:
                    output.append(distances[v])
                else:
                    output.append(float('inf'))
            return output
        else:
            # if it is the initial run, initializes the cost
            if not visited:
                distances[src] = 0
            # visit the neighbors
            for neighbor in graph[src]:
                if neighbor not in visited:
                    new_distance = distances[src] + self.c[theta][i][self.E.index((neighbor,src))]
                    if new_distance < distances.get(neighbor,float('inf')):
                        distances[neighbor] = new_distance
                        predecessors[neighbor] = src
            # mark as visited
            if src not in visited:
                visited.append(src)
            # now that all neighbors have been visited: recurse
            # select the non visited node with lowest distance 'x'
            # run Dijkstra with src='x'
            unvisited = {}
            for k in graph:
                if k not in visited:
                    unvisited[k] = distances.get(k,float('inf'))
            x=min(unvisited, key=unvisited.get)
            return self.dijkstra(graph, x, dest, i, theta, visited, distances, predecessors)

    def add_virtual_queue(self, theta, kanten_queue, start_queue, ende_queue, y0_queue, kap):
        """
        fügt virtuelle Spieler zu den jeweiligen Warteschlangen hinzu. Diese werden über die Liste "self.deques" ver-
        waltet. Diese Liste ist nach der Indizierung von "self.E", also entsprechend der Kantenindizes, indiziert.
        :param theta: aktueller Zeitpunkt
        :param kanten_queue: Liste aller Kanten, die irgendwann einen virtuellen Einfluss haben
        :param start_queue: Liste aller Startzeitpunkte der virtuellen Einflüsse
        :param ende_queue: Liste aller Endzeitpunkte der virtuellen Einflüsse
        :param y0_queue: Liste aller (ganzzahligen) Einflusswerte der virtuellen Einflüsse
        :param kap: Liste der Restkapazitäten aller Kanten
        :return: kein Rückgabewert
        """
        for t in range(len(start_queue)):
            if start_queue[t] <= theta and theta <= ende_queue[t]:
                for y in range(1, y0_queue[t] + 1):
                    # füge "theta" als virtuellen 'Spieler' in Queue ein
                    self.deques[self.E.index(kanten_queue[t])].appendleft(theta)
                    # speichere virtuelle Spieler in dict., um diese aus Plot entfernen zu können
                    self.spielerV[self.E.index(kanten_queue[t])].appendleft(
                        patches.Rectangle((0.91*self.posit[kanten_queue[t][0]][0] +
                                           0.09*self.posit[kanten_queue[t][1]][0],0.91 *
                                           self.posit[kanten_queue[t][0]][1] + 0.09*self.posit[kanten_queue[t][1]][1] +
                                           self.rec_height*
                                           (len(self.deques[self.E.index(kanten_queue[t])])
                                            -kap[self.E.index(kanten_queue[t])] -1)), self.rec_width, self.rec_height,
                                          linewidth=1, edgecolor='black', facecolor= 'white', alpha=0.5))
                    # füge virtuellen Spieler zu Plot hinzu
                    self.ax.add_patch(self.spielerV[self.E.index(kanten_queue[t])][0])
        return

    def run(self, theta):
        """
        Führt einen Vorwärtsschritt aus, zeigt also die Positionen aller Spieler zum Zeitpunkt "theta". Dazu wird
        erst bestimmt, ob die Funktion "run" bereits mit genau diesem "theta" zuvor schon aufgerufen wurde
        ("rerun"=True), oder nicht ("rerun"=False"). Im Falle eines reruns müssen nämlich die Distanzen und die Werte
        für f^+, f^- nicht erneut berechnet werden, da diese bereits im ersten Aufruf von "run(theta)" bestimmt wurden.
        In beiden Fällen werden die Warteschlangen an den Zeitpunkt "theta" angepasst, sowie der Plot aktualisiert.
        :param theta: aktueller Zeitpunkt
        :return: kein Rückgabewert
        """
        newDists = []
        # erste Liste in "first" enthält alle vorkommenden Positionen, zweite Liste den Spieler mit kleinstem Index in
        # dieser Position
        first = [[], []]
        residualCap = self.nu.copy()
        if len(self.label) < theta or theta == 0:
            rerun = False
        else:
            rerun = True

        if theta > 0:
            # hinzufügen virtueller Warteschlangen (festgelegt durch: "kanten_queue", "start_queue", "ende_queue",
            # "y0_queue")
            self.add_virtual_queue(theta -1, self.kanten_queue, self.start_queue, self.ende_queue, self.y0_queue,
                                 residualCap)

        if not rerun:
            if theta == 0:
                for i in self.I:
                    self.currentPos[i] = self.pos(i,0)
            for i in self.I:
                if self.z[i] == -1:
                    if self.currentPos[i] in self.V:
                        if theta > 0:
                            newDists.append(self.dijkstra(self.graphReversed, "t{}".format(self.ST[i][1]),
                                                          self.currentPos[i],i,theta -1,visited = [], distances={}))
                        else:
                            newDists.append(self.dijkstra(self.graphReversed, "t{}".format(self.ST[i][1]),
                                                          self.currentPos[i],i,theta,visited = [], distances={}))
                    elif self.currentPos[i][1] == 0:
                        # befindet sich Spieler "i" momentan auf einer Kante, so werden seine Labels nicht benötigt und
                        # daher auch nicht berechnet, um Rechenzeit zu sparen
                        newDists.append(self.n*[None])
                    else:
                        # befindet sich Spieler "i" momentan auf einer Kante, so werden seine Labels nicht benötigt und
                        # daher auch nicht berechnet, um Rechenzeit zu sparen
                        newDists.append(self.n*[None])
                else:
                    # ist Spieler "i" bereits bei seiner Senke angekommen, so werden die Labels nicht benötigt
                    newDists.append(self.n*[None])


        if not rerun and theta > 0:
            # füge weitere Zeile in "self.label" ein, diese enthält die self.labels aller Spieler zum Zeitpunkt "theta"
            self.label.append(newDists)

            J = self.n*[0]  # Hilfsvariable für Variante "B"

            # verwendet den Wert von "self.currentPos[i]" vom VORHERIGEN Zeitpunkt um "fp" für diesen Zeitpunkt
            # bestimmen zu können
            for i in self.I:
                if theta -1 < self.R[i]:
                    continue
                # überprüfe, ob Spieler "i" sich in einem Knoten befindet, den er zur nächsten Zeiteinheit verlassen
                # muss
                if self.currentPos[i] in self.V:
                    if self.variante == 'B':
                        # Anzahl Spieler, die sich gerade in einem Knoten befinden wird gespeichert (auch in welchem
                        # Knoten sie sich befinden) für Variante B
                        J[self.V.index(self.currentPos[i])] += 1
                    if self.currentPos[i] == "t{}".format(self.ST[i][1]):
                        continue
                    currentMin = float('inf')
                    currentNode = None
                    # bestimme Nachbarknoten mit niedrigstem Label
                    for v in self.G[self.currentPos[i]].keys():
                        if self.label[theta -1][i][self.V.index(v)] + \
                                self.c[theta -1][i][self.E.index((self.currentPos[i], v))] < currentMin:
                            currentMin = self.label[theta -1][i][self.V.index(v)] + \
                                         self.c[theta -1][i][self.E.index((self.currentPos[i], v))]
                            currentNode = v
                    # betritt Kante mit niedrigstem Label
                    self.fp[theta -1][i] = (self.currentPos[i], currentNode)
                    if None in self.fp[theta -1][i]:
                        raise TypeError('Knoten {} von {} nicht erreichbar!'.format("t{}".format(self.ST[i][1]),
                                                                                    self.currentPos[i]))

        for e in self.E:  # bestimme aktuelle Queue für jede Kante
            if not rerun:
                self.leaveQueue[theta -1].append([])
            residualCap[self.E.index(e)] -= min(len(self.deques[self.E.index(e)]), self.nu[self.E.index(e)])
            for out in range(min(len(self.deques[self.E.index(e)]), self.nu[self.E.index(e)])):
                nextPlayer = self.deques[self.E.index(e)].pop()
                if not rerun:
                    self.leaveQueue[theta -1][self.E.index(e)].append(nextPlayer)
                try:
                    int(nextPlayer)  # TypeError, falls "nextPlayer" nicht-virtueller Spieler
                    if len(self.spielerV[self.E.index(e)]) > 0:
                        v = self.spielerV[self.E.index(e)].pop()  # entferne virtuellen Spieler
                        v.remove()  # entferne virtuellen Spieler aus plot
                except TypeError:
                    # umfärben von rot auf ursprüngliche Farbe
                    self.spieler[self.spieler.index(nextPlayer)].set_facecolor(
                        self.colors[self.ST[self.spieler.index(nextPlayer)][1]-1 % len(self.colors)])
            if theta > 0:
                for i in self.I:
                    if self.fp[theta-1][i] == e:
                        if residualCap[self.E.index(e)] == 0:
                            self.deques[self.E.index(e)].appendleft(self.spieler[i])  # hinzufügen Spieler "i" zu Queue
                        else:
                            if not rerun:
                                self.leaveQueue[theta-1][self.E.index(e)].append(self.spieler[i])
                            residualCap[self.E.index(e)] -= 1

        if not rerun and theta > 0:
            for i in self.I:
                # verwendet den Wert von "self.currentPos[i]" vom VORHERIGEN Zeitpunkt um "fm" für den aktuellen
                # Zeitpunkt bestimmen zu können
                if self.currentPos[i] in self.V:
                    if self.fp[theta-1][i] is not None:
                        if self.r[self.E.index(self.fp[theta-1][i])] == 1:
                            for spielerliste in self.leaveQueue[theta-1][self.E.index(self.fp[theta-1][i])]:
                                try:
                                    int(spielerliste)  # überspringe virtuelle Spieler
                                except TypeError:  # "spielerliste" ist echter Spieler
                                    if i == self.spieler.index(spielerliste):
                                        # falls die Kante "fp[theta-1][i]", welche Spieler "i" zum Zeitpunkt "theta-1"
                                        # betritt, die Reisedauer 1 hat und Länge der Warteschlange kleiner als die
                                        # Kantenkapazität ist, so wird Spieler "i" die Kante zum Zeitpunkt "theta"
                                        # wieder verlassen
                                        self.fm[theta][i] = self.fp[theta-1][i]
                else:
                    # überprüft, ob sich Spieler "i" zum vorherigen Zeitpunkt auf einer Kante, genau eine Zeiteinheit
                    # vor dem Ende der Kante, befunden hat
                    if self.currentPos[i][1] == self.r[self.E.index(self.currentPos[i][0])] - 1:
                        if self.r[self.E.index(self.currentPos[i][0])] != 1:
                            # falls ja, wird Spieler "i" Kante "currentPos[0]" nun verlassen
                            self.fm[theta][i] = self.currentPos[i][0]
                        else:
                            # falls die Kante die Reisedauer 1 besitzt, so könnte sich der Austritt aus der Kante durch
                            # die queue noch verzögern, dies wird nun überprüft
                            for spielerliste in self.leaveQueue[theta-1][self.E.index(self.currentPos[i][0])]:
                                try:
                                    int(spielerliste)  # überspringe virtuelle Spieler
                                except TypeError:
                                    if i == self.spieler.index(spielerliste):
                                        self.fm[theta][i] = self.currentPos[i][0]

        for i in self.I:
            if self.z[i] == -1 or self.z[i] >= theta:
                self.currentPos[i] = self.pos(i, theta)
                # print("Position Spieler " + str(i+1) + " zum Zeitpunkt " + str(theta) + " : " +
                # str(self.currentPos[i]))
                if self.currentPos[i] in self.V:
                    # Spieler werden nicht mehr angezeigt, wenn Senke erreicht
                    if self.currentPos[i] == "t{}".format(self.ST[i][1]):
                        self.spieler[i].set_visible(False)
                        self.numbers[i].set_visible(False)
                        if self.z[i] == -1:
                            self.z[i] = theta  # setze Ankunftszeit
                        if self.z[i] == theta:  # aktualisiere Abbruchkriterium
                            self.ankunft += 1
                    elif self.currentPos[i] not in first[0]:
                        first[0].append(self.currentPos[i])
                        first[1].append(i)
                        # neue Position Spieler
                        self.spieler[i].set_xy(self.posit[self.currentPos[i]])
                        # neue Position Nummer
                        self.numbers[i].set_x(self.posit[self.currentPos[i]][0]+self.rec_width/3)
                        self.numbers[i].set_y(self.posit[self.currentPos[i]][1]+self.rec_height/4)
                    else:
                        # zähle Spieler in gleicher Position mit kleinerem Index
                        count = self.currentPos[:i].count(self.currentPos[i])
                        # ziehe Spieler die an ihrem Zielknoten angekommen (und nicht sichtbar) sind, ab
                        for sp in self.I[:i]:
                            if self.z[sp] != -1 and self.z[sp] <= theta and self.currentPos[sp] == self.currentPos[i]:
                                count -= 1
                        x, y = self.spieler[first[1][first[0].index(self.currentPos[i])]].get_xy()
                        # setze Rechtecke der Spieler in gleicher Position übereinander
                        self.spieler[i].set_xy((x, y+count*self.rec_height))
                        self.numbers[i].set_x(x + self.rec_width/3)
                        self.numbers[i].set_y(y+count*self.rec_height + self.rec_height/4)
                elif self.currentPos[i][1] != 0:
                    if self.currentPos[i] not in first[0]:
                        first[0].append(self.currentPos[i])
                        first[1].append(i)
                        self.spieler[i].set_xy((1-self.currentPos[i][1]/self.r[self.E.index(self.currentPos[i][0])]) *
                                               self.posit[self.currentPos[i][0][0]] +
                                               (self.currentPos[i][1]/self.r[self.E.index(self.currentPos[i][0])]) *
                                               self.posit[self.currentPos[i][0][1]])  # neue Position Spieler auf Kante
                        self.numbers[i].set_x(((1-self.currentPos[i][1]/self.r[self.E.index(self.currentPos[i][0])]) *
                                               self.posit[self.currentPos[i][0][0]] +
                                               (self.currentPos[i][1]/self.r[self.E.index(self.currentPos[i][0])]) *
                                               self.posit[self.currentPos[i][0][1]])[0]+self.rec_width/3)
                        self.numbers[i].set_y(((1-self.currentPos[i][1]/self.r[self.E.index(self.currentPos[i][0])]) *
                                               self.posit[self.currentPos[i][0][0]] +
                                               (self.currentPos[i][1]/self.r[self.E.index(self.currentPos[i][0])]) *
                                               self.posit[self.currentPos[i][0][1]])[1]+self.rec_height/4)
                    else:
                        count = self.currentPos[:i].count(self.currentPos[i])
                        # Koordinaten des ersten Spielers
                        x, y = self.spieler[first[1][first[0].index(self.currentPos[i])]].get_xy()
                        # setze Rechtecke der Spieler in gleicher Position übereinander
                        self.spieler[i].set_xy((x,y+count*self.rec_height))
                        self.numbers[i].set_x(x+self.rec_width/3)
                        self.numbers[i].set_y(y+count*self.rec_height+self.rec_height/4)

        for e in self.E:  # setzen der Positionen von Spielern in Warteschlange
            j = len(self.deques[self.E.index(e)]) -1  # erster Index in der deque
            vcount = 0
            while j >= 0:
                try:
                    # TypeError, falls Spieler "self.deques[self.E.index(e)][j]" KEIN virtueller Spieler
                    int(self.deques[self.E.index(e)][j])
                    # setze Position virtueller Spieler
                    self.spielerV[self.E.index(e)][vcount].set_xy((0.91*self.posit[e[0]][0] + 0.09*self.posit[e[1]][0],
                                                                   0.91*self.posit[e[0]][1] + 0.09*self.posit[e[1]][1] +
                                                                   self.rec_height *
                                                                   (len(self.deques[self.E.index(e)]) -1 -j)))
                    vcount += 1
                except TypeError:
                    self.deques[self.E.index(e)][j].set_facecolor('red')
                    # setze Position von Spieler "self.deques[self.E.index(e)][j]"
                    self.spieler[self.spieler.index(self.deques[self.E.index(e)][j])].set_xy((0.91*self.posit[e[0]][0] +
                                                                                              0.09*self.posit[e[1]][0],
                                                                                              0.91*self.posit[e[0]][1] +
                                                                                              0.09*self.posit[e[1]][1] +
                                                    self.rec_height*(len(self.deques[self.E.index(e)]) -1 -j)))
                    self.numbers[self.spieler.index(self.deques[self.E.index(e)][j])].set_x(0.91*self.posit[e[0]][0] +
                                                                                            0.09*self.posit[e[1]][0] +
                                                                                            self.rec_width/3)
                    self.numbers[self.spieler.index(self.deques[self.E.index(e)][j])].set_y(0.91*self.posit[e[0]][1] +
                                                                                            0.09*self.posit[e[1]][1] +
                                                                                            self.rec_height/4 +
                                                                                            self.rec_height *
                                                                            (len(self.deques[self.E.index(e)]) -1 -j))
                j -= 1

        plt.title('Theta = {}'.format(theta))  # setze neuen Titel des Plots
        self.fig.canvas.draw()  # redraw

        if not rerun and theta > 0:
            self.c.append(self.varianten(self.alpha, self.G, var=self.variante, anz=J))
        return

    def runback(self, theta):
        """
        wie "run"-Methode zum Zeitpunkt "theta", angepasst für den Fall, dass "Zurück" - Button gedrückt wurde. Im
        Unterschied zu einem Aufruf von "run(theta)" mit "rerun"=True, werden hier die Warteschlangen nicht verändert,
        da dies bereits beim Klicken des Buttons "Zurück" in der Funktion "ButtonWindowFrame.back()" geschieht.
        :param theta: aktueller Zeitpunkt
        :return: kein Rückgabewert
        """
        # erste Liste enthält alle vorkommenden Positionen, zweite Liste den Spieler mit kleinstem Index in jeder dieser
        # Positionen
        first = [[], []]

        for i in self.I:
            if self.z[i] == -1 or self.z[i] > theta:
                self.currentPos[i] = self.pos(i, theta)
                # print("Position Spieler " + str(i+1) + " zum Zeitpunkt " + str(theta) + " : " + str(currentPos[i]))
                if self.currentPos[i] in self.V:
                    if self.currentPos[i] != "t{}".format(self.ST[i][1]) and self.currentPos[i] not in first[0]:
                        first[0].append(self.currentPos[i])
                        first[1].append(i)
                        # neue Position Spieler
                        self.spieler[i].set_xy(self.posit[self.currentPos[i]])
                        # neue Position Nummer
                        self.numbers[i].set_x(self.posit[self.currentPos[i]][0]+self.rec_width/3)
                        self.numbers[i].set_y(self.posit[self.currentPos[i]][1]+self.rec_height/4)
                    else:
                        # zähle Spieler in gleicher Position mit kleinerem Index
                        count = self.currentPos[:i].count(self.currentPos[i])
                        # ziehe Spieler die an ihrem Zielknoten angekommen (und nicht sichtbar) sind, ab
                        for sp in self.I[:i]:
                            if self.z[sp] != -1 and self.z[sp] <= theta and self.currentPos[sp] == self.currentPos[i]:
                                count -= 1
                        x,y = self.spieler[first[1][first[0].index(self.currentPos[i])]].get_xy()
                        # setze Rechtecke der Spieler in gleicher Position übereinander
                        self.spieler[i].set_xy((x, y+count*self.rec_height))
                        self.numbers[i].set_x(x + self.rec_width/3)
                        self.numbers[i].set_y(y+count*self.rec_height + self.rec_height/4)
                elif self.currentPos[i][1] != 0:
                    if self.currentPos[i] not in first[0]:
                        first[0].append(self.currentPos[i])
                        first[1].append(i)
                        # neue Position Spieler auf Kante
                        self.spieler[i].set_xy((1-self.currentPos[i][1]/self.r[self.E.index(self.currentPos[i][0])]) *
                                               self.posit[self.currentPos[i][0][0]] +
                                               (self.currentPos[i][1]/self.r[self.E.index(self.currentPos[i][0])]) *
                                               self.posit[self.currentPos[i][0][1]])
                        self.numbers[i].set_x(((1-self.currentPos[i][1]/self.r[self.E.index(self.currentPos[i][0])]) *
                                               self.posit[self.currentPos[i][0][0]] +
                                               (self.currentPos[i][1]/self.r[self.E.index(self.currentPos[i][0])]) *
                                               self.posit[self.currentPos[i][0][1]])[0]+self.rec_width/3)
                        self.numbers[i].set_y(((1-self.currentPos[i][1]/self.r[self.E.index(self.currentPos[i][0])]) *
                                               self.posit[self.currentPos[i][0][0]] +
                                               (self.currentPos[i][1]/self.r[self.E.index(self.currentPos[i][0])]) *
                                               self.posit[self.currentPos[i][0][1]])[1]+self.rec_height/4)
                    else:
                        count = self.currentPos[:i].count(self.currentPos[i])
                        x, y = self.spieler[first[1][first[0].index(self.currentPos[i])]].get_xy()
                        # setze Rechtecke der Spieler in gleicherPosition übereinander
                        self.spieler[i].set_xy((x,y+count*self.rec_height))
                        self.numbers[i].set_x(x+self.rec_width/3)
                        self.numbers[i].set_y(y+count*self.rec_height+self.rec_height/4)

        for e in self.E:  # setzen der Positionen von Spielern in Warteschlange
            j = len(self.deques[self.E.index(e)]) -1  # erster Index in der deque
            vcount = 0
            while j >= 0:
                try:
                    # TypeError, falls Spieler "self.deques[self.E.index(e)][j]" KEIN virtueller Spieler
                    int(self.deques[self.E.index(e)][j])
                    # setze Position virtueller Spieler
                    self.spielerV[self.E.index(e)][vcount].set_xy((0.91*self.posit[e[0]][0] + 0.09*self.posit[e[1]][0],
                                                                   0.91*self.posit[e[0]][1] + 0.09*self.posit[e[1]][1] +
                                                    self.rec_height*(len(self.deques[self.E.index(e)]) -1 -j)))
                    vcount += 1
                except TypeError:
                    self.deques[self.E.index(e)][j].set_facecolor('red')
                    #setze Position von Spieler "self.deques[self.E.index(e)][j]"
                    self.spieler[self.spieler.index(self.deques[self.E.index(e)][j])].set_xy((0.91*self.posit[e[0]][0] +
                                                                                              0.09*self.posit[e[1]][0],
                                                                                              0.91*self.posit[e[0]][1] +
                                                                                              0.09*self.posit[e[1]][1] +
                                                                                              self.rec_height *
                                                                    (len(self.deques[self.E.index(e)]) -1 -j)))
                    self.numbers[self.spieler.index(self.deques[self.E.index(e)][j])].set_x(0.91*self.posit[e[0]][0] +
                                                                                            0.09*self.posit[e[1]][0] +
                                                                                            self.rec_width/3)
                    self.numbers[self.spieler.index(self.deques[self.E.index(e)][j])].set_y(0.91*self.posit[e[0]][1] +
                                                                                            0.09*self.posit[e[1]][1] +
                                                                                            self.rec_height/4 +
                                                                                            self.rec_height *
                                                                              (len(self.deques[self.E.index(e)]) -1 -j))
                j -= 1

        plt.title('Theta = {}'.format(theta))  # setze neuen Titel des Plots
        self.fig.canvas.draw()  # redraw
        return
