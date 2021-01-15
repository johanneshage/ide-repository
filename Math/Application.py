from Graphics.ButtonWindowFrame import ButtonWindowFrame
import numpy as np
from collections import deque


class Application:
    """
    Objekt der Klasse "Application" stellt eine Anwendung des Hauptprogramms dar
    """

    def __init__(self, G, R, ST, alpha, posit, variante='A', kanten_queue=[], start_queue=[],
                 ende_queue=[], y0_queue=[]):
        """
        Initialisiert Graph und alle Variablen, sowie Kantenkosten und erzeugt Plot zum Zeitpunkt 0 (vor eventuellem
        Hinzufügen weiterer virtueller Spieler durch "AbfrageVirtuelleSpieler.add_parameter()" ). Alle Eingabeparameter
        außer "posit" und "variante" werden in "data.py" spezifiziert, oder gegebenenfalls in "Main.py", und haben
        folgende Funktion:
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
         eingelesen, falls in "data.py" kein Graph angegeben
        :param variante: gibt Variante zur Kostenberechnung an; Möglichkeiten: 'A', 'B', 'C', 'D'. Standardwert: 'A'.
        :param kanten_queue: Liste, die alle Kanten mit virtueller Warteschlange, als Tupel der Form ('v','w'), enthält
        :param start_queue: Liste die zu den Einträgen in "kanten_queue" die entsprechenden Startzeitpunkte des
         virtuellen Einflusses enthält (i-ter Eintrag in "start_queue" bezieht sich auf i-ten Eintrag in "kanten_queue"
        :param ende_queue: Analog zu "start_queue" ist dies eine Liste, die die Endzeitpunkte des virtuellen Einflusses
         enthält
        :param y0_queue: Liste mit den Einflussgrößen des virtuellen Flusses, indiziert wie "kanten_queue",
         "start_queue", "ende_queue"
        """
        self.kanten_queue = kanten_queue  # Liste der Kanten mit virtuellen Spielern
        self.start_queue = start_queue  # Liste der Startzeitpunkte der virtuellen Zuflüsse
        self.ende_queue = ende_queue  # Liste der Endzeitpunkte der virtuellen Zuflüsse
        self.y0_queue = y0_queue  # Liste der virtuellen Zuflussmengen
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
        self.fm.append(self.num * [None])
        self.fp.append(self.num * [None])
        self.ST = ST  # Zuordnung Start- und Zielknoten zu Spieler
        self.alpha = alpha
        self.R = R  # Startzeitpunkte
        self.variante = variante
        self.G = G  # Graph
        self.ankunft = 0  # zählt angekommene Spieler
        self.zeitpunkt = 0
        self.unpaused = True  # Programm zu Beginn nicht pausiert
        self.deques = []  # Liste aller deques (= Warteschlangen)
        # Liste "leaveQueue" enthält für jeden Zeitpunkt für jede Kante eine Liste der Spieler, die zum gegebenen
        # Zeitpunkt die Warteschlange dieser Kante verlassen
        self.leaveQueue = [[]]
        self.items = G.items()
        self.keys = G.keys()

        for delta in self.items:
            for w in list(delta[1].keys()):
                self.E.append((delta[0], w))  # speichere alle Kanten in E
                self.r.append(list(self.items)[list(self.keys).index(delta[0])][1][w][1])  # speichere Reisezeiten in r
                self.nu.append(list(self.items)[list(self.keys).index(delta[0])][1][w][0])  # Kapazitäten in nu
                deq = deque()
                # "deques" enthält für jede Kante die entsprechende Warteschlange !!RECHTS IST VORNE, LINKS IST HINTEN!!
                self.deques.append(deq)
        self.m = len(self.E)  # Anz. Kanten

        self.graphReversed = self.reverse_graph(G)

        self.fm.append(self.num*[None])  # Initialisierung "self.fm" für alle Spieler
        self.currentPos = []
        for i in self.I:  # Initialisierung "self.currentPos" -> Spieler befinden sich in Quellen
            self.currentPos.append('s{}'.format(self.ST[i][0]))

        self.c.append(np.zeros((self.num, self.m)))
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
                    j[startknoten.index(s)] /= float(len(self.button_win.graph['s{}'.format(s)].keys()))
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

        self.button_win = ButtonWindowFrame(self.E, self.V, self.ST, self.r, posit)  # erzeuge Buttonleiste
        # weise Buttons Funktionen zu
        self.button_win.prev.configure(command=self.back)
        self.button_win.pause.configure(command=self.pause)
        self.button_win.nex.configure(command=self.weiter)

    def runner(self):
        """
        ist das Programm nicht pausiert, so wird jede Sekunde "run" aufgerufen, und "runner" erneut aufgerufen
        :return: 0, falls alle Spieler ihr Ziel erreicht haben,
                 1, falls Programm nicht pausiert. In diesem Fall wird "runner" erneut aufgerufen,
                 -1, falls Programm pausiert. In diesem Fall wird "runner" nicht erneut aufgerufen.
        """
        while self.zeitpunkt < self.maxiter:
            if self.ankunft >= self.num:   # Abbruchkriterium
                self.button_win.nex.config(state="disabled")
                return 0
            if self.unpaused:
                if self.button_win.prev['state'] == 'disabled':
                    self.button_win.prev.config(state="normal")
                # Vergrößere Liste "self.leaveQueue", falls nötig (also falls "self.run(self.button_win.get_zeit() + 1)"
                # vorher noch nicht aufgerufen wurde)
                if len(self.leaveQueue) < self.zeitpunkt + 2:
                    self.fm.append(self.num * [None])
                    self.fp.append(self.num * [None])
                    self.leaveQueue.append([])
                self.zeitpunkt += 1
                self.run(self.zeitpunkt)
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
            if (i, 'r') in self.leaveQueue[t][self.E.index(e)]:
                return e, theta - t  # Fall 2
        return e, 0  # Fall 3

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
                wartezeit[self.E.index(e)] = np.ceil(len(self.deques[self.E.index(e)]) /float(self.nu[self.E.index(e)]))
            else:
                raise TypeError('Ungültige Variante')

        for i in self.I:
            for e in self.E:
                # berechne Kosten aus Wartezeiten
                kosten[i][self.E.index(e)] = 1/(1 - alpha[i]) * ((1 - alpha[i]) * self.r[self.E.index(e)] +
                                                                 alpha[i]*wartezeit[self.E.index(e)])
        return kosten

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
                    self.button_win.add_virtual_player(kanten_queue[t], len(self.deques[self.E.index(kanten_queue[t])])
                                                                         -kap[self.E.index(kanten_queue[t])] -1)
        return

    def pause(self):
        """
        Bei Aufruf wird Pause-Status geändert (Fortsetzen falls pausiert und andersrum)
        :return: Kein Rückgabewert
        """
        if self.unpaused:
            self.unpaused = False
            self.button_win.set_pause_text("Fortsetzen")
        else:
            self.unpaused = True
            self.button_win.set_pause_text("Pause")
        return

    def back(self):
        """
        Aufruf bei Klicken des Buttons "Zurück" aus "ButtonWindowFrame". Bei Aufruf wird  ein Zeitschritt zurück
         gegangen. Dazu werden Warteschlangen entsprechend angepasst, wobei dazu virtuelle und reale Spieler gesondert
         betrachtet werden. Weiter werden Spieler, die bereits an ihrem Ziel angekommen sind und dieses nun wieder
         verlassen, wieder sichtbar gemacht. Anschließend wird die Funktion "self.runback()" aufgerufen, die alle
         Positionen neu berechnet und entsprechend anpasst.
        :return: Kein Rückgabewert
        """
        positionenV = []  # für graphischen Teil
        red_list = []  # für graphischen Teil
        for e in self.E:
            positionenV.append([])
            # füge Spieler wieder zu Queue hinzu
            for out in reversed(self.leaveQueue[self.zeitpunkt -1][self.E.index(e)]):
                try:
                    int(out)  # TypeError, falls "out" KEIN virtueller Spieler
                    self.deques[self.E.index(e)].append(out)  # virtuelle Spieler wieder zur Queue hinzufügen
                    positionenV[self.E.index(e)].append(len(self.deques[self.E.index(e)]) -1)  # merke Positionen
                    # virtueller Spieler in Warteschlange -> werden an "self.button_win" übergeben
                    continue
                except TypeError:  # nicht-virtueller Spieler
                    position = self.pos(out[0], self.zeitpunkt -1)
                    if position not in self.V:  # prüfe, ob Spieler zuvor in Warteschlange war
                        # falls ja, wird er wieder hinzugefügt
                        self.deques[self.E.index(position[0])].append(out)
                        red_list.append(out)  # merke Spieler, deren Farbe auf rot gesetzt werden soll

        self.button_win.restore_players(positionenV, red_list)  # graphische Umsetzung des Obigen

        # speichere Spieler und deren neue Farbe in Liste, zur Übergabe an 'self.button_win'
        player_recolor = []
        t_index = []
        for e in self.E:
            while True:
                try:
                    # IndexError falls deque leer. Entferne Spieler aus Queue, falls sie diese zum Zeitpunkt
                    # "self.zeitpunkt" betreten haben
                    latest = self.deques[self.E.index(e)].popleft()
                    # TypeError, falls "latest" KEIN virtueller Spieler ist, ansonsten: prüft, ob virtueller Spieler aus
                    # Plot und "self.button_win.spielerV" entfernt werden muss
                    if int(latest) == self.zeitpunkt or int(latest) == self.zeitpunkt -1:
                        # entferne Spieler auch aus "self.button_win.spielerV" und aus Plot
                        self.button_win.remove_artist(self.E.index(e))
                        continue
                    # latest ist virtueller Spieler und schon vor "self.zeitpunkt" in deque gewesen
                    else:
                        # "latest" wird wieder zur deque hinzugefügt und es wird abgebrochen
                        self.deques[self.E.index(e)].appendleft(latest)
                        break
                except IndexError:  # deque leer
                    break
                except TypeError:  # "latest" ist realer Spieler
                    # prüft, ob Spieler "latest" sich im Startknoten von "e" befindet, also die deque genau zu
                    # "self.zeitpunkt" betreten hat und entfernt werden muss
                    if self.pos(latest[0], self.zeitpunkt) == e[0]:
                        continue
                    # prüft, ob Spieler "latest" Kante "e" zum Zeitpunkt "self.zeitpunkt-1" betreten hat
                    elif self.fp[self.zeitpunkt-1][latest[0]] == e:
                        # wenn ja, umfärben von rot auf ursprüngliche Farbe und weiter mit nächstem Spieler
                        # merke Spieler
                        player_recolor.append(latest)
                        # merke Index des Zielknoten dieses Spielers, zur Bestimmung der Farbe
                        t_index.append(self.ST[latest[0]][1]-1)
                        continue
                    else:
                        self.deques[self.E.index(e)].appendleft(latest)  # wenn nein, Abbruch
                        break

        # umfärben
        self.button_win.change_color(player_recolor, t_index)
        # Informationen für graphischen Teil
        visibility = []
        nex_config = False
        prev_config = False
        # prüfe, ob Spieler zum Zeitpunkt "self.zeitpunkt" ihren Zielknoten erreicht haben
        for i in self.I:
            if self.zeitpunkt == self.z[i]:
                self.ankunft -= 1  # aktualisiere Abbruchkriterium
                visibility.append(i)
        if self.button_win.nex['state'] == 'disabled':
            nex_config = True
            self.button_win.after(1000, self.runner)  # runner neustarten, da Programm bereits terminiert hat
        self.zeitpunkt -= 1
        if self.zeitpunkt == 0:
            prev_config = True

        # visuelle updates
        self.button_win.updates(visibility, nex_config, prev_config)

        self.runback(self.zeitpunkt)

        return

    def weiter(self):
        """
        Bei Aufruf wird ein Zeitschritt weiter gegangen, d.h. 'self.run()' aufgerufen.
        :return: Kein Rückgabewert
        """
        if self.ankunft >= self.num:
            self.button_win.set_nex_status("disabled")
            return 0
        if self.button_win.prev['state'] == 'disabled':
            self.button_win.set_prev_status("normal")
        # Vergrößere Liste "self.leaveQueue", falls nötig (also falls "self.run(self.zeitpunkt + 1)" vorher noch
        # nicht aufgerufen wurde)
        if len(self.leaveQueue) < self.zeitpunkt + 2:
            self.fm.append(self.num * [None])
            self.fp.append(self.num * [None])
            self.leaveQueue.append([])
        self.zeitpunkt += 1
        self.run(self.zeitpunkt)
        return

    def run(self, theta):
        """
        Führt einen Vorwärtsschritt aus, zeigt also die Positionen aller Spieler zum Zeitpunkt "theta". Dazu wird
        erst bestimmt, ob die Funktion "run" bereits mit genau diesem "theta" zuvor schon aufgerufen wurde
        ("rerun"=True), oder nicht ("rerun"=False). Im Falle eines reruns müssen nämlich die Distanzen und die Werte
        für f^+, f^- nicht erneut berechnet werden, da diese bereits im ersten Aufruf von "run(theta)" bestimmt wurden.
        In beiden Fällen werden die Warteschlangen an den Zeitpunkt "theta" angepasst, sowie der Plot aktualisiert.
        :param theta: aktueller Zeitpunkt
        :return: kein Rückgabewert
        """
        newDists = []
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
            for i in self.I:
                if self.z[i] == -1:
                    if self.currentPos[i] in self.V:
                        if theta > 0:
                            newDists.append(self.dijkstra(self.graphReversed, "t{}".format(self.ST[i][1]),
                                                          self.currentPos[i],i,theta -1,visited=[], distances={}))
                        else:
                            newDists.append(self.dijkstra(self.graphReversed, "t{}".format(self.ST[i][1]),
                                                          self.currentPos[i],i,theta,visited=[], distances={}))
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

        # speichere Spieler und deren neue Farbe in Liste, zur Übergabe an 'self.button_win'
        player_recolor = []
        t_index = []
        for e in self.E:  # bestimme aktuelle Queue für jede Kante
            if not rerun and theta > 0:
                self.leaveQueue[theta -1].append([])
            residualCap[self.E.index(e)] -= min(len(self.deques[self.E.index(e)]), self.nu[self.E.index(e)])
            for out in range(min(len(self.deques[self.E.index(e)]), self.nu[self.E.index(e)])):
                nextPlayer = self.deques[self.E.index(e)].pop()
                if not rerun:
                    self.leaveQueue[theta -1][self.E.index(e)].append(nextPlayer)
                try:
                    int(nextPlayer)  # TypeError, falls "nextPlayer" nicht-virtueller Spieler
                    self.button_win.pop_right(self.E.index(e))
                except TypeError:
                    # merke umzufärbende Spieler
                    player_recolor.append(nextPlayer)
                    # merke Nummer des Zielknoten dieses Spielers, zur Bestimmung der Farbe
                    t_index.append(self.ST[nextPlayer[0]][1]-1)
            if theta > 0:
                for i in self.I:
                    if self.fp[theta-1][i] == e:
                        if residualCap[self.E.index(e)] == 0:
                            self.deques[self.E.index(e)].appendleft((i, 'r'))  # hinzufügen Spieler "i" zu Queue , 'r'
                            # markiert diesen als realen Spieler
                        else:
                            if not rerun:
                                self.leaveQueue[theta-1][self.E.index(e)].append((i, 'r'))  # 'r': realer Spieler
                            residualCap[self.E.index(e)] -= 1

        # Umfärben in Plot
        self.button_win.change_color(player_recolor, t_index)

        if not rerun and theta > 0:
            for i in self.I:
                # verwendet den Wert von "self.currentPos[i]" vom VORHERIGEN Zeitpunkt um "fm" für den aktuellen
                # Zeitpunkt bestimmen zu können
                if self.currentPos[i] in self.V:
                    if self.fp[theta-1][i] is not None and self.r[self.E.index(self.fp[theta-1][i])] == 1:
                        for spielerliste in self.leaveQueue[theta-1][self.E.index(self.fp[theta-1][i])]:
                            try:
                                int(spielerliste)  # überspringe virtuelle Spieler
                            except TypeError:  # "spielerliste" repräsentiert realen Spieler, ist also Tupel der Form
                                # (j, 'r'), wobei j der Index des entsprechenden Spielers ist
                                if (i, 'r') == spielerliste:
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
                                except TypeError:  # "spielerliste" repräsentiert realen Spieler, ist also Tupel der
                                    # Form (j, 'r'), wobei j der Index des entsprechenden Spielers ist
                                    if (i, 'r') == spielerliste:
                                        self.fm[theta][i] = self.currentPos[i][0]

        for i in self.I:
            if self.z[i] == -1 or self.z[i] >= theta:
                self.currentPos[i] = self.pos(i, theta)
                # print("Position Spieler " + str(i+1) + " zum Zeitpunkt " + str(theta) + " : " +
                # str(self.currentPos[i]))
                if self.currentPos[i] == "t{}".format(self.ST[i][1]):
                    if self.z[i] == -1:
                        self.z[i] = theta  # setze Ankunftszeit
                    if self.z[i] == theta:  # aktualisiere Abbruchkriterium
                        self.ankunft += 1

        # Aktualisierung Plot
        self.button_win.draw_new_positions(self.currentPos)
        self.button_win.draw_new_queue_positions(self.deques)
        self.button_win.redraw(theta)

        if not rerun and theta > 0:
            self.c.append(self.varianten(self.alpha, self.G, var=self.variante, anz=J))

        return

    def runback(self, theta):
        """
        wie "run"-Methode zum Zeitpunkt "theta", angepasst für den Fall, dass "Zurück" - Button gedrückt wurde. Im
        Unterschied zu einem Aufruf von "run(theta)" mit "rerun"=True, werden hier die Warteschlangen nicht verändert,
        da dies bereits beim Klicken des Buttons "Zurück" in der Funktion "self.back()" geschieht.
        :param theta: aktueller Zeitpunkt
        :return: kein Rückgabewert
        """

        # Aktualisierung Positionen
        for i in self.I:
            if self.z[i] == -1 or self.z[i] > theta:
                self.currentPos[i] = self.pos(i, theta)
                # print("Position Spieler " + str(i+1) + " zum Zeitpunkt " + str(theta) + " : " + str(currentPos[i]))

        # Aktualisierung Plot
        self.button_win.draw_new_positions(self.currentPos)
        self.button_win.draw_new_queue_positions(self.deques)
        self.button_win.redraw(theta)

        return
