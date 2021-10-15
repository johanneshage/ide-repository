import numpy as np
from Graphics.OutputTable import OutputTable


class ContApp:
    """
    IDE - Konstruktion für allgemeine single-sink Netzwerke
    """

    def __init__(self, G, u):
        """

        :param G:
        :param u: Liste, welche für jeden Knoten aus 'G' eine leere Liste (-> kein externer Einfluss in diesen Knoten)
                  oder eine Liste von 2-Tupeln der folgenden Form enthält:
                  (a_i,x_i): beschreibt den Einfluss von x_i Flusseinheiten innerhalb des Zeitintervalls [a_i,a_{i+1})
                  in den entsprechenden Knoten. Beachte: nach Voraussetzung sind alle Einflüsse endlich, d.h. für jeden
                  dieser Knoten ist der letzte Eintrag in dieser Liste ein Tupel der Form (t, 0) für einen Zeitpunkt
                  t < infty.
        """
        self.E = []  # Kanten
        self.nu = []  # Kapazitaeten
        self.r = []  # Reisezeiten
        self.c = []  # Kosten für alle Kanten zu Beginn jeder Phase
        self.labels = []  # Knotenlabels
        self.V = list(G.keys())  # Liste der Knoten
        self.n = len(self.V)  # Anzahl Knoten
        self.G = G
        self.u = u
        self.items = G.items()
        self.keys = G.keys()
        self.eps = 10**(-8)  # Für Rundungsfehler
        self.flow_vol = []  # merke Flusswerte in den einzelnen Knoten für OutputTable

        for delta in self.items:
            for w in list(delta[1].keys()):
                self.E.append((delta[0], w))  # speichere alle Kanten in E
                self.r.append(list(self.items)[list(self.keys).index(delta[0])][1][w][1])  # speichere Reisezeiten in r
                self.nu.append(list(self.items)[list(self.keys).index(delta[0])][1][w][0])  # Kapazitäten in nu
        self.m = len(self.E)  # Anz. Kanten
        self.global_phase = [0]  # speichere Startzeitpunkte der "globalen" Phasen, also alle Zeitpunkte, zu denen sich
        # mindestens eine lokale Phase ändert
        self.q_global = [self.m * [0]]  # speichere für jede globale Phase alle momentanen Warteschlangenlängen zu
        # Beginn der Phase

        self.fp = []  # f^+
        self.fp_ind = []  # Liste der Indizes der Phasen
        self.fm = []  # f^-
        for i in range(self.m):
            self.fp.append([(0, 0)])
            self.fp_ind.append([])
            self.fm.append([(0, 0)])

        self.c.append(np.zeros(self.m))
        for e in self.E:  # Initialisierung "self.c" (Kosten)
            self.c[0][self.E.index(e)] = self.r[self.E.index(e)]

        # Zeitpunkte, zu denen sich der Zufluss in mindestens einem Quellknoten ändert
        self.u_start = set()
        for s_list in self.u:
            for t in s_list:
                self.u_start.add(t[0])

        self.graphReversed = self.reverse_graph(G)

        self.E_active = [np.ones(self.m)]
        self.labels.append(self.dijkstra(self.graphReversed, "t1", 0,visited=[], distances={}))

        self.E_active = [np.zeros(self.m)]
        for v in self.V:
            v_ind = self.V.index(v)
            outneighbors = self.G[v].keys()
            for w in outneighbors:
                w_ind = self.V.index(w)
                edge = self.E.index((v, w))
                if abs(self.labels[0][v_ind] - self.labels[0][w_ind] - self.c[0][edge]) < self.eps:
                    self.E_active[0][edge] = 1

        self.del_plus_label = np.zeros(self.n)
        self.main()

    # Quelle:
    # https://www.geeksforgeeks.org/topological-sorting/#:~:text=Topological%20sorting%20for%20Directed%20Acyclic,4%202%203%201%200%E2%80%9D
    # A recursive function used by topologicalSort
    def topologicalSortUtil(self, v, visited, stack):

        # Mark the current node as visited.
        visited[self.V.index(v)] = True

        # Recur for all the vertices adjacent to this vertex
        for w in self.G[v]:
            i = self.V.index(w)
            if not visited[i]:
                self.topologicalSortUtil(w, visited, stack)

        # Push current vertex to stack which stores result
        stack.append(v)

    # The function to do Topological Sort. It uses recursive
    # topologicalSortUtil()
    def topologicalSort(self):
        # Mark all the vertices as not visited
        visited = [False]*self.n
        visited[self.V.index('t1')] = True
        stack = ['t1']

        # Call the recursive helper function to store Topological
        # Sort starting from all vertices one by one
        for i in range(self.n):
            if not visited[i]:
                self.topologicalSortUtil(self.V[i], visited, stack)
        return stack

    # Analog zu 'topologicalSortUtil', verwendet jedoch nur aktive Kanten (Kanten aus 'self.E_active')
    def topologicalSortUtilActive(self, v, visited, stack):
        # Mark the current node as visited.
        visited[self.V.index(v)] = True

        # Recur for all the vertices adjacent to this vertex
        outneighbors = self.G[v]
        for w in outneighbors:
            if self.E_active[-1][self.E.index((v,w))] and not visited[self.V.index(w)]:
                self.topologicalSortUtilActive(w, visited, stack)

        # Push current vertex to stack which stores result
        stack.append(v)

    # The function to do Topological Sort. It uses recursive topologicalSortUtil()
    # Analog zu 'topologicalSort', verwendet jedoch nur aktive Kanten (Kanten aus 'self.E_active')
    def topologicalSortActive(self):
        # Mark all the vertices as not visited
        visited = [False]*self.n
        visited[self.V.index('t1')] = True
        stack = ['t1']

        # Call the recursive helper function to store Topological
        # Sort starting from all vertices one by one
        for i in range(self.n):
            if not visited[i]:
                self.topologicalSortUtilActive(self.V[i], visited, stack)
        return stack

    # Quelle: http://www.gilles-bertrand.com/2014/03/dijkstra-algorithm-python-example-source-code-shortest-path.html
    def dijkstra(self, graph, src, phase_ind, visited=[], distances={}, predecessors={}):
        """
        Berechnet rekursiv kürzesten Weg und dessen Kosten von "src" zu jedem erreichbaren Knoten in "graph" zum
        Zeitpunkt "self.global_phase[phase_ind]"
        :param graph: Graph als Dictionary
        :param src: Startknoten
        :param phase_ind: Index des Zeitpunkts
        :param visited: Liste der bereits besuchten Knoten, anfangs leer, muss nicht beachtet werden, da nur intern für
         die Funktion benötigt
        :param distances: Liste der bereits berechneten Distanzen aller Knoten zu "src", anfangs leer, muss nicht
         beachtet werden, da nur intern für die Funktion benötigt
        :param predecessors: Liste der bereits ermittelten Vorfahren aller Knoten, anfangs leer, muss nicht beachtet
         werden, da nur intern für die Funktion benötigt
        :return: "output": Liste über alle Knoten mit deren Distanzen zur ersten(!) "src", also zu dem Knoten, der beim
                           ersten Aufruf von "dijkstra" als "src" übergeben wurde
                 "self.dijkstra(graph, x, theta, visited, distances, predecessors)": sorgt für rekursives
                           Aufrufen dieser Funktion, endet mit return von "output"
        """
        """ calculates a shortest path tree routed in src
        """
        # a few sanity checks
        if src not in graph:
            raise TypeError('The root of the shortest path tree cannot be found')
        # if it is the initial run, initializes the cost
        if not visited:
            distances[src] = 0
        # visit the neighbors
        for neighbor in graph[src]:
            if neighbor not in visited:
                new_distance = distances[src] + self.c[phase_ind][self.E.index((neighbor,src))]
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
        vals = unvisited.values()
        # check if no more nodes are reachable
        if np.all(len(vals) == 0 or vals == float('inf')):
            output = []
            for v in self.V:
                if v in distances:
                    output.append(distances[v])
                else:
                    output.append(float('inf'))
            return output
        x=min(unvisited, key=unvisited.get)
        return self.dijkstra(graph, x, phase_ind, visited, distances, predecessors)

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

    def waterfilling_algo(self, v, b, outneighbors, w_slope):
        """
        Berechnet Flussaufteilung des in 'v' vorhandenen Flussvolumens auf die ausgehenden aktiven Kanten
        :param v: Knoten
        :param b: Flussmenge im Knoten 'v'
        :param outneighbors: Liste aller Knoten, die über eine direkte aktive Kante von 'v' aus erreichbar sind
        :param w_slope: Liste über Änderung der Labels aller Knoten in 'outneighbors'
        :return: Flussaufteilung 'z' auf alle von 'v' ausgehenden aktiven Kanten in Form einer Liste
        """
        pk = len(outneighbors)
        if pk == 0:
            return []
        alpha = np.zeros(pk)
        beta = np.copy(w_slope)
        gamma = np.zeros(pk)
        old_sorting = np.zeros(pk)

        # setze 'beta' - Werte
        for i in range(pk):
            e_ind = self.E.index((v, outneighbors[i]))
            if self.q_global[-1][e_ind] > 0:
                beta[i] -= 1
            old_sorting[i] += i

        # sortiere 'beta' und 'outneighbors', sodass 'beta' - Werte aufsteigend
        zipped = sorted(zip(beta, outneighbors, old_sorting))
        # beta, outneighbors, old_sorting = zipped[0], zipped[1], zipped[2]
        beta = [element for element, _, _ in zipped]
        outneighbors = [element for _, element, _ in zipped]
        old_sorting = [element for _, _, element in zipped]

        # setze 'alpha'- und 'gamma' - Werte (mit der neuen Sortierung)
        for i in range(pk):
            e_ind = self.E.index((v, outneighbors[i]))
            alpha[i] += self.nu[e_ind]
            if self.q_global[-1][e_ind] == 0:
                gamma[i] += self.nu[e_ind]

        # h = lambda i, z: beta[i] + np.max([0, 1.0/alpha[i] * (z - gamma[i])])

        # Beginn des Algorithmus
        z = np.zeros(pk)
        # Bestimmt max{z | h_i(z) <= beta_r}, falls r >= i
        find_max = lambda i, r: alpha[i] * (beta[r] - beta[i]) + gamma[i]
        r = 0
        for r in range(1, pk + 1):
            sum = 0
            for i in range(r):
                sum += find_max(i, r - 1)
            if sum > b:
                r -= 1
                break
        if r == pk + 1:
            r -= 1

        else_case = False
        if r < pk:
            sum = 0
            for i in range(r):
                sum += find_max(i, r)
            if sum <= b:
                z_sum = 0
                for i in range(r):
                    z[i] += find_max(i, r)
                    z_sum += z[i]
                # z_sum -= z[r - 1]
                z[r] += b - z_sum
            else:
                else_case = True
        else:
            else_case = True

        if else_case:
            z_sum = 0
            for i in range(r):
                z[i] += find_max(i, r - 1)
                z_sum += z[i]
            b_prime = b - z_sum
            if r == 1:
                z[0] += b_prime
            else:
                alpha_sum = 0
                # for j in range(r - 1):
                for j in range(r):
                    alpha_sum += alpha[j]
                for i in range(r):
                    z[i] += b_prime * (alpha[i] + 0.0)/alpha_sum

        # sortiere 'z' nach der ursprünglichen Sortierung
        zipped = sorted(zip(old_sorting, z))
        z = [element for _, element in zipped]
        return z

    def calc_b(self, v, phase):
        """
        Berechnet zum Zeitpunkt 'phase' im Knoten 'v' vorhandene Flussmenge b_v^- (phase)
        :param v: Knoten
        :param phase: Zeitpunkt
        :return: b_v^- (phase)
        """
        b = 0
        v_ind = self.V.index(v)
        in_paths = self.get_ingoing_edges(v)
        for e_ind in in_paths:
            ind = self.last_fm_change(e_ind, phase)
            if len(self.fm[e_ind]) > ind + 1 and abs(self.fm[e_ind][ind + 1][0] - phase) < self.eps:
                ind += 1
            b += self.fm[e_ind][ind][1]
        u_v = self.u[v_ind]
        u_v_len = len(u_v)
        for tuple_ind in range(u_v_len - 1, -1, -1):
            if u_v[tuple_ind][0] <= phase:
                b += u_v[tuple_ind][1]
                break
        # speichere b - Wert für OutputTable
        self.flow_vol[-1][v_ind] = b
        return b

    def get_ingoing_edges(self, v):
        """
        bestimmt alle in 'v' eingehenden Kanten
        :param v: Knoten
        :return: Liste der Indizes der Kanten
        """
        preds = self.graphReversed[v].keys()
        return [self.E.index((u,v)) for u in preds]

    def get_ingoing_active_edges(self, v):
        """
        bestimmt alle in 'v' eingehende, momentan aktive Kanten
        :param v: Knoten
        :return: Liste der Indizes der aktiven Kanten
        """
        preds = self.graphReversed[v].keys()
        delta_m = [self.E.index((u,v)) for u in preds]
        return [e for e in delta_m if self.E_active[-1][e]]

    def get_outgoing_edges(self, v):
        """
        bestimmt alle aus 'v' ausgehenden Kanten
        :param v: Knoten
        :return: Liste der Indizes der Kanten
        """
        return [self.E.index((v,u)) for u in self.G[v].keys()]

    def get_outgoing_active_edges(self, v):
        """
        bestimmt alle aus 'v' ausgehende, momentan aktive Kanten
        :param v: Knoten
        :return: Liste der aktiven Kanten
        """
        delta_p = [self.E.index((v,u)) for u in self.G[v].keys()]
        return [e for e in delta_p if self.E_active[-1][e]]

    def change_of_label(self, e, phase, z):
        """
        Gibt die momentane Änderung des Labels des Startknotens von 'e' unter Einfluss 'z' in 'e' an.
        :param e: Betrachtete Kante
        :param phase: Betrachteter Zeitpunkt
        :param z: Einflussrate in Kante 'e'
        :return: Änderungsrate des Labels des Startknotens von 'e'
        """
        e_ind = self.E.index(e)
        tar_ind = self.V.index(e[1])
        if self.q_global[self.global_phase.index(phase)][e_ind] > self.eps:
            return (z - self.nu[e_ind])/self.nu[e_ind] + self.del_plus_label[tar_ind]
        return np.max([(z - self.nu[e_ind])/self.nu[e_ind], 0]) + self.del_plus_label[tar_ind]

    def change_of_cost(self, e_ind, phase, z):
        """
        Gibt die momentane Kostenänderung der Kante mit Index 'e_ind' bei Zufluss 'z' an.
        :param e_ind: Index der betrachteten Kante
        :param phase: Betrachteter Zeitpunkt
        :param z: Einflussrate in Kante 'e'
        :return: Änderungsrate der Kosten
        """
        if self.q_global[self.global_phase.index(phase)][e_ind] > self.eps:
            return (z - self.nu[e_ind])/self.nu[e_ind]
        return np.max([(z - self.nu[e_ind])/self.nu[e_ind], 0])

    def last_fm_change(self, e_ind, theta):
        """
        Bestimmt den größten Zeitpunkt <= 'theta', zu dem sich der f^- -Wert von Kante 'e_ind' ändert
        :param e_ind: Kantenindex
        :param theta: aktueller Zeitpunkt
        :return: Index des gesuchten Zeitpunkts
        """
        fm_len = len(self.fm[e_ind])
        for t in range(fm_len - 1, 0, -1):
            if self.fm[e_ind][t][0] < theta + self.eps:
                return t
        return 0

    def main(self):
        """
        Hauppteil: Konstuiert schrittweise einen kontinuierlichen IDE-Fluss zu den Eingabedaten aus 'Math.cont_data'. Erzeugt Tabelle der
        Klasse 'Graphics.OutputTable' mit den entsprechenden Daten des Flusses.
        :return: 0
        """
        theta = 0
        # Obergrenze für theta
        T = 10000
        # Aufteilung des Flusses
        x_total = np.zeros(self.m)
        stop_outflow = []
        top_ord = self.topologicalSort()
        while theta < T:
            # in der Zukunft liegende Zeitpunkte aus der Liste 'self.u_start'
            start_points = [t for t in self.u_start if t > theta]
            # in der Zukunft liegende Zeitpunkte, zu denen der f^- -Wert von mindestens einer Kante auf 0 springt (wird
            # während der Laufzeit aktualisiert)
            stop_outflow = [t for t in stop_outflow if t > theta]
            top_ord_act = self.topologicalSortActive()
            theta_ind = self.global_phase.index(theta)
            # Liste aller Kanten, deren Warteschlange in der aktuellen Phase 0 wird, und deren f^- -Werte auf 0
            # gesetzt werden müssen
            fm_to_zero = []
            # 'next_phase' bestimmt am Ende der Schleife, wie lange aktuelle Phase andauert
            if len(start_points) > 0 or len(stop_outflow) > 0:
                next_phase = np.min(start_points + stop_outflow) - theta
            else:
                next_phase = T
            self.flow_vol.append(np.zeros(self.n))
            for v in top_ord:
                # Flussaufteilung des im Knoten 'v' vorhandenen Flussvolumens
                active_paths = self.get_outgoing_active_edges(v)
                active_neighbors = [self.E[e][1] for e in active_paths]
                w_slope = [self.del_plus_label[self.V.index(node)] for node in active_neighbors]
                x = self.waterfilling_algo(v, self.calc_b(v, theta), active_neighbors, w_slope)
                firstedge = True
                # betrachte aktive Kanten
                for e_ind in active_paths:
                    e = self.E[e_ind]
                    x_total[e_ind] = x[active_paths.index(e_ind)]

                    if theta == 0:
                        start_ind = self.V.index(e[0])
                        if len(self.u[start_ind]) > 0 and self.u[start_ind][0][0] == 0:
                            self.fp[e_ind][0] = (0, x_total[e_ind])
                            self.fp_ind[e_ind].append(0)
                    # falls sich f^+ -Wert in dieser Phase ändert, aktualisiere 'self.fp'
                    elif abs(self.fp[e_ind][-1][1] - x_total[e_ind]) > self.eps:
                        self.fp[e_ind].append((theta, x_total[e_ind]))
                        self.fp_ind[e_ind].append(theta_ind)
                    if self.q_global[theta_ind][e_ind] > self.eps:
                        outflow = self.nu[e_ind]
                    else:
                        outflow = np.min([x_total[e_ind], self.nu[e_ind]])
                    fm_ind = self.last_fm_change(e_ind, theta + self.r[e_ind])
                    # falls sich f^- -Wert durch die Flussaufteilung in dieser Phase ändert, aktualisiere 'self.fm'
                    if abs(self.fm[e_ind][fm_ind][1] - outflow) > self.eps:
                        self.fm[e_ind].append((theta + self.r[e_ind], outflow))

                    # bestimme, ob sich vor dem Ende der aktuellen Phase ein f^- -Wert ändert -> verkürze Phase
                    last_fm_ind = self.last_fm_change(e_ind, theta)
                    if len(self.fm[e_ind]) > last_fm_ind + 1 and \
                            self.eps < self.fm[e_ind][last_fm_ind + 1][0] - theta < next_phase - self.eps:
                        next_phase = self.fm[e_ind][last_fm_ind + 1][0] - theta
                        # wird zurückgesetzt, da durch Verkürzung der Phase der f^- -Wert erst in einer späteren Phase
                        # neu gesetzt wird
                        fm_to_zero = []

                    change_of_c = self.change_of_cost(e_ind, theta, x[active_paths.index(e_ind)])
                    change_of_q = change_of_c * self.nu[e_ind]
                    new_del_plus = change_of_c + self.del_plus_label[self.V.index(e[1])]
                    # prüfe, ob Änderung des Labels von Knoten 'v' angepasst werden muss
                    if firstedge or new_del_plus < self.del_plus_label[self.V.index(v)]:
                        self.del_plus_label[self.V.index(v)] = new_del_plus
                        firstedge = False
                    # Falls die Warteschlange von 'e' unter aktuellem Fluss abgebaut wird, bestimme Zeitpunkt, zu
                    # dem diese vollständig abgebaut ist (bei gleich bleibendem Fluss)
                    if change_of_q < -self.eps:
                        # 'phase_length': Dauer bis Warteschlangenlänge gleich 0
                        phase_length = - self.q_global[theta_ind][e_ind] / change_of_q
                        if phase_length < next_phase - self.eps:
                            next_phase = phase_length
                            # prüfe ob Zufluss 0 und somit 'fm' auf 0 gesetzt werden muss
                            if change_of_q + self.nu[e_ind] < self.eps:
                                fm_to_zero = [e_ind]
                        # phase_length == next_phase
                        elif max([abs(phase_length - next_phase), change_of_q + self.nu[e_ind]]) < self.eps:
                            fm_to_zero.append(e_ind)

            # erneuter Durchlauf der vorherigen Schleife, diesmal werden die inaktiven Kanten betrachtet. Diese Aufteilung ist notwendig,
            # damit die Änderungen der Knotenlabels erst richtig gesetzt (siehe Schleife zuvor), und dann für weitere Rechnungen
            # (siehe nachstehende Schleife) verwendet werden.
            for v in top_ord:
                active_paths = self.get_outgoing_active_edges(v)
                delta_p = self.get_outgoing_edges(v)
                inactive_paths = [e for e in delta_p if e not in active_paths]
                for e_ind in inactive_paths:
                    x_total[e_ind] = 0
                    # falls f^+ -Wert vorher > 0 war, so wird dieser hier auf 0 gesetzt, da Kante inaktiv
                    if abs(self.fp[e_ind][-1][1]) > self.eps:
                        self.fp[e_ind].append((theta, 0))
                        self.fp_ind[e_ind].append(theta_ind)

                    if self.q_global[theta_ind][e_ind] > self.eps:
                        outflow = self.nu[e_ind]
                    else:
                        outflow = 0
                    fm_ind = self.last_fm_change(e_ind, theta + self.r[e_ind])
                    # falls sich f^- -Wert durch die Flussaufteilung in dieser Phase ändert, aktualisiere 'self.fm'
                    if abs(self.fm[e_ind][fm_ind][1] - outflow) > self.eps:
                        self.fm[e_ind].append((theta + self.r[e_ind], outflow))

                    # bestimme, ob sich vor dem Ende der aktuellen Phase ein f^- -Wert ändert -> verkürze Phase
                    last_fm_ind = self.last_fm_change(e_ind, theta)
                    if len(self.fm[e_ind]) > last_fm_ind + 1 and \
                            self.eps < self.fm[e_ind][last_fm_ind + 1][0] - theta < next_phase - self.eps:
                        next_phase = self.fm[e_ind][last_fm_ind + 1][0] - theta
                        # wird zurückgesetzt, da durch Verkürzung der Phase der f^- -Wert erst in einer späteren Phase
                        # neu gesetzt wird
                        fm_to_zero = []

                    change_of_q = self.change_of_cost(e_ind, theta, 0) * self.nu[e_ind]
                    # Falls Warteschlange existiert (und somit abgebaut wird da inaktiv), bestimme Zeitpunkt, zu dem
                    # diese vollständig abgebaut ist
                    if change_of_q < -self.eps:
                        phase_length = -self.q_global[theta_ind][e_ind] / change_of_q
                        # phase_length < next_phase
                        if phase_length < next_phase - self.eps:
                            next_phase = phase_length
                            # prüfe ob Zufluss 0 und somit 'fm' auf 0 gesetzt werden muss
                            if change_of_q + self.nu[e_ind] < self.eps:
                                fm_to_zero = [e_ind]
                        # phase_length == next_phase
                        elif max([abs(phase_length - next_phase), change_of_q + self.nu[e_ind]]) < self.eps:
                            fm_to_zero.append(e_ind)

                len_act = len(active_paths)
                if len_act > 0:
                    for i in range(len_act):
                        active_ind = active_paths[i]
                        if x_total[active_ind] > 0 or i == len_act - 1:
                            # active_ind = active_paths[i]  # aktive Kante, die während der gesamten Phase aktiv bleibt
                            active_change = self.change_of_cost(active_ind, theta, x_total[active_ind])  # Änderung der Kosten dieser Kante
                            break

                    for e_ind in inactive_paths:
                        change = self.change_of_cost(e_ind, theta, 0)
                        # prüfe, wann inaktive Kanten unter momentanem Einfluss aktiv werden
                        tar_ind = self.V.index(self.E[e_ind][1])
                        act_ind = self.V.index(self.E[active_ind][1])
                        if self.labels[theta_ind][tar_ind] + self.q_global[-1][e_ind]/self.nu[e_ind] + self.r[e_ind] + \
                                (change + self.del_plus_label[tar_ind]) * next_phase < \
                                self.labels[theta_ind][act_ind] + self.q_global[-1][active_ind]/self.nu[active_ind] + \
                                self.r[active_ind] + (active_change + self.del_plus_label[act_ind]) * next_phase and \
                                abs(change + self.del_plus_label[tar_ind] - active_change - self.del_plus_label[act_ind])\
                                > self.eps:
                            time_ub = np.abs((self.labels[theta_ind][tar_ind] + self.q_global[-1][e_ind]/self.nu[e_ind]
                                              + self.r[e_ind] - self.labels[theta_ind][act_ind]
                                              - self.q_global[-1][active_ind]/self.nu[active_ind] - self.r[active_ind])
                                             / (active_change + self.del_plus_label[act_ind] - change - self.del_plus_label[tar_ind]))
                            if time_ub < next_phase:
                                next_phase = time_ub
                                # wird zurückgesetzt, da durch Verkürzung der Phase f^- -Wert erst in einer späteren
                                # Phase neu gesetzt wird
                                fm_to_zero = []

            if next_phase != T:
                # aktualisiere Warteschlangenlängen und Kosten
                new_q_global = []
                self.c.append(np.zeros(self.m))
                for e_ind in range(self.m):
                    next_q_len = self.q_global[theta_ind][e_ind] + \
                                 self.change_of_cost(e_ind, theta, x_total[e_ind]) * self.nu[e_ind] * next_phase
                    if next_q_len < self.eps:
                        next_q_len = 0
                    new_q_global.append(next_q_len)
                    self.c[-1][e_ind] = new_q_global[-1] / self.nu[e_ind] + self.r[e_ind]
                # speichere aktuelle Warteschlangenlängen
                self.q_global.append(new_q_global)

            print("Kanten mit positivem Einfluss zum Zeitpunkt", theta, ":")
            for v_ind in range(self.n):
                out_neighbors = self.get_outgoing_edges(top_ord_act[v_ind])
                for e_ind in out_neighbors:
                    if x_total[e_ind] > 0:
                        print(self.E[e_ind], x_total[e_ind])

            theta += next_phase
            if next_phase != T:
                # speichere Phase
                self.global_phase.append(theta)
                self.labels.append(self.dijkstra(self.graphReversed, "t1", len(self.global_phase) - 1,
                                                 visited=[], distances={}))

                self.E_active.append(np.zeros(self.m))
                for v in self.V:
                    v_ind = self.V.index(v)
                    outneighbors = self.G[v].keys()
                    for w in outneighbors:
                        w_ind = self.V.index(w)
                        edge = self.E.index((v, w))
                        if abs(self.labels[-1][v_ind] - self.labels[-1][w_ind] - self.c[-1][edge]) < self.eps:
                            self.E_active[-1][edge] = 1

                for e_ind in fm_to_zero:
                    self.fm[e_ind].append((theta + self.r[e_ind], 0))
                    stop_outflow.append(theta + self.r[e_ind])

        # am Ende sind alle f^+ -Werte 0
        for e in range(self.m):
            if abs(self.fp[e][-1][1]) > self.eps:
                self.fp[e].append((theta - next_phase, 0))
                self.fp_ind[e].append(theta_ind)

        # erzeuge Ausgabe
        OutputTable(self.V, self.E, self.nu, self.fp, self.fp_ind, self.fm, self.q_global, self.global_phase, self.c, self.labels,
                    self.flow_vol)
        return 0

