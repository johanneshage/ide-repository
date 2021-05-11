import numpy as np
import scipy.integrate as integrate
from scipy.optimize import minimize, LinearConstraint, Bounds



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
        self.fm = []  # f^-
        for i in range(self.m):
            self.fp.append([(0,0)])
            self.fm.append([(0,0)])

        self.c.append(np.zeros(self.m))
        for e in self.E:  # Initialisierung "self.c" (Kosten)
            self.c[0][self.E.index(e)] = self.r[self.E.index(e)]

        self.u_start = set()
        for s_list in self.u:
            for t in s_list:
                self.u_start.add(t[0])

        self.graphReversed = self.reverse_graph(G)

        self.E_active = np.ones(self.m)
        self.labels.append(self.dijkstra(self.graphReversed, "t1", 0,visited=[], distances={}))

        self.E_active = np.zeros(self.m)
        for v in self.V:
            for w in self.G[v].keys():
                v_ind = self.V.index(v)
                w_ind = self.V.index(w)
                edge = self.E.index((v, w))
                if self.labels[0][v_ind] == self.labels[0][w_ind] + self.c[0][edge]:
                    self.E_active[edge] = 1

        self.del_plus_label = np.zeros(self.n)
        self.main()

    # Quelle: https://www.geeksforgeeks.org/topological-sorting/#:~:text=Topological%20sorting%20for%20Directed%20Acyclic,4%202%203%201%200%E2%80%9D.
    # A recursive function used by topologicalSort
    def topologicalSortUtil(self, v, visited, stack):
        # Mark the current node as visited.
        visited[self.V.index(v)] = True

        # Recur for all the vertices adjacent to this vertex
        outneighbors = self.G[v]
        for w in outneighbors:
            if self.E_active[self.E.index((v,w))] and not visited[self.V.index(w)]:
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

    def q(self, e, global_theta, epsilon):
        """
        berechnet Warteschlangenlänge für Kante 'e' zum Zeitpunkt 'global_theta' + 'epsilon', unter der Voraussetzung,
        dass Einflussrate zwischen den globalen Zeitpunkten linear
        :param e: Kantenindex
        :param global_theta: Startzeitpunkt der aktuellen globalen Phase
        :param epsilon: Zeit seit beginn der aktuellen globalen Phase (insbesondere klein genug, dass keine neue
        globale Phase begonnen hat)
        :return: Länge der Warteschlange
        """
        if epsilon == 0:
            return self.q_global[self.global_phase.index(global_theta)][e]

        if self.q_global[self.global_phase.index(global_theta)][e] == 0:
            if self.fp[e][global_theta] <= self.nu[e]:
                return 0
        # beachte: folgender Wert ist größer gleich 0, da sonst 'epsilon' ungültig, d.h. 'epsilon' ist so groß, dass
        # 'global_theta' + 'epsilon' bereits in neuer globaler Phase liegt
        return self.q_global[self.global_phase.index(global_theta)][e] + epsilon*(self.fp[e][global_theta] - self.nu[e])

    def opt(self, v, phase):
        """
        Optimierungsproblem für Knoten 'v' zur Bestimmung der Verteilung des in der aktuellen 'phase' eintreffenden
        Flusses.
        :param v:
        :param phase:
        :return: Vektor x_e mit Eintrag für jede von 'v' ausgehende aktive Kante
        """
        neighbors = self.G[v].keys()
        delta_p = [self.E.index((v,w)) for w in neighbors]
        active_paths = [self.E[e] for e in delta_p if self.E_active[e]]
        if len(active_paths) == 0:
            return 0
        x0 = np.zeros(len(active_paths))
        bound = Bounds(np.zeros(len(active_paths)), len(active_paths) * [np.inf])

        def constraint(x):
            """
            Berechnung des Vektors b_v^-(phase) und Nebenbedingung \sum_{i=1}^{p_k} x_{vw_i} = b_v^-(phase) für
            Optimierungsproblem OPT-b_v^-(phase)
            :param x: Optimierungsvariable
            :param v: Knoten
            :param phase: Zeit
            :return: Gleichheits - Nebenbedingung
            """
            b = 0
            active_paths = self.get_ingoing_active_edges(v)
            for e in active_paths:
                e_ind = self.E.index(e)
                phase_ind = 0
                if phase > self.fm[e_ind][-1][0]:
                    phase_ind = len(self.fm[e_ind]) - 1
                else:
                    while phase > self.fm[e_ind][phase_ind][0]:
                        phase_ind += 1
                    if phase < self.fm[e_ind][phase_ind][0]:
                        phase_ind -= 1
                b += self.fm[e_ind][phase_ind][1]
            u_v = self.u[self.V.index(v)]
            u_v_len = len(u_v)
            for tuple_ind in range(u_v_len - 1, -1, -1):
                if u_v[tuple_ind][0] <= phase:
                    b += u_v[tuple_ind][1]
                    break
            return sum(x[:]) - b
            #return LinearConstraint(np.ones(len(x)), b, b)

        con = {'type':'eq', 'fun': constraint}

        def objective(x):
            return sum([integrate.quad(lambda z: self.change_of_label(e, phase, z), 0, x[active_paths.index(e)])[0] for
                        e in active_paths])

        return minimize(objective, x0, method='SLSQP',bounds=bound, constraints=con)

    def get_ingoing_edges(self, v):
        """
        bestimmt alle in 'v' eingehenden Kanten
        :param v: Knoten
        :return: Liste der Kanten
        """
        preds = self.graphReversed[v].keys()
        return [(u,v) for u in preds]

    def get_ingoing_active_edges(self, v):
        """
        bestimmt alle in 'v' eingehende, momentan aktive Kanten
        :param v: Knoten
        :return: Liste der aktiven Kanten
        """
        preds = self.graphReversed[v].keys()
        delta_m = [self.E.index((u,v)) for u in preds]
        return [self.E[e] for e in delta_m if self.E_active[e]]

    def get_outgoing_edges(self, v):
        """
        bestimmt alle aus 'v' ausgehenden Kanten
        :param v: Knoten
        :return: Liste der Kanten
        """
        return [(v,u) for u in self.G[v].keys()]

    def get_outgoing_active_edges(self, v):
        """
        bestimmt alle aus 'v' ausgehende, momentan aktive Kanten
        :param v: Knoten
        :return: Liste der aktiven Kanten
        """
        delta_p = [self.E.index((v,u)) for u in self.G[v].keys()]
        return [self.E[e] for e in delta_p if self.E_active[e]]

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
        if self.q_global[self.global_phase.index(phase)][e_ind] > 0:
            return (z - self.nu[e_ind])/self.nu[e_ind] + self.del_plus_label[tar_ind]
        return np.max([(z - self.nu[e_ind])/self.nu[e_ind], 0]) + self.del_plus_label[tar_ind]

    def change_of_cost(self, e, phase, z):
        """
        Gibt die momentane Änderung der Kostenänderung der Kante 'e' bei Zufluss 'z' an.
        :param e: Betrachtete Kante
        :param phase: Betrachteter Zeitpunkt
        :param z: Einflussrate in Kante 'e'
        :return: Änderungsrate der Kosten
        """
        e_ind = self.E.index(e)
        if self.q_global[self.global_phase.index(phase)][e_ind] > 0:
            return (z - self.nu[e_ind])/self.nu[e_ind]
        return np.max([(z - self.nu[e_ind])/self.nu[e_ind], 0])

    def main(self):
        theta = 0
        top_ord = self.topologicalSort()
        # Aufteilung des Flusses
        x_total = np.zeros(self.m)
        T = 100
        while theta < T:
            start_points = [t for t in self.u_start if t > theta]
            # 'next_phase' bestimmt am Ende der Schleife, wie lange aktuelle Phase andauert
            if len(start_points) > 0:
                next_phase = np.min(start_points)
            else:
                next_phase = T
            for v in top_ord:
                x = self.opt(v, theta)
                if x == 0:
                    continue
                active_paths = self.get_outgoing_active_edges(v)
                for e in active_paths:
                    e_ind = self.E.index(e)
                    e_act_ind = active_paths.index(e)
                    x_total[e_ind] = x.x[e_act_ind]

                    if theta == 0:
                        start_ind = self.V.index(self.E[e_ind][0])
                        if len(self.u[start_ind]) > 0 and self.u[start_ind][0][0] == 0:
                            self.fp[e_ind][0] = (0, x_total[e_ind])
                    if self.fp[e_ind][-1][1] != x_total[e_ind]:
                        self.fp[e_ind].append((theta, x_total[e_ind]))
                    outflow = np.min([x_total[e_ind], self.nu[e_ind]])
                    if self.fm[e_ind][-1][1] != outflow:
                        self.fm[e_ind].append((theta + self.r[e_ind], outflow))

                    change = self.change_of_cost(e, theta, x.x[e_act_ind])
                    self.del_plus_label[self.V.index(e[0])] = change + self.del_plus_label[self.V.index(e[1])]
                    # Falls die Warteschlange von 'e' unter aktuellem Fluss abgebaut wird, bestimme Zeitpunkt, zu dem
                    # diese vollständig abgebaut ist (bei gleich bleibendem Fluss)
                    if change < 0:
                        # 'phase_length': Dauer bis Warteschlangenlänge gleich 0
                        phase_length = - self.q_global[self.global_phase.index(theta)][self.E.index(e)] / change
                        if phase_length < next_phase:
                            next_phase = phase_length

                delta_p = self.get_outgoing_edges(v)
                inactive_paths = [e for e in delta_p if e not in active_paths]
                active = active_paths[0]  # aktive Kante
                active_change = self.change_of_cost(active, theta, x.x[0])  # Änderung der Kosten dieser Kante
                for e in inactive_paths:
                    e_ind = self.E.index(e)
                    if self.fp[e_ind][-1][1] != 0:
                        self.fp[e_ind].append((theta, 0))

                    change = self.change_of_cost(e, theta, 0)
                    # Falls Warteschlange existiert (und somit abgebaut wird da inaktiv), bestimme Zeitpunkt, zu dem
                    # diese vollständig abgebaut ist
                    if change < 0:
                        phase_length = -self.q_global[self.global_phase.index(theta)][self.E.index(e)]/change
                        if phase_length < next_phase:
                            next_phase = phase_length
                    # prüfe, wann inaktive Kanten unter momentanem Einfluss aktiv werden
                    tar_ind = self.V.index(e[1])
                    act_ind = self.V.index(active[1])
                    theta_ind = self.global_phase.index(theta)
                    if self.labels[theta_ind][tar_ind] + change * next_phase < self.labels[theta_ind][act_ind] +\
                            active_change * next_phase and change != active_change:
                        time_ub = np.abs((self.labels[theta_ind][act_ind] - self.labels[theta_ind][tar_ind]) /
                                        (change - active_change))
                        if time_ub < next_phase:
                            next_phase = time_ub

            new_q_global = []
            self.c.append(np.zeros(self.m))
            theta_ind = self.global_phase.index(theta)
            for e in self.E:
                ind = self.E.index(e)
                new_q_global.append(self.q_global[theta_ind][ind] + self.change_of_cost(e, theta, x_total[ind]) *
                                    self.nu[ind] * next_phase)
                # überprüfe, ob Warteschlange von 'e' in dieser Phase vollständig abgebaut wird und 'e' inaktiv ist,
                # dann muss f^-(e) auf 0 gesetzt werden
                if self.q_global[theta_ind][ind] > 0 and new_q_global[-1] == 0 and x_total[ind] == 0:
                    self.fm[ind].append((theta + next_phase, 0))
                self.c[-1][ind] = new_q_global[-1] / self.nu[ind] + self.r[ind]
            # speichere aktuelle Warteschlangenlängen
            self.q_global.append(new_q_global)

            print("Kanten mit positivem Einfluss zum Zeitpunkt", theta, " :")
            for v_ind in range(self.n):
                out_neighbors = self.get_outgoing_edges(top_ord[v_ind])
                for e in out_neighbors:
                    e_ind = self.E.index(e)
                    if x_total[e_ind] > 0:
                        print(e, x_total[self.E.index(e)])

            theta += next_phase
            # speichere Phase
            self.global_phase.append(theta)
            self.labels.append(self.dijkstra(self.graphReversed, "t1", len(self.global_phase) - 1,
                                             visited=[], distances={}))
        return 0

