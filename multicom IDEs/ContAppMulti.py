import numpy as np
# from OutputTableMulti import OutputTableMulti
import scipy.sparse
import copy
from collections import defaultdict
import time
# import pickle
# import matplotlib as mpl
# import matplotlib.pyplot  as  plt
# import os


class ContAppMulti:
    """
    IDE - Konstruktion für allgemeine single-sink Netzwerke
    """

    def __init__(self, G, u):
        """
        :param G:
        :param u: Liste, welche für jede Senke eine Liste von Einflussraten für jede Quelle enthält. Es gilt also: Der
                  Eintrag u[j][i] ist eine Liste mit den Einflussraten in Knoten s_{i+1} von Gut j+1. Diese Liste besteht für
                  ein festes Quelle -Senke -Paar aus Tupeln der folgenden Form:
                  (a_i,x_i): beschreibt den Einfluss von x_i Flusseinheiten innerhalb des Zeitintervalls [a_i,a_{i+1})
                  in den entsprechenden Knoten. Beachte: nach Voraussetzung sind alle Einflüsse endlich, d.h. für jeden
                  dieser Knoten ist der letzte Eintrag in dieser Liste ein Tupel der Form (t, 0) für einen Zeitpunkt
                  t < infty.
        """
        self.start_time = time.time()
        self.E = []  # Kanten
        self.nu = []  # Kapazitaeten
        self.r = []  # Reisezeiten
        self.G = G
        self.u = u
        self.I = len(self.u)
        self.labels = []  # Knotenlabels
        self.V = list(G.keys())  # Liste der Knoten
        self.n = len(self.V)  # Anzahl Knoten
        # self.b_ot = []
        self.b = []
        self.items = G.items()
        self.keys = G.keys()
        self.eps = 10**(-12)  # Für Rundungsfehler
        # self.eps = 10**(-5)
        self.bd_tol = 10**(-6)
        # self.flow_vol = []  # merke Flusswerte in den einzelnen Knoten für OutputTable # NUR FÜR OutputTable BENÖTIGT UND KOSTET SPEICHER
        self.delta_p = []

        for delta in self.items:
            for w in list(delta[1].keys()):
                self.E.append((delta[0], w))  # speichere alle Kanten in E
                self.r.append(list(self.items)[list(self.keys).index(delta[0])][1][w][1])  # speichere Reisezeiten in r
                self.nu.append(list(self.items)[list(self.keys).index(delta[0])][1][w][0])  # Kapazitäten in nu
        self.m = len(self.E)  # Anz. Kanten
        self.c = np.zeros(self.m)  # Kantenkosten für alle Kanten (Reisezeit + Warteschlangenlänge / Kapazität)
        self.global_phase = [0]  # speichere Startzeitpunkte der "globalen" Phasen, also alle Zeitpunkte, zu denen sich mindestens eine lokale Phase ändert
        # self.q_global = [scipy.sparse.csr_matrix(np.zeros(self.m))]  # speichere für jede globale Phase alle momentanen Warteschlangenlängen zu Beginn der Phase

        self.fp = [[] for _ in range(self.I)]  # f^+    [[],[]]
        self.fp_ind = [[] for _ in range(self.I)]  # Liste der Indizes der Phasen # NUR FÜR OutputTable BENÖTIGT UND KOSTET SPEICHER (?)
        self.fm = [[] for _ in range(self.I)]  # f^-
        self.q_global = [[0] for _ in range(self.m)]  # enthält für jede Kante alle 'piecewise_linear' Eckpunkte der Warteschlangen
        self.q_ind = [[0] for _ in range(self.m)]  # enthält die zu diesen Eckpunkten zugehörigen Phasen
        self.q = np.zeros(self.m)  # enthält alle aktuellen Warteschalngenlängen
        for i in range(self.I):
            for k in range(self.m):
                self.fp[i].append([(0, 0)])
                self.fp_ind[i].append([])
                self.fm[i].append([(0, 0)])

        for e_ind in range(self.m):  # Initialisierung "self.c" (Kosten)
            self.c[e_ind] = self.r[e_ind]

        # Zeitpunkte, zu denen sich der Zufluss in mindestens einem Quellknoten ändert
        self.u_start = set()
        for com in self.u:
            for s_list in com:
                for t in s_list:
                    self.u_start.add(t[0])

        self.graphReversed = self.reverse_graph(G)

        self.E_active = np.ones((self.I, self.m))
        for i in range(self.I):
            self.labels.append(self.dijkstra(self.graphReversed, 't{}'.format(i+1), visited=[], distances={}))

        self.E_active = scipy.sparse.lil_matrix((self.I, self.m))
        for (v_ind, v) in enumerate(self.V):
            outneighbors = self.G[v].keys()
            self.delta_p.append([])
            self.delta_p[v_ind] = outneighbors
            for w in outneighbors:
                w_ind = self.V.index(w)
                edge = self.E.index((v, w))
                for i in range(self.I):
                    if abs(self.labels[i][v_ind] - self.labels[i][w_ind] - self.c[edge]) < self.eps:
                        self.E_active[i, edge] = 1

        self.time_vor_main = time.time()
        self.main()

    # Quelle:
    # https://www.geeksforgeeks.org/topological-sorting/#:~:text=Topological%20sorting%20for%20Directed%20Acyclic,4%202%203%201%200%E2%80%9D
    # A recursive function used by topologicalSort
    def topologicalSortUtil(self, v_ind, i, visited, stack, dp_act_i):

        # Mark the current node as visited.
        visited[v_ind] = True

        # Recur for all the vertices adjacent to this vertex
        for e_ind in dp_act_i[v_ind]:
            w_ind = self.V.index(self.E[e_ind][1])
            if not visited[w_ind]:
                self.topologicalSortUtil(w_ind, i, visited, stack, dp_act_i)

        # Push current vertex to stack which stores result
        stack.append(v_ind)

    # The function to do Topological Sort. It uses recursive
    # topologicalSortUtil()
    def topologicalSort(self, i, dp_act_i):
        # Mark all the vertices as not visited
        ti_ind = self.V.index('t{}'.format(i+1))
        visited = [False]*self.n
        visited[ti_ind] = True
        stack = [ti_ind]

        # Call the recursive helper function to store Topological
        # Sort starting from all vertices one by one
        for k in range(self.n):
            if not visited[k]:
                self.topologicalSortUtil(k, i, visited, stack, dp_act_i)
        return stack

    # Quelle: http://www.gilles-bertrand.com/2014/03/dijkstra-algorithm-python-example-source-code-shortest-path.html
    def dijkstra(self, graph, src, visited=[], distances={}, predecessors={}):
        """
        Berechnet rekursiv kürzesten Weg und dessen Kosten von "src" zu jedem erreichbaren Knoten in "graph"
        :param graph: Graph als Dictionary
        :param src: Startknoten
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
        while True:
            # visit the neighbors
            for neighbor in graph[src]:
                if neighbor not in visited:
                    new_distance = distances[src] + self.c[self.E.index((neighbor,src))]
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
            src=min(unvisited, key=unvisited.get)

    # Quelle: http://www.gilles-bertrand.com/2014/03/dijkstra-algorithm-python-example-source-code-shortest-path.html
    def dijkstraalt(self, graph, src, visited=[], distances={}, predecessors={}):
        """
        Berechnet rekursiv kürzesten Weg und dessen Kosten von "src" zu jedem erreichbaren Knoten in "graph"
        :param graph: Graph als Dictionary
        :param src: Startknoten
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
                new_distance = distances[src] + self.c[self.E.index((neighbor,src))]
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
        return self.dijkstraalt(graph, x, phase_ind, visited, distances, predecessors)

    @staticmethod
    def reverse_graph(graph):
        """
        ersetzt alle Kanten in "graph" (gegeben als Dictionary) durch die jeweilige Rückwärtskante, Kosten bleiben
        gleich
        :param graph: Eingabegraph, dessen Kantenorientierungen vertauscht werden sollen
        :return: "graphReversed", entspricht "graph" mit vertauschten Kantenorientierungen als Dictionary
        """
        graphReversed = {}
        nodes = graph.keys()
        for v in nodes:
            graphReversed[v] = {}
        for v in nodes:
            for w in graph[v]:
                graphReversed[w][v] = graph[v][w]
        return graphReversed

    def waterfilling_algo(self, v, flow_vol, bi, outneighbors, wi_slope):
        """
        Berechnet Flussaufteilung des in 'v' vorhandenen Flussvolumens auf die ausgehenden aktiven Kanten
        :param v: Knoten
        :param flow_vol: gesamtes (alle Güter) momentanes Flussvolumen im Knoten 'v'
        :param bi: Flussmenge (eines bestimmten Guts) im Knoten 'v'
        :param outneighbors: Liste aller Knoten, die über eine direkte aktive (für das entsprechende Gut) Kante von 'v' aus erreichbar sind
        :param wi_slope: Liste über Änderung der Labels (für dieses Gut) aller Knoten in 'outneighbors'
        :return: Flussaufteilung 'z' auf alle von 'v' ausgehenden aktiven Kanten in Form einer Liste
        """
        pk = len(outneighbors)
        if pk == 0:
            return []
        if bi < self.eps:
            return np.zeros(pk)
        flow_ratio = float(flow_vol)/bi
        alpha = np.zeros(pk)
        beta = np.copy(wi_slope)
        gamma = np.zeros(pk)
        old_sorting = np.zeros(pk)

        # setze 'beta' - Werte
        for i in range(pk):
            e_ind = self.E.index((v, outneighbors[i]))
            if self.q[e_ind] > 0:
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
            alpha[i] += self.nu[e_ind] * flow_ratio
            if self.q[e_ind] == 0:
                gamma[i] += self.nu[e_ind] * flow_ratio

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
            if sum > bi:
                r -= 1
                break
        if r == pk + 1:
            r -= 1

        else_case = False
        if r < pk:
            sum = 0
            for i in range(r):
                sum += find_max(i, r)
            if sum <= bi:
                z_sum = 0
                for i in range(r):
                    z[i] += find_max(i, r)
                    z_sum += z[i]
                # z_sum -= z[r - 1]
                z[r] += bi - z_sum
            else:
                else_case = True
        else:
            else_case = True

        if else_case:
            z_sum = 0
            for i in range(r):
                z[i] += find_max(i, r - 1)
                z_sum += z[i]
            b_prime = bi - z_sum
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

    def calc_b(self, i, v_ind, phase):
        """
        Berechnet zum Zeitpunkt 'phase' im Knoten 'v' vorhandene Flussmenge b_{i,v}^- (phase) von Gut 'i'
        :param i: Gut
        :param v_ind: Index des Knotens
        :param phase: Zeitpunkt
        :return: b_{i,v}^- (phase)
        """
        b = 0
        if self.V[v_ind] == 't{}'.format(i+1):
            return 0
        in_paths = self.get_ingoing_edges(v_ind)
        for e_ind in in_paths:
            ind = self.last_fm_change(i, e_ind, phase)
            b += self.fm[i][e_ind][ind][1]
        u_v = self.u[i][v_ind]
        u_v_len = len(u_v)
        for tuple_ind in range(u_v_len - 1, -1, -1):
            if u_v[tuple_ind][0] <= phase:
                b += u_v[tuple_ind][1]
                break
        return b

    def get_ingoing_edges(self, v_ind):
        """
        bestimmt alle in 'v_ind' eingehenden Kanten
        :param v_ind: Knotenindex
        :return: Liste der Indizes der Kanten
        """
        v = self.V[v_ind]
        preds = self.graphReversed[v].keys()
        return [self.E.index((u,v)) for u in preds]

    def get_ingoing_active_edges(self, v, i):
        """
        bestimmt alle in 'v' eingehende, momentan für Gut 'i' aktive Kanten
        :param v: Knoten
        :param i: Gut
        :return: Liste der Indizes der aktiven Kanten
        """
        preds = self.graphReversed[v].keys()
        delta_m = [self.E.index((u,v)) for u in preds]
        return [e for e in delta_m if self.E_active[i, e]]

    def get_outgoing_edges(self, v_ind):
        """
        bestimmt alle aus 'v_ind' ausgehenden Kanten
        :param v_ind: Knoten
        :return: Liste der Indizes der Kanten
        """
        v = self.V[v_ind]
        gv_keys = self.G[v].keys()
        return [self.E.index((v,u)) for u in gv_keys]

    def get_outgoing_active_edges(self, i, v_ind, delta_p=None):
        """
        bestimmt alle aus 'v' ausgehende, momentan für Gut 'i' aktive Kanten
        :param i: Index des Guts
        :param v_ind: Knoten
        :param delta_p: Menge aller von 'v' ausgehenden Kanten. Optional: wird berechnet, falls nicht angegeben
        :return: Liste der Indizes der aktiven Kanten
        """
        if not delta_p:
            v = self.V[v_ind]
            outneighbors = self.G[v].keys()
            delta_p = [self.E.index((v,u)) for u in outneighbors]
        return [e for e in delta_p if self.E_active[i, e]]

    def change_of_cost(self, e_ind, z):
        """
        Gibt die momentane Kostenänderung der Kante mit Index 'e_ind' bei Zufluss 'z' an.
        :param e_ind: Index der betrachteten Kante
        :param z: Einflussrate in Kante 'e'
        :return: Änderungsrate der Kosten
        """
        dif = z - self.nu[e_ind]
        if abs(dif) < self.bd_tol * self.I:
            return 0
        if self.q[e_ind] > self.eps:
            return dif / self.nu[e_ind]
        return np.max([dif / self.nu[e_ind], 0])

    def last_fm_change(self, i, e_ind, theta):
        """
        Bestimmt den größten Zeitpunkt <= 'theta', zu dem sich der f_i^- -Wert von Kante 'e_ind' ändert.
        :param i: Index des Guts
        :param e_ind: Kantenindex
        :param theta: aktueller Zeitpunkt
        :return: Index des gesuchten Zeitpunkts
        """
        fm_len = len(self.fm[i][e_ind])
        for t in range(fm_len - 1, 0, -1):
            if self.fm[i][e_ind][t][0] < theta + 2 * self.bd_tol:
                return t
        return 0

    def g(self, e_ind, x):
        """
        Bestimmt absolute Änderung der Warteschlange von Kante 'e_ind' zum Zeitpunkt 'theta_ind' bei Einfluss 'x'.
        :param e_ind: Kantenindex
        :param x: momentaner Einfluss in Kante 'e_ind'
        :return: momentane Änderungsrate der Warteschlange
        """
        dif = x - self.nu[e_ind]
        if abs(dif) < self.I * self.bd_tol:
            return 0
        if self.q[e_ind] > self.eps:
            return dif
        return max([dif, 0])

    def gamma(self, x_sum, a, delta_p_act, coms):
        """
        Berechnet EIN Element aus der Menge Gamma(x) der Gamma Funktion aus paper. Im Fall dass diese Menge einen Fixpunkt besitzt (mit Konvention eindeutig), wird dieser auch
        ausgegeben.
        :param x:
        :param x_sum:
        :param a:
        :param delta_p_act:
        :param coms:
        :param theta_ind:
        :return:
        """
        # gx = scipy.sparse.lil_matrix((self.I, self.m))
        forced_zeros = []
        for i in range(self.I):
            forced_zeros.append({})
        for v_ind in coms:
            for i in coms[v_ind]:
                forced_zeros[i][v_ind] = []
                # x_sum_nz = 0
                # cor = 0
                for e_ind in delta_p_act[i][v_ind]:
                    w_ind = self.V.index(self.E[e_ind][1])
                    if abs(a[i, v_ind] - self.g(e_ind, x_sum[e_ind]) / self.nu[e_ind] - a[i, w_ind]) > self.eps:
                        # 'gx[i, e_ind]' = 0 muss gelten
                        forced_zeros[i][v_ind].append(e_ind)
                        """# merke Unterschied zwischen 'x[i, e_ind]' und 'gx[i, e_ind]' um später Flusserhaltung wiederherzustellen
                        cor += x[i, e_ind]
                    else:
                        # merke bereits festgelegte Flussmenge
                        x_sum_nz += x[i, e_ind]
                dp_nz = list(set(delta_p_act[i][v_ind]) - set(forced_zeros[i][v_ind]))
                if x_sum_nz > self.eps:
                    for e_ind in dp_nz:
                        # Verteile übriges Flussvolumen ('cor' viel) auf aktive, nicht-0-Kanten im Verhältnis der bereits zugeteilten Flusswerte
                        gx[i, e_ind] = x[i, e_ind] * (1 + cor / x_sum_nz)
                else:
                    dp_nz_len = len(dp_nz)
                    for e_ind in dp_nz:
                        gx[i, e_ind] = cor / dp_nz_len"""
        return forced_zeros

    def get_items(self, s):
        """
        Bestimmt Indizes der gesetzen Werte der (dünnbesetzten) Matrix 's'
        :param s: Matrix vom Typ 'scipy.sparse.lil_matrix'
        :return: Menge der 2-Tupel mit gesuchten Indizes
        """
        s_coo = s.tocoo()
        return set(zip(s_coo.row, s_coo.col))

    def fp_approx(self, theta_ind):
        theta = self.global_phase[theta_ind]
        xk = scipy.sparse.lil_matrix((self.I, self.m))
        xfp = scipy.sparse.lil_matrix((self.I, self.m))
        x_sum = np.zeros(self.m)
        x_sum_fix = np.zeros(self.m)
        delta_p_act = []
        delta_p_inact = []
        nu_sum_act = []
        coms = defaultdict(list)
        coms_lens = np.zeros(self.n)
        # self.b_ot.append(scipy.sparse.lil_matrix((self.I, self.n)))
        self.b = scipy.sparse.lil_matrix((self.I, self.n))
        # self.flow_vol.append(scipy.sparse.lil_matrix((self.n, 1)))  # NUR FÜR OutputTable BENÖTIGT UND KOSTET SPEICHER

        def calc_flow_by_bounds(xk_old=None):
            """
            :param xk_old:
            :return:
            """

            xk = scipy.sparse.lil_matrix((self.I, self.m))
            x_sum = np.zeros(self.m)

            for v_ind in coms:
                # enthält für jedes Gut in 'v_ind' Menge des bereits einer Kante zugewiesenen Flussvolumens
                flow_sent = {}
                for i in coms[v_ind]:
                    bds = bounds[v_ind][i].keys()
                    if not bds:
                        for e_ind in delta_p_act_nf[i][v_ind]:
                            xk[i, e_ind] = xk_old[i, e_ind]
                            x_sum[e_ind] += xk[i, e_ind]
                        continue
                    e_bds_nf = list(set(bds) - set(bounds_f[v_ind][i]))
                    flow_sent[i] = 0
                    for e_ind in e_bds_nf:
                        xk[i, e_ind] = 0.5 * (bounds[v_ind][i][e_ind][0] + bounds[v_ind][i][e_ind][1])
                        x_sum[e_ind] += xk[i, e_ind]
                        flow_sent[i] += xk[i, e_ind]
                    for e_ind in bounds_f[v_ind][i]:
                        xk[i, e_ind] = 0.5 * (bounds[v_ind][i][e_ind][0] + bounds[v_ind][i][e_ind][1])
                        x_sum[e_ind] += xk[i, e_ind]
                    # in 2 Fällen muss korrigiert werden:
                    # Fall 1: es wurde mehr Fluss von Gut 'i' verschickt, als vorhanden ist -> skaliere Flusswerte nach unten
                    # Fall 2: es wurde nicht das gesamte Flussvolumen von Gut 'i' verschickt und es gibt keine weiteren aktiven Kanten -> skaliere Flusswerte nach oben
                    if flow_sent[i] > b_nf[i, v_ind] + self.eps:
                        flow_diff = flow_sent[i] - b_nf[i, v_ind]
                        int_diff = 0
                        for e_ind in e_bds_nf:
                            int_diff += xk[i, e_ind] - bounds[v_ind][i][e_ind][0]
                        for e_ind in e_bds_nf:
                            flow_diff_part = flow_diff * (xk[i, e_ind] - bounds[v_ind][i][e_ind][0]) / int_diff
                            xk[i, e_ind] -= flow_diff_part
                            x_sum[e_ind] -= flow_diff_part

                    elif flow_sent[i] < b_nf[i, v_ind] - self.eps:
                        flow_diff = b_nf[i, v_ind] - flow_sent[i]
                        int_diff = 0
                        for e_ind in e_bds_nf:
                            int_diff += bounds[v_ind][i][e_ind][1] - xk[i, e_ind]
                        for e_ind in e_bds_nf:
                            flow_diff_part = flow_diff * (bounds[v_ind][i][e_ind][1] - xk[i, e_ind]) / int_diff
                            xk[i, e_ind] += flow_diff_part
                            x_sum[e_ind] += flow_diff_part
            return xk, x_sum

        def calc_a():
            a = scipy.sparse.lil_matrix((self.I, self.n))
            a_min2 = scipy.sparse.lil_matrix((self.I, self.n))
            sparse_items = []
            for i in range(self.I):
                t_ind = self.V.index('t{}'.format(i+1))
                sparse_items.append((i, t_ind))
            argmin = []
            argmin_nf = []
            for i in range(self.I):
                argmin.append({})
                argmin_nf.append({})
                for v_ind in top_ords[i]:
                    for e_ind in delta_p_act[i][v_ind]:
                        w_ind = self.V.index(self.E[e_ind][1])
                        a_e = self.g(e_ind, x_sum[e_ind] + x_sum_fix[e_ind]) / self.nu[e_ind] + a[i, w_ind]
                        if (i, v_ind) not in sparse_items:
                            if abs(a_e) > 1 * self.bd_tol * np.max([1, coms_lens[v_ind]]):
                                a[i, v_ind] = a_e
                            sparse_items.append((i, v_ind))
                            argmin[i][v_ind] = [e_ind]
                            if e_ind in delta_p_act_nf[i][v_ind]:
                                argmin_nf[i][v_ind] = [e_ind]
                            else:
                                # min wird in bereits fixierter Kante angenommen. 'argmin_nf' enthält jedoch nur nicht-fixierte, minimale Kanten
                                argmin_nf[i][v_ind] = []
                        elif a_e < a[i, v_ind] + 1 * self.bd_tol * np.max([1, coms_lens[v_ind]]):  # '* coms_lens[v_ind]' -> Weil: der Flusswert jedes Guts kann 'self.bd_tol'/2 "daneben" liegen
                            if a_e < a[i, v_ind]:
                                a_min2[i, v_ind] = a[i, v_ind]
                                if abs(a_e) > 1 * self.bd_tol * np.max([1, coms_lens[v_ind]]):
                                    a[i, v_ind] = a_e
                                else:
                                    a[i, v_ind] = 0
                                argmin[i][v_ind] = [e_ind]
                                if e_ind in delta_p_act_nf[i][v_ind]:
                                    argmin_nf[i][v_ind] = [e_ind]
                                else:
                                    # min wird in bereits fixierter Kante angenommen. 'argmin_nf' enthält jedoch nur nicht-fixierte, minimale Kanten
                                    argmin_nf[i][v_ind] = []
                            elif e_ind in delta_p_act_nf[i][v_ind]:
                                argmin_nf[i][v_ind].append(e_ind)
                                argmin[i][v_ind].append(e_ind)
                            else:
                                argmin[i][v_ind].append(e_ind)
                        elif a_min2[i, v_ind] == 0:
                            a_min2[i, v_ind] = a_e

            return a, a_min2, argmin, argmin_nf

        def relax_bounds():
            """
            Überprüft für jedes relevante (i,v) - Paar ob einer von zwei Fällen eintritt:
                1.: Es sind alle aktiven, nichtminmalen, von 'v' ausgehenden Kanten fixiert. Da Aufteilung noch nicht im Sinne eines IDE-Flusses und müssen die Flusswerte zwischen
                den einzelnen Kanten noch umverteilt werden, dies ist aber nicht möglich, da nur Kante im 'argmin' umverteilt werden können, was
                2.:
            :return:
            """
            x_sum_nm = np.zeros(self.m)
            for v_ind in coms:
                for i in coms[v_ind]:
                    arg = list(set(delta_p_act[i][v_ind]) - set(argmin[i][v_ind]))
                    for e_ind in arg:
                        x_sum_nm[e_ind] += xk[i, e_ind]
            for i in range(self.I):
                for v_ind in top_ords[i]:
                    if v_ind in coms and i in coms[v_ind]:
                        arg = list(set(delta_p_act[i][v_ind]) - set(argmin[i][v_ind]))
                        arg_nf = list(set(arg) - set(bounds_f[v_ind][i]))
                        if not arg_nf:
                            # In diesem Fall sind alle aktiven nichtminimalen Kanten bereits fixiert. Damit können keine Fortschritte mehr gemacht werden -> relaxiere bds
                            for e_ind in arg:
                                if xk[i, e_ind] < self.bd_tol * (np.max([1, coms_lens[v_ind]]) + self.I+1):
                                    continue
                                w_ind = self.V.index(self.E[e_ind][1])
                                if self.q[e_ind] < self.eps:
                                    # Falls q = 0 -> g_e nicht-negativ -> kleinster a-Wert kann nicht erreicht werden -> wähle kleinstmögliche untere Schranke: 0
                                    bounds[v_ind][i][e_ind] = (0, xk[i, e_ind])
                                else:
                                    # relaxiere untere Schranke, sodass der minimale a-Wert erreicht werden kann
                                    bounds[v_ind][i][e_ind] = (bounds[v_ind][i][e_ind][0] - (x_sum[e_ind] + x_sum_fix[e_ind] - (1 + a[i, v_ind] - a[i, w_ind]) * self.nu[e_ind]) * xk[i, e_ind] / x_sum_nm[e_ind], xk[i, e_ind])
                                bounds_f[v_ind][i].remove(e_ind)
                                delta_p_act_nf[i][v_ind].append(e_ind)
                                b_nf[i, v_ind] += xk[i, e_ind]
                        if not argmin_nf[i][v_ind]:
                            # Dagegen sind in diesem Fall alle minimalen Kanten fixiert. Damit können ebenfalls keine Forschritte gemacht werden -> relaxiere bds
                            for e_ind in argmin[i][v_ind]:
                                # if xk[i, e_ind] > self.b_ot[-1][i, v_ind] - self.bd_tol:
                                if xk[i, e_ind] > self.b[i, v_ind] - self.bd_tol:
                                    continue
                                w_ind = self.V.index(self.E[e_ind][1])
                                xe = x_sum[e_ind] + x_sum_fix[e_ind]
                                x2 = (a_min2[i, v_ind] - a[i, w_ind]) * self.nu[e_ind] - self.g(e_ind, xe)
                                if self.q[e_ind] < self.eps and xe < self.nu[e_ind] - self.eps:
                                    x2 += self.nu[e_ind] - xe
                                bounds_f[v_ind][i].remove(e_ind)
                                delta_p_act_nf[i][v_ind].append(e_ind)
                                # relaxierte obere Schranke mit x2: Dies ist der Flusswert, sodass der zweitkleinste a-Wert (= 'a_min2[i, v_ind]') auch für Kante 'e_ind'
                                # erreicht werden kann.
                                # bounds[v_ind][i][e_ind] = (xk[i, e_ind], np.min([x2 + xk[i, e_ind], self.b_ot[-1][i, v_ind]]))
                                bounds[v_ind][i][e_ind] = (xk[i, e_ind], np.min([x2 + xk[i, e_ind], self.b[i, v_ind]]))
                                b_nf[i, v_ind] += xk[i, e_ind]
            return

        def fix_nodes():
            coms_keys = list(coms)
            for v_ind in coms_keys:
                # Prüfe, ob alle aktiven Kanten den gleichen a-Wert liefern
                # ??? Eigentlich nur für Güter mit nichtdisjunkten Mengen delta_p_act_nf[i][v_ind] ?
                if np.all([v_ind not in argmin_nf[i] or len(argmin_nf[i][v_ind]) == len(delta_p_act_nf[i][v_ind]) for i in coms[v_ind]]):
                    # alle Kanten, für die ein Flusswert >0 fixiert wurde, müssen ebenfalls in 'argmin[i][v_ind]' liegen
                    if np.any([[e_ind for e_ind in bounds_f[v_ind][i] if xk[i, e_ind] > self.bd_tol and e_ind not in argmin[i][v_ind]] for i in coms[v_ind]]):
                        # sonst -> Abbruch
                        continue
                    outneighbors = list(set(np.concatenate([[self.V.index(self.E[e_ind][1]) for e_ind in argmin_nf[i][v_ind]] for i in coms[v_ind]])))
                    # sind zusätzlich für alle Endknoten 'w_ind' die Flussaufteilungen bereits fixiert, kann auch die Aufteilung in Knoten 'v_ind' fixiert werden
                    if np.all([w_ind not in coms for w_ind in outneighbors]):
                        # if v_ind == 0 and theta_ind == 0:
                        #    continue
                        for i in coms[v_ind]:
                            for e_ind in delta_p_act[i][v_ind]:
                                if xk[i, e_ind] < self.bd_tol * 2:  # (np.max([1, coms_lens[v_ind]]) + self.I+1):
                                    xfp[i, e_ind] = 0
                                # elif xk[i, e_ind] > self.b_ot[-1][i, v_ind] - self.bd_tol * 2:
                                elif xk[i, e_ind] > self.b[i, v_ind] - self.bd_tol * 2:
                                    # xfp[i, e_ind] = self.b_ot[-1][i, v_ind]
                                    xfp[i, e_ind] = self.b[i, v_ind]
                                    x_sum_fix[e_ind] += xfp[i, e_ind]
                                else:
                                    xfp[i, e_ind] = xk[i, e_ind]
                                    x_sum_fix[e_ind] += xfp[i, e_ind]
                                x_sum[e_ind] -= xk[i, e_ind]
                                xk[i, e_ind] = 0
                            delta_p_act_nf[i][v_ind] = []
                        del coms[v_ind]
            return

        # Initialisiere x0
        for i in range(self.I):
            delta_p_act.append({})
            delta_p_inact.append({})
            nu_sum_act.append({})
            for v_ind in range(self.n):
                delta_p = self.get_outgoing_edges(v_ind)
                delta_p_act[i][v_ind] = self.get_outgoing_active_edges(i, v_ind, delta_p=delta_p)
                delta_p_inact[i][v_ind] = [e for e in delta_p if e not in delta_p_act[i][v_ind]]
                # self.b_ot[-1][i, v_ind] = self.calc_b(i, v_ind, theta)
                self.b[i, v_ind] = self.calc_b(i, v_ind, theta)
                # if self.b_ot[-1][i, v_ind] > 0:
                if self.b[i, v_ind] > 0:
                    coms_lens[v_ind] += 1
                    # self.flow_vol[-1][v_ind, 0] += self.b_ot[-1][i, v_ind]
                    if len(delta_p_act[i][v_ind]) == 1:
                        e_ind = delta_p_act[i][v_ind][0]
                        # xfp[i, e_ind] = self.b_ot[-1][i, v_ind]
                        xfp[i, e_ind] = self.b[i, v_ind]
                        x_sum_fix[e_ind] += xfp[i, e_ind]
                        continue
                    coms[v_ind].append(i)
                    nu_sum_act[i][v_ind] = 0
                    for e_ind in delta_p_act[i][v_ind]:
                        nu_sum_act[i][v_ind] += self.nu[e_ind]
                    for e_ind in delta_p_act[i][v_ind]:
                        # xk[i, e_ind] = self.b_ot[-1][i, v_ind] * self.nu[e_ind] / nu_sum_act[i][v_ind]
                        xk[i, e_ind] = self.b[i, v_ind] * self.nu[e_ind] / nu_sum_act[i][v_ind]
                        x_sum[e_ind] += xk[i, e_ind]

        # Teilmenge von 'delta_p_act', welche für jedes Gut nur die Kanten mit nicht bereits fixiertem Flusswert enthält
        delta_p_act_nf = copy.deepcopy(delta_p_act)

        # Bestimme a_i,v - Werte zu x0
        a = scipy.sparse.lil_matrix((self.I, self.n))
        top_ords = []
        sparse_items = []
        for i in range(self.I):
            t_ind = self.V.index('t{}'.format(i+1))
            sparse_items.append((i, t_ind))
        argmin = []
        for i in range(self.I):
            argmin.append({})
            top_ords.append(self.topologicalSort(i, delta_p_act[i]))
            for v_ind in top_ords[i]:
                for e_ind in delta_p_act[i][v_ind]:
                    w_ind = self.V.index(self.E[e_ind][1])
                    a_e = self.change_of_cost(e_ind, x_sum[e_ind] + x_sum_fix[e_ind]) + a[i, w_ind]
                    if (i, v_ind) not in sparse_items:
                        if abs(a_e) > 1 * self.bd_tol * np.max([1, coms_lens[v_ind]]):
                            a[i, v_ind] = a_e
                        sparse_items.append((i, v_ind))
                        argmin[i][v_ind] = [e_ind]
                    elif a_e < a[i, v_ind] + 1 * self.bd_tol * np.max([1, coms_lens[v_ind]]):  # '* coms_lens[v_ind]' -> Weil: der Flusswert jedes Guts kann 'self.bd_tol'/2 "daneben" liegen
                        if a_e < a[i, v_ind]:
                            if abs(a_e) > 1 * self.bd_tol * np.max([1, coms_lens[v_ind]]):
                                a[i, v_ind] = a_e
                            else:
                                a[i, v_ind] = 0
                            argmin[i][v_ind] = [e_ind]
                        else:
                            argmin[i][v_ind].append(e_ind)

        coms_keys = list(coms)
        for v_ind in coms_keys:
            # für x0 ist 'argmin' = 'argmin_nf' und 'delta_p_act' = 'delta_p_act_nf'
            if np.all([len(argmin[i][v_ind]) == len(delta_p_act[i][v_ind]) for i in coms[v_ind]]):
                outneighbors = list(set(np.concatenate([[self.V.index(self.E[e_ind][1]) for e_ind in argmin[i][v_ind]] for i in coms[v_ind]])))
                # sind zusätzlich für alle Endknoten 'w_ind' die Flussaufteilungen bereits fixiert, kann auch die Aufteilung in Knoten 'v_ind' fixiert werden
                if np.all([w_ind not in coms for w_ind in outneighbors]):
                    for i in coms[v_ind]:
                        for e_ind in delta_p_act[i][v_ind]:
                            # wirklich so?
                            if xk[i, e_ind] < self.bd_tol * (np.max([1, coms_lens[v_ind]]) + self.I+1):
                                xfp[i, e_ind] = 0
                            else:
                                xfp[i, e_ind] = xk[i, e_ind]
                                x_sum_fix[e_ind] += xfp[i, e_ind]
                        delta_p_act_nf[i][v_ind] = []
                    del coms[v_ind]
                    if not coms:
                        # Abbruchbedingung: keine Güter mehr übrig, die noch auf Kanten aufgeteilt werden müssen
                        return xfp, x_sum_fix, a, top_ords, coms_lens, delta_p_act, delta_p_inact

        bounds = {}  # enthält untere und obere Schranken für 'xfp[i, e_ind]' - Werte
        bounds_f = {}
        for v_ind in coms:
            bounds[v_ind] = {}
            bounds_f[v_ind] = {}
            for i in coms[v_ind]:
                bounds[v_ind][i] = {}
                bounds_f[v_ind][i] = []
        # b_nf = copy.deepcopy(self.b_ot[-1])
        b_nf = copy.deepcopy(self.b)
        fp_comp = []

        while True:
            # print("XK", xk)
            # print(theta)
            # nur im ersten Schritt leer, danach enthält 'fp_comp' Tupel der Form '(i, v_ind)', was impliziert, dass die Aufteilung von Gut 'i' im Knoten 'v_ind'
            # vollständig fixiert wurde (d.h. die Boundsintervalle aller aktiven ausgehenden Kanten sind hinreichend klein)
            for (i, v_ind) in fp_comp:
                for e_ind in bounds[v_ind][i]:
                    if bounds[v_ind][i][e_ind][1] - bounds[v_ind][i][e_ind][0] > self.bd_tol:
                        # In diesem Fall ist geforderte Toleranz 'bd_diff' nicht erreicht, da '(i, v_ind)' aber trotzdem in 'fp_comp', ist bereits 'b_nf[i, v_ind]' 0, d.h.
                        # eine solche Kante erhält Flusswert 0
                        # ??? tritt dieser Fall überhaupt auf, oder wird er schon in 'fix_nodes()' behandelt ???
                        xfp[i, e_ind] = 0
                        continue
                    xfp[i, e_ind] = (bounds[v_ind][i][e_ind][0] + bounds[v_ind][i][e_ind][1]) * 0.5
                    if xfp[i, e_ind] < self.bd_tol * (np.max([1, coms_lens[v_ind]]) + self.I+1):
                        xfp[i, e_ind] = 0
                    # elif xfp[i, e_ind] > self.b_ot[-1][i, v_ind] - 2 * self.bd_tol:
                    elif xfp[i, e_ind] > self.b[i, v_ind] - 2 * self.bd_tol:
                        # xfp[i, e_ind] = self.b_ot[-1][i, v_ind]
                        xfp[i, e_ind] = self.b[i, v_ind]
                        x_sum_fix[e_ind] += xfp[i, e_ind]
                    else:
                        x_sum_fix[e_ind] += xfp[i, e_ind]
                del bounds[v_ind][i]
                del bounds_f[v_ind][i]
                coms[v_ind].remove(i)
                if not coms[v_ind]:
                    del coms[v_ind]

            if not coms:  # Abbruchbedingung: keine Güter mehr übrig, die noch auf Kanten aufgeteilt werden müssen
                return xfp, x_sum_fix, a, top_ords, coms_lens, delta_p_act, delta_p_inact

            if fp_comp:  # nur im 1. Schritt nicht erfüllt. In diesem wird x0 verwendet
                xk, x_sum = calc_flow_by_bounds()

            a, a_min2, argmin, argmin_nf = calc_a()

            relax_bounds()

            fix_nodes()

            if not coms:  # Abbruchbedingung: keine Güter mehr übrig, die noch auf Kanten aufgeteilt werden müssen
                return xfp, x_sum_fix, a, top_ords, coms_lens, delta_p_act, delta_p_inact

            while True:
                print(theta_ind, theta)
                forced_zeros = self.gamma(x_sum + x_sum_fix, a, delta_p_act_nf, coms)
                no_fz = True
                for v_ind in coms:
                    for i in coms[v_ind]:
                        if forced_zeros[i][v_ind]:
                            no_fz = False
                            break
                    if not no_fz:
                        break
                # all_fz = [[forced_zeros[i][v_ind] for i in coms[v_ind]] for v_ind in coms]
                if no_fz:
                    sol = xfp + xk
                    sparse_items = self.get_items(xk)
                    for (i, e_ind) in sparse_items:
                        if xk[i, e_ind] < self.bd_tol * (np.max([1, coms_lens[self.V.index(self.E[e_ind][0])]]) + self.I+1):
                            x_sum[e_ind] -= xk[i, e_ind]
                            sol[i, e_ind] -= xk[i, e_ind]
                    return sol, x_sum + x_sum_fix, a, top_ords, coms_lens, delta_p_act, delta_p_inact

                """ide_err = {}
                ide_delta_err = {}
                for v_ind in coms:
                    ide_err[v_ind] = {}
                    ide_delta_err[v_ind] = {}
                    for i in coms[v_ind]:
                        max_cost = -np.Inf
                        min_cost = np.Inf
                        maxes = []
                        mins = []
                        for e in delta_p_act[i][v_ind]:
                            w = self.V.index(self.E[e][1])
                            if self.labels[i][w] + self.c[e] < min_cost:
                                min_cost = self.labels[i][w] + self.c[e]
                                mins = [e]
                            elif self.labels[i][w] + self.c[e] == min_cost:
                                mins.append(e)
                            if self.labels[i][w] + self.c[e] > max_cost:
                                max_cost = self.labels[i][w] + self.c[e]
                                maxes = [e]
                            elif self.labels[i][w] + self.c[e] == max_cost:
                                maxes.append(e)
                        ide_err[v_ind][i] = max_cost - min_cost
                        ide_delta_err[v_ind][i] = [[], []]
                        for e in mins:
                            w = self.V.index(self.E[e][1])
                            ide_delta_err[v_ind][i][0].append(a[i, w] + self.g(e, x_sum[e] + x_sum_fix[e]))
                        for e in maxes:
                            w = self.V.index(self.E[e][1])
                            ide_delta_err[v_ind][i][1].append(a[i, w] + self.g(e, x_sum[e] + x_sum_fix[e]))"""

                fp_comp = []
                for v_ind in coms:
                    coms_len = len(coms[v_ind])  # beachtet Löschen der bereits fixierten Güter (anders als 'coms_lens[v_ind]')
                    for i in coms[v_ind]:
                        fz_nf = list(set(forced_zeros[i][v_ind]) - set(bounds_f[v_ind][i]))
                        for e_ind in fz_nf:
                            if e_ind in bounds[v_ind][i]:
                                bounds[v_ind][i][e_ind] = (bounds[v_ind][i][e_ind][0], xk[i, e_ind])
                            else:
                                bounds[v_ind][i][e_ind] = (0, xk[i, e_ind])
                            bd_diff = abs(bounds[v_ind][i][e_ind][0] - bounds[v_ind][i][e_ind][1])
                            if bd_diff < self.bd_tol and e_ind not in bounds_f[v_ind][i]:
                                if bd_diff == 0 and xk[i, e_ind] > self.eps:
                                    # bd wird relaxiert
                                    bounds[v_ind][i][e_ind] = (0.5 * xk[i, e_ind], xk[i, e_ind])
                                else:
                                    bounds_f[v_ind][i].append(e_ind)
                                    b_nf[i, v_ind] -= xk[i, e_ind]
                                    delta_p_act_nf[i][v_ind].remove(e_ind)
                                    if b_nf[i, v_ind] < self.bd_tol * coms_len:
                                        fp_comp.append((i, v_ind))
                        if v_ind not in argmin_nf[i] or len(argmin_nf[i][v_ind]) == len(delta_p_act_nf[i][v_ind]):
                            # Anpassung der bounds nicht notwendig
                            continue
                        for e_ind in argmin_nf[i][v_ind]:
                            if e_ind in bounds[v_ind][i]:
                                bounds[v_ind][i][e_ind] = (xk[i, e_ind], bounds[v_ind][i][e_ind][1])
                            else:
                                # bounds[v_ind][i][e_ind] = (xk[i, e_ind], self.b_ot[-1][i, v_ind])
                                bounds[v_ind][i][e_ind] = (xk[i, e_ind], self.b[i, v_ind])
                            bd_diff = abs(bounds[v_ind][i][e_ind][0] - bounds[v_ind][i][e_ind][1])
                            if bd_diff < self.bd_tol and e_ind not in bounds_f[v_ind][i]:
                                # if bd_diff == 0 and xk[i, e_ind] < self.b_ot[-1][i, v_ind] - self.eps:
                                if bd_diff == 0 and xk[i, e_ind] < self.b[i, v_ind] - self.eps:
                                    # bd wird relaxiert
                                    # bounds[v_ind][i][e_ind] = (xk[i, e_ind], np.min([2 * xk[i, e_ind], self.b_ot[-1][i, v_ind]]))
                                    bounds[v_ind][i][e_ind] = (xk[i, e_ind], np.min([2 * xk[i, e_ind], self.b[i, v_ind]]))
                                else:
                                    bounds_f[v_ind][i].append(e_ind)
                                    b_nf[i, v_ind] -= xk[i, e_ind]
                                    delta_p_act_nf[i][v_ind].remove(e_ind)
                                    if b_nf[i, v_ind] < self.bd_tol * coms_len:
                                        fp_comp.append((i, v_ind))
                                        break

                if fp_comp:
                    # es existieren Knoten in 'fp_comp', welche fixiert werden können
                    break

                xk, x_sum = calc_flow_by_bounds(xk_old=xk)

                a, a_min2, argmin, argmin_nf = calc_a()

                relax_bounds()

                fix_nodes()

                if not coms:
                    # Abbruchbedingung: keine Güter mehr übrig, die noch auf Kanten aufgeteilt werden müssen
                    return xfp, x_sum_fix, a, top_ords, coms_lens, delta_p_act, delta_p_inact

    def main(self):
        """
        Hauptteil: Konstruiert schrittweise einen kontinuierlichen IDE-Fluss zu den Eingabedaten aus 'Math.cont_data'.
        :return: 0
        """
        theta = 0
        theta_ind = -1
        # Obergrenze für theta
        T = 14
        # s_err = [[], [], []]
        while theta < T:
            # in der Zukunft liegende Zeitpunkte aus der Liste 'self.u_start'
            start_points = [t for t in self.u_start if t > theta]
            theta_ind += 1
            # Liste aller Kanten, deren Warteschlange in der aktuellen Phase 0 wird, und deren f^- -Werte auf 0
            # gesetzt werden müssen
            fm_to_zero = {}
            # 'next_phase' bestimmt am Ende der Schleife, wie lange aktuelle Phase andauert
            if start_points:
                next_phase = [np.min(start_points) - theta]
            else:
                next_phase = [T]
            x_total, x_sum, a, top_ords, coms_lens, delta_p_act, delta_p_inact = self.fp_approx(theta_ind)

            """for i in range(self.I):
                max_err = 0
                for v_ind in range(self.n):
                    delta_p_v_act = self.get_outgoing_active_edges(i, v_ind)
                    vmax = -np.Inf
                    vmin = np.Inf
                    for e_ind in delta_p_v_act:
                        w_ind = self.V.index(self.E[e_ind][1])
                        val = a[i, w_ind] + self.g(e_ind, x_sum[e_ind]) / self.nu[e_ind]
                        if x_total[i, e_ind] > 0:
                            # if x_total[i, e_ind] < self.b_ot[-1][i, v_ind] - self.bd_tol * 2:
                            if x_total[i, e_ind] < self.b[i, v_ind] - self.bd_tol * 2:
                                if val > vmax:
                                    vmax = val
                                if val < vmin:
                                    vmin = val
                            else:
                                # fc_diff = self.b_ot[-1][i, v_ind] - x_total[i, e_ind]
                                fc_diff = self.b[i, v_ind] - x_total[i, e_ind]
                                if fc_diff > max_err:
                                    max_err = fc_diff
                    if vmax - vmin > max_err:
                        max_err = vmax - vmin
                if not s_err[i] or max_err != s_err[i][-1][1]:
                    s_err[i].append((theta, max_err))"""

            for ti in range(self.I):
                for v_ind in top_ords[ti]:
                    # betrachte aktive Kanten
                    for e_ind in delta_p_act[ti][v_ind]:
                        if theta == 0:
                            if x_total[ti, e_ind] > self.bd_tol:
                                self.fp[ti][e_ind][0] = (0, x_total[ti, e_ind])
                                self.fp_ind[ti][e_ind].append(0)
                        # falls sich f^+ -Wert in dieser Phase ändert, aktualisiere 'self.fp'
                        elif abs(self.fp[ti][e_ind][-1][1] - x_total[ti, e_ind]) > self.bd_tol:
                            self.fp[ti][e_ind].append((theta, x_total[ti, e_ind]))
                            self.fp_ind[ti][e_ind].append(theta_ind)
                            if self.q_ind[e_ind][-1] != theta_ind and (self.q[e_ind] > self.eps or self.q_global[e_ind][-1] > self.eps or x_sum[e_ind] > self.nu[e_ind]):
                                self.q_ind[e_ind].append(theta_ind)
                        if x_sum[e_ind] > self.bd_tol * coms_lens[v_ind]:  # coms_lens ???
                            flow_ratio = x_total[ti, e_ind] / x_sum[e_ind]
                        else:
                            flow_ratio = 0
                        if self.q[e_ind] > self.eps:
                            outflow_time = theta + self.q[e_ind]/self.nu[e_ind] + self.r[e_ind]
                            outflow = self.nu[e_ind] * flow_ratio
                        else:
                            outflow_time = theta + self.r[e_ind]
                            if x_sum[e_ind] <= self.nu[e_ind]:
                                outflow = x_total[ti, e_ind]
                            else:
                                outflow = self.nu[e_ind] * flow_ratio
                        fm_ind = self.last_fm_change(ti, e_ind, outflow_time)
                        # falls sich f_{ti}^- -Wert durch die Flussaufteilung in dieser Phase ändert, aktualisiere 'self.fm'
                        if abs(self.fm[ti][e_ind][fm_ind][1] - outflow) > self.bd_tol:
                            if abs(self.fm[ti][e_ind][fm_ind][0] - outflow_time) < self.eps:
                                del self.fm[ti][e_ind][fm_ind]
                            else:
                                self.fm[ti][e_ind].append((outflow_time, outflow))
            for ti in range(self.I):
                # überspringen jeweilige Senke, da diese hier uninteressant
                for v_ind in top_ords[ti][1:]:
                    for e_ind in delta_p_act[ti][v_ind]:
                        last_fm_ind = self.last_fm_change(ti, e_ind, theta)
                        # bestimme, ob sich vor dem Ende der aktuellen Phase ein f^- -Wert ändert -> verkürze Phase
                        if len(self.fm[ti][e_ind]) > last_fm_ind + 1 and self.eps < self.fm[ti][e_ind][last_fm_ind + 1][0] - theta < next_phase[0] + 2 * self.bd_tol * self.I:
                            next_fm = self.fm[ti][e_ind][last_fm_ind + 1][0] - theta
                            if next_fm < next_phase[0]:
                                nph_rev = copy.deepcopy(next_phase)
                                nph_rev.reverse()
                                for p in nph_rev:
                                    if p > next_fm + 2 * self.bd_tol * self.I:
                                        next_phase.pop()
                                        if p in fm_to_zero:
                                            # wird zurückgesetzt, da durch Verkürzung der Phase der f^- -Wert erst in einer späteren Phase neu gesetzt wird
                                            del fm_to_zero[p]
                                    else:
                                        break
                                next_phase.insert(0, next_fm)
                            else:
                                last_ind = len(next_phase) - 1
                                for p in range(last_ind, -1, -1):
                                    if next_fm >= next_phase[p]:
                                        if next_fm > next_phase[p]:
                                            next_phase.insert(p+1, next_fm)
                                        break

                    for e_ind in delta_p_inact[ti][v_ind]:
                        # falls f^+ -Wert vorher > 0 war, so wird dieser hier auf 0 gesetzt, da Kante inaktiv
                        if abs(self.fp[ti][e_ind][-1][1]) > self.eps:
                            self.fp[ti][e_ind].append((theta, 0))
                            self.fp_ind[ti][e_ind].append(theta_ind)
                            if self.q_ind[e_ind][-1] != theta_ind and (self.q[e_ind] > self.eps or self.q_global[e_ind][-1] > self.eps or x_sum[e_ind] > self.nu[e_ind]):
                                self.q_ind[e_ind].append(theta_ind)
                            if self.q[e_ind] > self.eps:
                                outflow_time = theta + self.q[e_ind]/self.nu[e_ind] + self.r[e_ind]
                            else:
                                outflow_time = theta + self.r[e_ind]
                            fm_ind = self.last_fm_change(ti, e_ind, outflow_time)
                            # setze 'self.fm'- Wert auf 0, falls dies noch nicht geschehen ist
                            if abs(self.fm[ti][e_ind][fm_ind][1]) > self.eps:
                                self.fm[ti][e_ind].append((outflow_time, 0))

                        # bestimme, ob sich vor dem Ende der aktuellen Phase ein f^- -Wert ändert -> verkürze Phase
                        last_fm_ind = self.last_fm_change(ti, e_ind, theta)
                        if len(self.fm[ti][e_ind]) > last_fm_ind + 1 and self.eps < self.fm[ti][e_ind][last_fm_ind + 1][0] - theta < next_phase[0] + 2 * self.bd_tol * self.I:
                            next_fm = self.fm[ti][e_ind][last_fm_ind + 1][0] - theta
                            if next_fm < next_phase[0]:
                                nph_rev = copy.deepcopy(next_phase)
                                nph_rev.reverse()
                                for p in nph_rev:
                                    if p > next_fm + 2 * self.bd_tol * self.I:
                                        next_phase.pop()
                                        if p in fm_to_zero:
                                            # wird zurückgesetzt, da durch Verkürzung der Phase der f^- -Wert erst in einer späteren Phase neu gesetzt wird
                                            del fm_to_zero[p]
                                    else:
                                        break
                                next_phase.insert(0, next_fm)
                            else:
                                last_ind = len(next_phase) - 1
                                for p in range(last_ind, -1, -1):
                                    if next_fm >= next_phase[p]:
                                        if next_fm > next_phase[p]:
                                            next_phase.insert(p+1, next_fm)
                                        break

                    len_act = len(delta_p_act[ti][v_ind])
                    for ei in range(len_act):
                        active_ind = delta_p_act[ti][v_ind][ei]
                        # bestimme Kante, die während der gesamten Phase für Gut 'ti' aktiv bleibt
                        if x_total[ti, active_ind] > 0 or ei == len_act - 1:
                            # Änderung der Kosten dieser Kante
                            active_change = self.change_of_cost(active_ind, x_sum[active_ind])
                            break
                    if len_act:
                        for e_ind in delta_p_inact[ti][v_ind]:
                            change = self.change_of_cost(e_ind, x_sum[e_ind])
                            # prüfe, wann inaktive Kanten unter momentanem Einfluss aktiv werden
                            tar_ind = self.V.index(self.E[e_ind][1])
                            act_ind = self.V.index(self.E[active_ind][1])
                            if self.labels[ti][tar_ind] + self.q[e_ind]/self.nu[e_ind] + self.r[e_ind] + \
                                    (change + a[ti, tar_ind]) * (next_phase[0] + 2 * self.bd_tol * self.I) < \
                                    self.labels[ti][act_ind] + self.q[active_ind]/self.nu[active_ind] + \
                                    self.r[active_ind] + (active_change + a[ti, act_ind]) * (next_phase[0] + 2 * self.bd_tol * self.I) and \
                                    abs(change + a[ti, tar_ind] - active_change - a[ti, act_ind]) \
                                    > self.eps:
                                time_ub = (self.labels[ti][act_ind] + self.q[active_ind]/self.nu[active_ind] + self.r[active_ind] - \
                                          self.labels[ti][tar_ind] - self.q[e_ind]/self.nu[e_ind] - self.r[e_ind]) \
                                          / (change + a[ti, tar_ind] - active_change - a[ti, act_ind])

                                if time_ub < next_phase[0]:
                                    nph_rev = copy.deepcopy(next_phase)
                                    nph_rev.reverse()
                                    for p in nph_rev:
                                        if p > time_ub + 2 * self.bd_tol * self.I:
                                            next_phase.pop()
                                            if p in fm_to_zero:
                                                # wird zurückgesetzt, da durch Verkürzung der Phase der f^- -Wert erst in einer späteren Phase neu gesetzt wird
                                                del fm_to_zero[p]
                                        else:
                                            break
                                    next_phase.insert(0, time_ub)
                                else:
                                    last_ind = len(next_phase) - 1
                                    for p in range(last_ind, -1, -1):
                                        if time_ub >= next_phase[p]:
                                            if time_ub > next_phase[p]:
                                                next_phase.insert(p+1, time_ub)
                                            break

            for e_ind in range(self.m):
                change_of_q = self.change_of_cost(e_ind, x_sum[e_ind]) * self.nu[e_ind]
                # Falls die Warteschlange von 'e' unter aktuellem Fluss abgebaut wird, bestimme Zeitpunkt, zu dem diese vollständig abgebaut ist (bei gleich bleibendem
                # Fluss)
                if change_of_q < -self.eps:
                    # 'phase_length': Dauer bis Warteschlangenlänge gleich 0
                    phase_length = - self.q[e_ind] / change_of_q
                    if phase_length < next_phase[0] + 2 * self.bd_tol * self.I:
                        if phase_length < next_phase[0]:
                            nph_rev = copy.deepcopy(next_phase)
                            nph_rev.reverse()
                            for p in nph_rev:
                                if p > phase_length + 2 * self.bd_tol * self.I:
                                    next_phase.pop()
                                    if p in fm_to_zero:
                                        # wird zurückgesetzt, da durch Verkürzung der Phase der f^- -Wert erst in einer späteren Phase neu gesetzt wird
                                        del fm_to_zero[p]
                                else:
                                    break
                                next_phase.insert(0, phase_length)
                            # prüfe ob Zufluss 0 und somit 'fm' auf 0 gesetzt werden muss
                            if change_of_q + self.nu[e_ind] < self.eps:
                                fm_to_zero[phase_length] = [e_ind]
                        else:
                            last_ind = len(next_phase) - 1
                            for p in range(last_ind, -1, -1):
                                if phase_length >= next_phase[p]:
                                    if phase_length > next_phase[p]:
                                        next_phase.insert(p+1, phase_length)
                                    break
                            if change_of_q + self.nu[e_ind] < self.eps:
                                if phase_length in fm_to_zero:
                                    fm_to_zero[phase_length].append(e_ind)
                                else:
                                    fm_to_zero[phase_length] = [e_ind]

            if theta_ind > 0:
                for e_ind in range(self.m):
                    if np.any([theta_ind in self.fp_ind[i][e_ind] for i in range(self.I)]) and (self.q[e_ind] > self.eps or self.q_global[e_ind][-1] > self.eps or x_sum[e_ind] > self.nu[e_ind]):
                        self.q_global[e_ind].append(self.q[e_ind])

            theta += next_phase[-1]
            if next_phase[-1] != T:
                # aktualisiere Warteschlangenlängen und Kosten
                new_q = []
                for e_ind in range(self.m):
                    if theta_ind == 2 and e_ind == 1:
                        print(5)
                    next_q_len = self.q[e_ind] + self.g(e_ind, x_sum[e_ind]) * next_phase[-1]
                    # if next_q_len < next_phase[-1] * self.bd_tol * self.I:
                    if next_q_len < 2 * self.bd_tol * self.I:
                        next_q_len = 0
                        if self.q[e_ind] > self.eps:
                            self.q_global[e_ind].append(0)
                            self.q_ind[e_ind].append(theta_ind + 1)
                    new_q.append(next_q_len)
                    self.c[e_ind] = new_q[e_ind] / self.nu[e_ind] + self.r[e_ind]
                # speichere aktuelle Warteschlangenlängen
                self.q = new_q

                # speichere Phase
                self.global_phase.append(theta)
                self.E_active = scipy.sparse.lil_matrix((self.I, self.m))
                qi = []
                for e_ind in range(self.m):
                    qi.append({})
                new_labels = []
                for i in range(self.I):
                    new_labels.append([])
                    for v_ind in range(self.n):
                        new_labels[i].append(self.labels[i][v_ind] + next_phase[-1] * a[i, v_ind])
                    # ??? 'top_ords' an dieser Stelle veraltet, wieso wird das hier verwendet?
                    for v_ind in top_ords[i][1:]:
                        if np.isinf(new_labels[i][v_ind]):
                            continue
                        v = self.V[v_ind]
                        outneighbors = self.G[v].keys()
                        dv_act = []
                        for w in outneighbors:
                            w_ind = self.V.index(w)
                            edge = self.E.index((v, w))
                            label_dif = new_labels[i][v_ind] - new_labels[i][w_ind] - self.c[edge]
                            if label_dif > 0 or abs(label_dif) < self.bd_tol * 4 * self.I * theta:
                                # ???
                                # or abs(label_dif) < self.bd_tol * (np.max([1, coms_lens[v_ind]]) + self.I+1):
                                # 'label_dif' > 0: geschieht nur aufgrund von Approximationsfehler, für Kanten die aktiv sein sollten
                                self.E_active[i, edge] = 1
                                dv_act.append(edge)
                        if len(dv_act) == 1:
                            # labels werden angepasst, damit diese durch Approximationsfehler nicht zu stark auseinanderdriften.
                            e = self.E[dv_act[0]]
                            w_ind = self.V.index(e[1])
                            self.labels[i][v_ind] = new_labels[i][w_ind] + self.c[dv_act[0]]
                            qi[dv_act[0]][i] = self.q[dv_act[0]]
                        else:
                            l_sum = 0
                            for e_ind in dv_act:
                                w_ind = self.V.index(self.E[e_ind][1])
                                l_sum += new_labels[i][w_ind] + self.c[e_ind]
                            self.labels[i][v_ind] = l_sum / len(dv_act)
                            for e_ind in dv_act:
                                w_ind = self.V.index(self.E[e_ind][1])
                                # 'qi[e_ind][i]': Warteschlangenlönge von Kante 'e_ind', passend zu den labels von Gut 'i'
                                qi[e_ind][i] = (self.labels[i][v_ind] - new_labels[i][w_ind] - self.r[e_ind]) * self.nu[e_ind]
                                # if qi[e_ind][i] < self.bd_tol * self.I * next_phase[-1]:
                                if qi[e_ind][i] < 2 * self.I * self.bd_tol:
                                    qi[e_ind][i] = 0

                for e_ind in range(self.m):
                    len_qie = len(qi[e_ind])
                    if len_qie > 1:
                        mean_qi = 0
                        for i in qi[e_ind]:
                            mean_qi += qi[e_ind][i]
                        mean_qi /= len_qie
                        # if mean_qi > next_phase[-1] * self.bd_tol * self.I:
                        if mean_qi > 2 * self.bd_tol * self.I:
                            self.q[e_ind] = mean_qi
                            self.c[e_ind] = self.q[e_ind] / self.nu[e_ind] + self.r[e_ind]

                for p in fm_to_zero:
                    for e_ind in fm_to_zero[p]:
                        for i in range(self.I):
                            if self.fm[i][e_ind][-1][1] > self.eps:
                                self.fm[i][e_ind].append((theta + self.r[e_ind], 0))

        """for i in range(self.I):
            # am Ende sind alle f^+ -Werte 0
            for e in range(self.m):
                if abs(self.fp[i][e][-1][1]) > self.eps:
                    self.fp[i][e].append((theta - next_phase[-1], 0))
                    self.fp_ind[i][e].append(theta_ind)"""

        """s_times0 = [t for (t, v) in s_err[0]]
        s_end_times0 = [t for (t, v) in s_err[0][1:]] + [self.global_phase[-1]]
        nz_times0 = []
        nz_vals0 = []
        for tvi in range(len(s_err[0])):
            if s_err[0][tvi][1] > self.eps:
                nz_times0.append(s_err[0][tvi][0])
                nz_times0.append(s_err[0][tvi+1][0])
                nz_vals0.append(s_err[0][tvi][1])
                nz_vals0.append(s_err[0][tvi][1])
        s_vals0 = [v for (t, v) in s_err[0]]
        plt.hlines(s_vals0, s_times0, s_end_times0, colors='r')
        plt.vlines(nz_times0, np.zeros(len(nz_vals0)), nz_vals0, colors='r', linestyles='dotted')
        s_times1 = [t for (t, v) in s_err[1]]
        s_end_times1 = [t for (t, v) in s_err[1][1:]] + [self.global_phase[-1]]
        nz_times1 = []
        nz_vals1 = []
        for tvi in range(len(s_err[1])):
            if s_err[1][tvi][1] > self.eps:
                nz_times1.append(s_err[1][tvi][0])
                nz_times1.append(s_err[1][tvi+1][0])
                nz_vals1.append(s_err[1][tvi][1])
                nz_vals1.append(s_err[1][tvi][1])
        s_vals1 = [v for (t, v) in s_err[1]]
        plt.hlines(s_vals1, s_times1, s_end_times1, colors='b')
        plt.vlines(nz_times1, np.zeros(len(nz_vals1)), nz_vals1, colors='b', linestyles='dotted')
        s_times2 = [t for (t, v) in s_err[2]]
        s_end_times2 = [t for (t, v) in s_err[2][1:]] + [self.global_phase[-1]]
        nz_times2 = []
        nz_vals2 = []
        for tvi in range(len(s_err[2])):
            if s_err[2][tvi][1] > self.eps:
                nz_times2.append(s_err[2][tvi][0])
                nz_times2.append(s_err[2][tvi+1][0])
                nz_vals2.append(s_err[2][tvi][1])
                nz_vals2.append(s_err[2][tvi][1])
        s_vals2 = [v for (t, v) in s_err[2]]
        plt.hlines(s_vals2, s_times2, s_end_times2, colors='g')
        plt.vlines(nz_times2, np.zeros(len(nz_vals2)), nz_vals2, colors='g', linestyles='dotted')
        plt.xlabel("Zeit")
        plt.ylabel("Fehler")
        plt.xlim(0, 6)
        plt.ylim(0, 4 * self.bd_tol)
        plt.show()"""
        # erzeuge Ausgabe
        # OutputTableMulti(self.G, self.V, self.E, self.I, self.r, self.nu, self.fp, self.fp_ind, self.fm, self.q_global, self.q_ind, self.global_phase, self.labels, self.flow_vol, self.bd_tol)

        end_time = time.time()
        timediff1 = end_time - self.time_vor_main
        timediff2 = end_time - self.start_time
        print("times")
        print(timediff1)
        print(timediff2)

        """with open('output_examples/holzkirchen.txt', 'wb') as f:
            pickle.dump(self.fp, f)
            f.close()
        with open('output_examples/holzkirchen.txt', 'ab') as f:
            pickle.dump(self.fm, f)
            pickle.dump(self.q_global, f)
            pickle.dump(self.global_phase, f)
            f.close()"""

        """with open('output_examples/no_term200-6.txt', 'wb') as f:
            pickle.dump(self.fp, f)
            f.close()
        with open('output_examples/no_term200-6.txt', 'ab') as f:
            pickle.dump(self.fm, f)
            pickle.dump(self.q_global, f)
            pickle.dump(self.q_ind, f)
            pickle.dump(self.global_phase, f)
            f.close()"""

        """output_txt.write('\n f^+:')
        for i in range(self.I):
            for e_ind in range(self.m):
                # output_txt.write('\n ({0}, {1}): '.format(i, e_ind))
                output_txt.write('\n')
                for val in self.fp[i][e_ind]:
                    output_txt.write('{0}; '.format(val))

        output_txt.write('\n f^-:')
        for i in range(self.I):
            for e_ind in range(self.m):
                # output_txt.write('\n ({0}, {1}): '.format(i, e_ind))
                output_txt.write('\n')
                for val in self.fm[i][e_ind]:
                    output_txt.write('{0}; '.format(val))

        output_txt.write('\n q:')
        q_len = len(self.q_global)
        for e_ind in range(self.m):
            q_val = 0
            # output_txt.write('\n {0}: '.format(e_ind))
            output_txt.write('\n')
            for t in range(q_len):
                if self.q_global[t][e_ind] != q_val:
                    q_val = self.q_global[t][e_ind]
                    output_txt.write('({0}, {1}); '.format(self.global_phase[t], q_val))"""

        """output_json = open("output_examples/output-flow2.json", "w")
        output_json.write('{"network": {\n "nodes": [')
        output_json.close()
        output_json = open("output_examples/output-flow2.json", "a")
        for v_ind in range(self.n):
            output_json.write(' {{"id": {0}, "x": {1}, "y": 0.0}},'.format(v_ind, 0.0 + v_ind))
        output_json.close()
        with open("output_examples/output-flow2.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow2.json", "a")
        output_json.write('], \n "edges": [')
        for e_ind in range(self.m):
            v_ind = self.V.index(self.E[e_ind][0])
            w_ind = self.V.index(self.E[e_ind][1])
            output_json.write(' {{"id": {0}, "from": {1}, "to": {2}, "capacity": {3}, "transitTime": {4} }},'.format(e_ind, v_ind, w_ind, self.nu[e_ind], self.r[e_ind]))
        output_json.close()
        with open("output_examples/output-flow2.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow2.json", "a")
        output_json.write('], \n "commodities": [')
        colors = ["red", "blue", "green"]
        for i in range(self.I):
            output_json.write(' {{ "id": {0}, "color": "{1}" }},'.format(i, colors[i]))
        output_json.close()
        with open("output_examples/output-flow2.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow2.json", "a")
        output_json.write('] }, \n "flow": { \n "inflow": [')
        for e_ind in range(self.m):
            output_json.write('{')
            for i in range(self.I):
                output_json.write('"{0}": {{ \n "times": ['.format(i))
                for val in self.fp[i][e_ind]:
                    output_json.write(' {},'.format(val[0]))
                output_json.close()
                with open("output_examples/output-flow2.json", 'rb+') as oj:
                    oj.seek(-1, os.SEEK_END)
                    oj.truncate()
                    oj.close()
                output_json = open("output_examples/output-flow2.json", "a")
                output_json.write('], \n "values": [')
                for val in self.fp[i][e_ind]:
                    output_json.write(' {},'.format(val[1]))
                output_json.close()
                with open("output_examples/output-flow2.json", 'rb+') as oj:
                    oj.seek(-1, os.SEEK_END)
                    oj.truncate()
                    oj.close()
                output_json = open("output_examples/output-flow2.json", "a")
                output_json.write('] },')
            output_json.close()
            with open("output_examples/output-flow2.json", 'rb+') as oj:
                oj.seek(-1, os.SEEK_END)
                oj.truncate()
                oj.close()
            output_json = open("output_examples/output-flow2.json", "a")
            output_json.write('},')
        output_json.close()
        with open("output_examples/output-flow2.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow2.json", "a")
        output_json.write('], \n "outflow": [')
        for e_ind in range(self.m):
            output_json.write('{')
            for i in range(self.I):
                output_json.write('"{0}": {{ \n "times": ['.format(i))
                for val in self.fm[i][e_ind]:
                    output_json.write(' {},'.format(val[0]))
                output_json.close()
                with open("output_examples/output-flow2.json", 'rb+') as oj:
                    oj.seek(-1, os.SEEK_END)
                    oj.truncate()
                    oj.close()
                output_json = open("output_examples/output-flow2.json", "a")
                output_json.write('], \n "values": [')
                for val in self.fm[i][e_ind]:
                    output_json.write(' {},'.format(val[1]))
                output_json.close()
                with open("output_examples/output-flow2.json", 'rb+') as oj:
                    oj.seek(-1, os.SEEK_END)
                    oj.truncate()
                    oj.close()
                output_json = open("output_examples/output-flow2.json", "a")
                output_json.write('] },')
            output_json.close()
            with open("output_examples/output-flow2.json", 'rb+') as oj:
                oj.seek(-1, os.SEEK_END)
                oj.truncate()
                oj.close()
            output_json = open("output_examples/output-flow2.json", "a")
            output_json.write('},')
        output_json.close()
        with open("output_examples/output-flow2.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow2.json", "a")
        output_json.write('], \n "queues": [')
        for e_ind in range(self.m):
            output_json.write('{ "times": [')
            for t in self.q_ind[e_ind]:
                output_json.write(' {},'.format(self.global_phase[t]))
            output_json.close()
            with open("output_examples/output-flow2.json", 'rb+') as oj:
                oj.seek(-1, os.SEEK_END)
                oj.truncate()
                oj.close()

            output_json = open("output_examples/output-flow2.json", "a")
            output_json.write('], "values": [')
            for val in self.q_global[e_ind]:
                output_json.write(' {},'.format(val))
            output_json.close()
            with open("output_examples/output-flow2.json", 'rb+') as oj:
                oj.seek(-1, os.SEEK_END)
                oj.truncate()
                oj.close()

            output_json = open("output_examples/output-flow2.json", "a")
            output_json.write('], \n "domain": ["-Infinity", "Infinity"], "lastSlope": 0.0, "firstSlope": 0.0},')

        output_json.close()
        with open("output_examples/output-flow2.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()
        output_json = open("output_examples/output-flow2.json", "a")
        output_json.write('] } }')
        output_json.close()"""

        """nB = 7 * 12
        esB = list(range(12 * 12))
        for ano in range(11):
              esB.remove(4 + ano*12)
              esB.remove(8 + ano*12)
              esB.remove(11 + ano*12)
        esB.remove(3 +11*12)
        esB.remove(4 +11*12)
        esB.remove(8 +11*12)
        esB.remove(11 +11*12)

        output_json = open("output_examples/output-flow.json", "w")
        output_json.write('{"network": {\n "nodes": [')
        output_json.close()
        output_json = open("output_examples/output-flow.json", "a")
        for ano in range(12):
            output_json.write(' {{"id": {0}, "x": 0.0, "y": {1}}},'.format(0 + ano*7, 0.0 - 3 * ano))
            output_json.write(' {{"id": {0}, "x": 0.0, "y": {1}}},'.format(1 + ano*7, -2.0 - 3*ano))
            output_json.write(' {{"id": {0}, "x": -1.0, "y": {1}}},'.format(2 + ano*7, -2.0 - 3*ano))
            output_json.write(' {{"id": {0}, "x": -2.0, "y": {1}}},'.format(3 + ano*7, -1.0 - 3*ano))
            output_json.write(' {{"id": {0}, "x": -1.0, "y": {1}}},'.format(4 + ano*7, 0.0 - 3*ano))
            output_json.write(' {{"id": {0}, "x": 1.0, "y": {1}}},'.format(5 + ano*7, -2.0 - 3*ano))
            output_json.write(' {{"id": {0}, "x": 1.0, "y": {1}}},'.format(6 + ano*7, 0.0 - 3*ano))
        output_json.close()
        with open("output_examples/output-flow.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow.json", "a")
        output_json.write('], \n "edges": [')
        for e_ind in esB:
            v_ind = self.V.index(self.E[e_ind][0])
            w_ind = self.V.index(self.E[e_ind][1])
            output_json.write(' {{"id": {0}, "from": {1}, "to": {2}, "capacity": {3}, "transitTime": {4} }},'.format(esB.index(e_ind), v_ind, w_ind, self.nu[e_ind], self.r[e_ind]))
        output_json.close()
        with open("output_examples/output-flow.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow.json", "a")
        output_json.write('], \n "commodities": [')
        colors = ["red", "blue", "green"]
        for i in range(1):
            output_json.write(' {{ "id": {0}, "color": "{1}" }},'.format(i, colors[i]))
        output_json.close()
        with open("output_examples/output-flow.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow.json", "a")
        output_json.write('] }, \n "flow": { \n "inflow": [')
        for e_ind in esB:
            output_json.write('{')
            for i in range(1):
                output_json.write('"{0}": {{ \n "times": ['.format(i))
                for val in self.fp[i][e_ind]:
                    output_json.write(' {},'.format(val[0]))
                output_json.close()
                with open("output_examples/output-flow.json", 'rb+') as oj:
                    oj.seek(-1, os.SEEK_END)
                    oj.truncate()
                    oj.close()
                output_json = open("output_examples/output-flow.json", "a")
                output_json.write('], \n "values": [')
                for val in self.fp[i][e_ind]:
                    output_json.write(' {},'.format(val[1]))
                output_json.close()
                with open("output_examples/output-flow.json", 'rb+') as oj:
                    oj.seek(-1, os.SEEK_END)
                    oj.truncate()
                    oj.close()
                output_json = open("output_examples/output-flow.json", "a")
                output_json.write('] },')
            output_json.close()
            with open("output_examples/output-flow.json", 'rb+') as oj:
                oj.seek(-1, os.SEEK_END)
                oj.truncate()
                oj.close()
            output_json = open("output_examples/output-flow.json", "a")
            output_json.write('},')
        output_json.close()
        with open("output_examples/output-flow.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow.json", "a")
        output_json.write('], \n "outflow": [')
        for e_ind in esB:
            output_json.write('{')
            for i in range(1):
                output_json.write('"{0}": {{ \n "times": ['.format(i))
                for val in self.fm[i][e_ind]:
                    output_json.write(' {},'.format(val[0]))
                output_json.close()
                with open("output_examples/output-flow.json", 'rb+') as oj:
                    oj.seek(-1, os.SEEK_END)
                    oj.truncate()
                    oj.close()
                output_json = open("output_examples/output-flow.json", "a")
                output_json.write('], \n "values": [')
                for val in self.fm[i][e_ind]:
                    output_json.write(' {},'.format(val[1]))
                output_json.close()
                with open("output_examples/output-flow.json", 'rb+') as oj:
                    oj.seek(-1, os.SEEK_END)
                    oj.truncate()
                    oj.close()
                output_json = open("output_examples/output-flow.json", "a")
                output_json.write('] },')
            output_json.close()
            with open("output_examples/output-flow.json", 'rb+') as oj:
                oj.seek(-1, os.SEEK_END)
                oj.truncate()
                oj.close()
            output_json = open("output_examples/output-flow.json", "a")
            output_json.write('},')
        output_json.close()
        with open("output_examples/output-flow.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow.json", "a")
        q_times = {}
        q_vals = {}
        for e_ind in esB:
            q_times[e_ind] = [0]
            q_vals[e_ind] = [0]
        len_q_global = len(self.q_global)
        for t in range(len_q_global):
            for e_ind in esB:
                if self.q_global[t][e_ind] != q_vals[e_ind][-1]:
                    q_times[e_ind].append(self.global_phase[t])
                    q_vals[e_ind].append(self.q_global[t][e_ind])
        output_json.write('], \n "queues": [')
        for e_ind in esB:
            output_json.write('{ "times": [')
            for t in q_times[e_ind]:
                output_json.write(' {},'.format(t))
            output_json.close()
            with open("output_examples/output-flow.json", 'rb+') as oj:
                oj.seek(-1, os.SEEK_END)
                oj.truncate()
                oj.close()

            output_json = open("output_examples/output-flow.json", "a")
            output_json.write('], "values": [')
            for val in q_vals[e_ind]:
                output_json.write(' {},'.format(val))
            output_json.close()
            with open("output_examples/output-flow.json", 'rb+') as oj:
                oj.seek(-1, os.SEEK_END)
                oj.truncate()
                oj.close()

            output_json = open("output_examples/output-flow.json", "a")
            output_json.write('], \n "domain": ["-Infinity", "Infinity"], "lastSlope": 0.0, "firstSlope": 0.0},')

        output_json.close()
        with open("output_examples/output-flow.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()
        output_json = open("output_examples/output-flow.json", "a")
        output_json.write('] } }')
        output_json.close()"""

        """nB = 7
        esB = list(range(12))
        esB.remove(3)
        esB.remove(4)
        esB.remove(8)
        esB.remove(11)

        output_json = open("output_examples/output-flow.json", "w")
        output_json.write('{"network": {\n "nodes": [')
        output_json.close()
        output_json = open("output_examples/output-flow.json", "a")
        for ano in range(1):
            output_json.write(' {{"id": {0}, "x": 0.0, "y": {1}}},'.format(0 + ano*7, 0.0 - 3 * ano))
            output_json.write(' {{"id": {0}, "x": 0.0, "y": {1}}},'.format(1 + ano*7, -2.0 - 3*ano))
            output_json.write(' {{"id": {0}, "x": -1.0, "y": {1}}},'.format(2 + ano*7, -2.0 - 3*ano))
            output_json.write(' {{"id": {0}, "x": -2.0, "y": {1}}},'.format(3 + ano*7, -1.0 - 3*ano))
            output_json.write(' {{"id": {0}, "x": -1.0, "y": {1}}},'.format(4 + ano*7, 0.0 - 3*ano))
            output_json.write(' {{"id": {0}, "x": 1.0, "y": {1}}},'.format(5 + ano*7, -2.0 - 3*ano))
            output_json.write(' {{"id": {0}, "x": 1.0, "y": {1}}},'.format(6 + ano*7, 0.0 - 3*ano))
        output_json.close()
        with open("output_examples/output-flow.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow.json", "a")
        output_json.write('], \n "edges": [')
        for e_ind in esB:
            v_ind = self.V.index(self.E[e_ind][0])
            w_ind = self.V.index(self.E[e_ind][1])
            output_json.write(' {{"id": {0}, "from": {1}, "to": {2}, "capacity": {3}, "transitTime": {4} }},'.format(esB.index(e_ind), v_ind, w_ind, self.nu[e_ind], self.r[e_ind]))
        output_json.close()
        with open("output_examples/output-flow.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow.json", "a")
        output_json.write('], \n "commodities": [')
        colors = ["red", "blue", "green"]
        for i in range(1):
            output_json.write(' {{ "id": {0}, "color": "{1}" }},'.format(i, colors[i]))
        output_json.close()
        with open("output_examples/output-flow.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow.json", "a")
        output_json.write('] }, \n "flow": { \n "inflow": [')
        for e_ind in esB:
            output_json.write('{')
            for i in range(1):
                output_json.write('"{0}": {{ \n "times": ['.format(i))
                for val in self.fp[i][e_ind]:
                    output_json.write(' {},'.format(val[0]))
                output_json.close()
                with open("output_examples/output-flow.json", 'rb+') as oj:
                    oj.seek(-1, os.SEEK_END)
                    oj.truncate()
                    oj.close()
                output_json = open("output_examples/output-flow.json", "a")
                output_json.write('], \n "values": [')
                for val in self.fp[i][e_ind]:
                    output_json.write(' {},'.format(val[1]))
                output_json.close()
                with open("output_examples/output-flow.json", 'rb+') as oj:
                    oj.seek(-1, os.SEEK_END)
                    oj.truncate()
                    oj.close()
                output_json = open("output_examples/output-flow.json", "a")
                output_json.write('] },')
            output_json.close()
            with open("output_examples/output-flow.json", 'rb+') as oj:
                oj.seek(-1, os.SEEK_END)
                oj.truncate()
                oj.close()
            output_json = open("output_examples/output-flow.json", "a")
            output_json.write('},')
        output_json.close()
        with open("output_examples/output-flow.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow.json", "a")
        output_json.write('], \n "outflow": [')
        for e_ind in esB:
            output_json.write('{')
            for i in range(1):
                output_json.write('"{0}": {{ \n "times": ['.format(i))
                for val in self.fm[i][e_ind]:
                    output_json.write(' {},'.format(val[0]))
                output_json.close()
                with open("output_examples/output-flow.json", 'rb+') as oj:
                    oj.seek(-1, os.SEEK_END)
                    oj.truncate()
                    oj.close()
                output_json = open("output_examples/output-flow.json", "a")
                output_json.write('], \n "values": [')
                for val in self.fm[i][e_ind]:
                    output_json.write(' {},'.format(val[1]))
                output_json.close()
                with open("output_examples/output-flow.json", 'rb+') as oj:
                    oj.seek(-1, os.SEEK_END)
                    oj.truncate()
                    oj.close()
                output_json = open("output_examples/output-flow.json", "a")
                output_json.write('] },')
            output_json.close()
            with open("output_examples/output-flow.json", 'rb+') as oj:
                oj.seek(-1, os.SEEK_END)
                oj.truncate()
                oj.close()
            output_json = open("output_examples/output-flow.json", "a")
            output_json.write('},')
        output_json.close()
        with open("output_examples/output-flow.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow.json", "a")
        q_times = {}
        q_vals = {}
        for e_ind in esB:
            q_times[e_ind] = [0]
            q_vals[e_ind] = [0]
        len_q_global = len(self.q_global)
        for t in range(len_q_global):
            for e_ind in esB:
                if self.q_global[t][e_ind] != q_vals[e_ind][-1]:
                    q_times[e_ind].append(self.global_phase[t])
                    q_vals[e_ind].append(self.q_global[t][e_ind])
        output_json.write('], \n "queues": [')
        for e_ind in esB:
            output_json.write('{ "times": [')
            for t in q_times[e_ind]:
                output_json.write(' {},'.format(t))
            output_json.close()
            with open("output_examples/output-flow.json", 'rb+') as oj:
                oj.seek(-1, os.SEEK_END)
                oj.truncate()
                oj.close()

            output_json = open("output_examples/output-flow.json", "a")
            output_json.write('], "values": [')
            for val in q_vals[e_ind]:
                output_json.write(' {},'.format(val))
            output_json.close()
            with open("output_examples/output-flow.json", 'rb+') as oj:
                oj.seek(-1, os.SEEK_END)
                oj.truncate()
                oj.close()

            output_json = open("output_examples/output-flow.json", "a")
            output_json.write('], \n "domain": ["-Infinity", "Infinity"], "lastSlope": 0.0, "firstSlope": 0.0},')

        output_json.close()
        with open("output_examples/output-flow.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()
        output_json = open("output_examples/output-flow.json", "a")
        output_json.write('] } }')
        output_json.close()"""

        """output_json = open("output_examples/output-flow2.json", "w")
        output_json.write('{"network": {\n "nodes": [')
        output_json.close()
        output_json = open("output_examples/output-flow2.json", "a")
        for v_ind in range(4):
            output_json.write(' {{"id": {0}, "x": {1}, "y": 0.0}},'.format(v_ind, 0.0 + v_ind))
        output_json.close()
        with open("output_examples/output-flow2.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow2.json", "a")
        output_json.write('], \n "edges": [')
        for e_ind in range(3):
            v_ind = self.V.index(self.E[e_ind][0])
            w_ind = self.V.index(self.E[e_ind][1])
            output_json.write(' {{"id": {0}, "from": {1}, "to": {2}, "capacity": {3}, "transitTime": {4} }},'.format(e_ind, v_ind, w_ind, self.nu[e_ind], self.r[e_ind]))
        output_json.close()
        with open("output_examples/output-flow2.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow2.json", "a")
        output_json.write('], \n "commodities": [')
        colors = ["red", "blue", "green"]
        for i in range(self.I):
            output_json.write(' {{ "id": {0}, "color": "{1}" }},'.format(i, colors[i]))
        output_json.close()
        with open("output_examples/output-flow2.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow2.json", "a")
        output_json.write('] }, \n "flow": { \n "inflow": [')
        for e_ind in range(3):
            output_json.write('{')
            for i in range(self.I):
                output_json.write('"{0}": {{ \n "times": ['.format(i))
                for val in self.fp[i][e_ind]:
                    output_json.write(' {},'.format(val[0]))
                output_json.close()
                with open("output_examples/output-flow2.json", 'rb+') as oj:
                    oj.seek(-1, os.SEEK_END)
                    oj.truncate()
                    oj.close()
                output_json = open("output_examples/output-flow2.json", "a")
                output_json.write('], \n "values": [')
                for val in self.fp[i][e_ind]:
                    output_json.write(' {},'.format(val[1]))
                output_json.close()
                with open("output_examples/output-flow2.json", 'rb+') as oj:
                    oj.seek(-1, os.SEEK_END)
                    oj.truncate()
                    oj.close()
                output_json = open("output_examples/output-flow2.json", "a")
                output_json.write('] },')
            output_json.close()
            with open("output_examples/output-flow2.json", 'rb+') as oj:
                oj.seek(-1, os.SEEK_END)
                oj.truncate()
                oj.close()
            output_json = open("output_examples/output-flow2.json", "a")
            output_json.write('},')
        output_json.close()
        with open("output_examples/output-flow2.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow2.json", "a")
        output_json.write('], \n "outflow": [')
        for e_ind in range(3):
            output_json.write('{')
            for i in range(self.I):
                output_json.write('"{0}": {{ \n "times": ['.format(i))
                for val in self.fm[i][e_ind]:
                    output_json.write(' {},'.format(val[0]))
                output_json.close()
                with open("output_examples/output-flow2.json", 'rb+') as oj:
                    oj.seek(-1, os.SEEK_END)
                    oj.truncate()
                    oj.close()
                output_json = open("output_examples/output-flow2.json", "a")
                output_json.write('], \n "values": [')
                for val in self.fm[i][e_ind]:
                    output_json.write(' {},'.format(val[1]))
                output_json.close()
                with open("output_examples/output-flow2.json", 'rb+') as oj:
                    oj.seek(-1, os.SEEK_END)
                    oj.truncate()
                    oj.close()
                output_json = open("output_examples/output-flow2.json", "a")
                output_json.write('] },')
            output_json.close()
            with open("output_examples/output-flow2.json", 'rb+') as oj:
                oj.seek(-1, os.SEEK_END)
                oj.truncate()
                oj.close()
            output_json = open("output_examples/output-flow2.json", "a")
            output_json.write('},')
        output_json.close()
        with open("output_examples/output-flow2.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow2.json", "a")
        q_times = []
        q_vals = []
        for e_ind in range(3):
            q_times.append([0])
            q_vals.append([0])
        len_q_global = len(self.q_global)
        for t in range(len_q_global):
            for e_ind in range(3):
                if self.q_global[t][e_ind] != q_vals[e_ind][-1]:
                    q_times[e_ind].append(self.global_phase[t])
                    q_vals[e_ind].append(self.q_global[t][e_ind])
        output_json.write('], \n "queues": [')
        for e_ind in range(3):
            output_json.write('{ "times": [')
            for t in q_times[e_ind]:
                output_json.write(' {},'.format(t))
            output_json.close()
            with open("output_examples/output-flow2.json", 'rb+') as oj:
                oj.seek(-1, os.SEEK_END)
                oj.truncate()
                oj.close()

            output_json = open("output_examples/output-flow2.json", "a")
            output_json.write('], "values": [')
            for val in q_vals[e_ind]:
                output_json.write(' {},'.format(val))
            output_json.close()
            with open("output_examples/output-flow2.json", 'rb+') as oj:
                oj.seek(-1, os.SEEK_END)
                oj.truncate()
                oj.close()

            output_json = open("output_examples/output-flow2.json", "a")
            output_json.write('], \n "domain": ["-Infinity", "Infinity"], "lastSlope": 0.0, "firstSlope": 0.0},')

        output_json.close()
        with open("output_examples/output-flow2.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()
        output_json = open("output_examples/output-flow2.json", "a")
        output_json.write('] } }')
        output_json.close()"""

        """nB = 7 * 12 * 5
        esB = list(range(12 * 12 * 5))
        for bno in range(5):
            for ano in range(12 * bno, 12 * bno + 11):
                esB.remove(4 + ano * 12)
                esB.remove(8 + ano * 12)
                esB.remove(11 + ano * 12)
            esB.remove(3 + 11 * 12 + 12 * 12 * bno)
            esB.remove(4 + 11 * 12 + 12 * 12 * bno)
            esB.remove(8 + 11 * 12 + 12 * 12 * bno)
            esB.remove(11 + 11 * 12 + 12 * 12 * bno)

        output_json = open("output_examples/output-flow.json", "w")
        output_json.write('{"network": {\n "nodes": [')
        output_json.close()
        output_json = open("output_examples/output-flow.json", "a")
        for bno in range(5):
            for ano in range(12):
                output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(0 + ano*7 + bno*12*7, 0.0 + bno * 5, 0.0 - 3 * ano))
                output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(1 + ano*7 + bno*12*7, 0.0 + bno * 5, -2.0 - 3*ano))
                output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(2 + ano*7 + bno*12*7, -1.0 + bno * 5, -2.0 - 3*ano))
                output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(3 + ano*7 + bno*12*7, -2.0 + bno * 5, -1.0 - 3*ano))
                output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(4 + ano*7 + bno*12*7, -1.0 + bno * 5, 0.0 - 3*ano))
                output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(5 + ano*7 + bno*12*7, 1.0 + bno * 5, -2.0 - 3*ano))
                output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(6 + ano*7 + bno*12*7, 1.0 + bno * 5, 0.0 - 3*ano))
        output_json.close()
        with open("output_examples/output-flow.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow.json", "a")
        output_json.write('], \n "edges": [')
        for e_ind in esB:
            v_ind = self.V.index(self.E[e_ind][0])
            w_ind = self.V.index(self.E[e_ind][1])
            output_json.write(' {{"id": {0}, "from": {1}, "to": {2}, "capacity": {3}, "transitTime": {4} }},'.format(esB.index(e_ind), v_ind, w_ind, self.nu[e_ind], self.r[e_ind]))
        output_json.close()
        with open("output_examples/output-flow.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow.json", "a")
        output_json.write('], \n "commodities": [')
        colors = ["red", "blue", "green"]
        for i in range(1):
            output_json.write(' {{ "id": {0}, "color": "{1}" }},'.format(i, colors[i]))
        output_json.close()
        with open("output_examples/output-flow.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow.json", "a")
        output_json.write('] }, \n "flow": { \n "inflow": [')
        for e_ind in esB:
            output_json.write('{')
            for i in range(1):
                output_json.write('"{0}": {{ \n "times": ['.format(i))
                for val in self.fp[i][e_ind]:
                    output_json.write(' {},'.format(val[0]))
                output_json.close()
                with open("output_examples/output-flow.json", 'rb+') as oj:
                    oj.seek(-1, os.SEEK_END)
                    oj.truncate()
                    oj.close()
                output_json = open("output_examples/output-flow.json", "a")
                output_json.write('], \n "values": [')
                for val in self.fp[i][e_ind]:
                    output_json.write(' {},'.format(val[1]))
                output_json.close()
                with open("output_examples/output-flow.json", 'rb+') as oj:
                    oj.seek(-1, os.SEEK_END)
                    oj.truncate()
                    oj.close()
                output_json = open("output_examples/output-flow.json", "a")
                output_json.write('] },')
            output_json.close()
            with open("output_examples/output-flow.json", 'rb+') as oj:
                oj.seek(-1, os.SEEK_END)
                oj.truncate()
                oj.close()
            output_json = open("output_examples/output-flow.json", "a")
            output_json.write('},')
        output_json.close()
        with open("output_examples/output-flow.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow.json", "a")
        output_json.write('], \n "outflow": [')
        for e_ind in esB:
            output_json.write('{')
            for i in range(1):
                output_json.write('"{0}": {{ \n "times": ['.format(i))
                for val in self.fm[i][e_ind]:
                    output_json.write(' {},'.format(val[0]))
                output_json.close()
                with open("output_examples/output-flow.json", 'rb+') as oj:
                    oj.seek(-1, os.SEEK_END)
                    oj.truncate()
                    oj.close()
                output_json = open("output_examples/output-flow.json", "a")
                output_json.write('], \n "values": [')
                for val in self.fm[i][e_ind]:
                    output_json.write(' {},'.format(val[1]))
                output_json.close()
                with open("output_examples/output-flow.json", 'rb+') as oj:
                    oj.seek(-1, os.SEEK_END)
                    oj.truncate()
                    oj.close()
                output_json = open("output_examples/output-flow.json", "a")
                output_json.write('] },')
            output_json.close()
            with open("output_examples/output-flow.json", 'rb+') as oj:
                oj.seek(-1, os.SEEK_END)
                oj.truncate()
                oj.close()
            output_json = open("output_examples/output-flow.json", "a")
            output_json.write('},')
        output_json.close()
        with open("output_examples/output-flow.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()

        output_json = open("output_examples/output-flow.json", "a")
        q_times = {}
        q_vals = {}
        for e_ind in esB:
            q_times[e_ind] = [0]
            q_vals[e_ind] = [0]
        len_q_global = len(self.q_global)
        for t in range(len_q_global):
            for e_ind in esB:
                if self.q_global[t][e_ind] != q_vals[e_ind][-1]:
                    q_times[e_ind].append(self.global_phase[t])
                    q_vals[e_ind].append(self.q_global[t][e_ind])
        output_json.write('], \n "queues": [')
        for e_ind in esB:
            output_json.write('{ "times": [')
            for t in q_times[e_ind]:
                output_json.write(' {},'.format(t))
            output_json.close()
            with open("output_examples/output-flow.json", 'rb+') as oj:
                oj.seek(-1, os.SEEK_END)
                oj.truncate()
                oj.close()

            output_json = open("output_examples/output-flow.json", "a")
            output_json.write('], "values": [')
            for val in q_vals[e_ind]:
                output_json.write(' {},'.format(val))
            output_json.close()
            with open("output_examples/output-flow.json", 'rb+') as oj:
                oj.seek(-1, os.SEEK_END)
                oj.truncate()
                oj.close()

            output_json = open("output_examples/output-flow.json", "a")
            output_json.write('], \n "domain": ["-Infinity", "Infinity"], "lastSlope": 0.0, "firstSlope": 0.0},')

        output_json.close()
        with open("output_examples/output-flow.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()
        output_json = open("output_examples/output-flow.json", "a")
        output_json.write('] } }')
        output_json.close()"""

        nachwrite = time.time()
        writetime = nachwrite - end_time
        print("writetime", writetime)
        print("phases", len(self.global_phase))
        return 0
