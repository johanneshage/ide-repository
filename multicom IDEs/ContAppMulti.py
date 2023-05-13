import numpy as np
import scipy.sparse
import copy
from collections import defaultdict
import time
import pickle
# import matplotlib as mpl
import matplotlib.pyplot  as  plt
import os


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
        self.b = scipy.sparse.lil_matrix((self.I, self.n))
        self.items = G.items()
        self.keys = G.keys()
        self.eps = 10**(-12)  # Rechengenauigkeit
        self.bd_tol = 10**(-6)
        # self.flow_vol = []  # merke Flusswerte in den einzelnen Knoten für OutputTable # NUR FÜR OutputTable BENÖTIGT UND KOSTET SPEICHER
        self.delta_p = {}

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

        """self.E_active = np.ones((self.I, self.m))
        for i in range(self.I):
            self.labels.append(self.dijkstra(self.graphReversed, 't{}'.format(i+1), visited=[], distances={}))"""

        """with open('output_examples/holzkirchen_tau_dijkstra.txt', 'wb') as f:
            pickle.dump(self.labels[0], f)
            pickle.dump(self.labels[1], f)
            f.close()"""

        """with open('output_examples/no_term_dijkstra.txt', 'wb') as f:
            pickle.dump(self.labels[0], f)
            pickle.dump(self.labels[1], f)
            f.close()"""

        def loadall(filename):
            with open(filename, "rb") as f:
                while True:
                    try:
                        yield pickle.load(f)
                    except EOFError:
                        break

        """items = loadall('output_examples/holzkirchen_tau_dijkstra.txt')
        for item in items:
            self.labels.append(item.copy())"""

        items = loadall('output_examples/no_term_dijkstra.txt')
        for item in items:
            self.labels.append(item.copy())

        self.E_active = scipy.sparse.lil_matrix((self.I, self.m))
        for (v_ind, v) in enumerate(self.V):
            outneighbors = self.G[v].keys()
            self.delta_p[v_ind] = [self.E.index((v,u)) for u in outneighbors]
            for e_ind in self.delta_p[v_ind]:
                w_ind = self.V.index(self.E[e_ind][1])
                for i in range(self.I):
                    if abs(self.labels[i][v_ind] - self.labels[i][w_ind] - self.c[e_ind]) < self.eps:
                        self.E_active[i, e_ind] = 1

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
            if np.all([val == float('inf') for val in vals]):
                output = []
                for v in self.V:
                    if v in distances:
                        output.append(distances[v])
                    else:
                        output.append(float('inf'))
                return output
            src=min(unvisited, key=unvisited.get)

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

    '''def waterfilling_algo(self, v, flow_vol, bi, outneighbors, wi_slope):
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
        return z'''

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

    def get_outgoing_active_edges(self, i, v_ind):
        """
        bestimmt alle aus 'v' ausgehende, momentan für Gut 'i' aktive Kanten
        :param i: Index des Guts
        :param v_ind: Knoten
        :return: Liste der Indizes der aktiven Kanten
        """
        return [e for e in self.delta_p[v_ind] if self.E_active[i, e]]

    def change_of_cost(self, e_ind, z, mld):
        """
        Gibt die momentane Kostenänderung der Kante mit Index 'e_ind' bei Zufluss 'z' an.
        :param e_ind: Index der betrachteten Kante
        :param z: Einflussrate in Kante 'e'
        :return: Änderungsrate der Kosten
        """
        dif = z - self.nu[e_ind]
        if abs(dif) < self.bd_tol:  # * t1:
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
            if self.fm[i][e_ind][t][0] < theta + self.I * self.bd_tol:
                return t
        return 0

    def g(self, e_ind, x, mld):
        """
        Bestimmt absolute Änderung der Warteschlange von Kante 'e_ind' zum Zeitpunkt 'theta_ind' bei Einfluss 'x'.
        :param e_ind: Kantenindex
        :param x: momentaner Einfluss in Kante 'e_ind'
        :return: momentane Änderungsrate der Warteschlange
        """
        dif = x - self.nu[e_ind]
        if abs(dif) < self.bd_tol:  # * t1:
            return 0
        if self.q[e_ind] > self.eps:
            return dif
        return max([dif, 0])

    def gamma(self, x_sum, a, argmin, delta_p_act, coms, coms_lens):
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
                    if e_ind not in argmin[i][v_ind]:
                    #if abs(a[i, v_ind] - self.g(e_ind, x_sum[e_ind]) / self.nu[e_ind] - a[i, w_ind]) > self.bd_tol * coms_lens[v_ind]:
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

    def fp_approx(self, theta_ind,max_ld, top_ords, delta_p_act):
        theta = self.global_phase[-1]
        xk = scipy.sparse.lil_matrix((self.I, self.m))
        xfp = scipy.sparse.lil_matrix((self.I, self.m))
        x_sum = np.zeros(self.m)
        x_sum_fix = np.zeros(self.m)
        nu_sum_act = []
        coms = defaultdict(list)
        init_coms = defaultdict(list)
        coms_lens = np.zeros(self.n)  # 'coms_lens[v_ind]' entspricht der Länge von 'init_coms[v_ind]'
        # self.flow_vol.append(scipy.sparse.lil_matrix((self.n, 1)))  # NUR FÜR OutputTable BENÖTIGT UND KOSTET SPEICHER
        print(theta_ind, theta)

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
                    if flow_sent[i] > b_nf[i, v_ind] + self.eps:  # ???
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
                        a_e = self.g(e_ind, x_sum[e_ind] + x_sum_fix[e_ind], max_ld[v_ind]) / self.nu[e_ind] + a[i, w_ind]
                        a_e_tol = (self.bd_tol * coms_lens[v_ind]) / self.nu[e_ind]
                        if (i, v_ind) not in sparse_items:
                            if abs(a_e) > a_e_tol + self.eps:
                                a[i, v_ind] = a_e
                            sparse_items.append((i, v_ind))
                            argmin[i][v_ind] = [e_ind]
                            if e_ind in delta_p_act_nf[i][v_ind]:
                                argmin_nf[i][v_ind] = [e_ind]
                            else:
                                # min wird in bereits fixierter Kante angenommen. 'argmin_nf' enthält jedoch nur nicht-fixierte, minimale Kanten
                                argmin_nf[i][v_ind] = []
                        elif a_e < a[i, v_ind] + a_e_tol:
                            if a_e < a[i, v_ind] - a_e_tol:
                                a_min2[i, v_ind] = a[i, v_ind]
                                if abs(a_e) > a_e_tol + self.eps:
                                    a[i, v_ind] = a_e
                                else:
                                    a[i, v_ind] = 0
                                argmin[i][v_ind] = [e_ind]
                                if e_ind in delta_p_act_nf[i][v_ind]:
                                    argmin_nf[i][v_ind] = [e_ind]
                                else:
                                    # min wird in bereits fixierter Kante angenommen. 'argmin_nf' enthält jedoch nur nicht-fixierte, minimale Kanten
                                    argmin_nf[i][v_ind] = []
                            elif a_e < a[i, v_ind]:
                                if abs(a_e) > a_e_tol + self.eps:
                                    a[i, v_ind] = a_e
                                else:
                                    a[i, v_ind] = 0
                                argminiv = copy.deepcopy(argmin[i][v_ind])
                                for e_arg in argminiv:
                                    a_e_arg = self.g(e_arg, x_sum[e_arg] + x_sum_fix[e_arg], max_ld[v_ind]) / self.nu[e_arg] + a[i, self.V.index(self.E[e_arg][1])]
                                    if a_e_arg > a[i, v_ind] + a_e_tol:
                                        argmin[i][v_ind].remove(e_arg)
                                        if e_arg not in bounds_f[v_ind][i]:
                                            argmin_nf[i][v_ind].remove(e_arg)
                                        if a_e_arg < a_min2[i, v_ind]:
                                            a_min2[i, v_ind] = a_e_arg
                                argmin[i][v_ind].append(e_ind)
                                if e_ind in delta_p_act_nf[i][v_ind]:
                                    argmin_nf[i][v_ind].append(e_ind)
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
            """x_sum_m = np.zeros(self.m)
            for v_ind in coms:
                for i in coms[v_ind]:
                    if not argmin_nf[i][v_ind]:
                        for e_ind in argmin[i][v_ind]:
                            if xk[i, e_ind] > self.b[i, v_ind] - self.bd_tol * 2:
                                continue
                            x_sum_m[e_ind] += xk[i, e_ind]
            x_sum_nm = np.zeros(self.m)
            for v_ind in coms:
                for i in coms[v_ind]:
                    arg = list(set(delta_p_act[i][v_ind]) - set(argmin[i][v_ind]))
                    arg_nf = list(set(arg) - set(bounds_f[v_ind][i]))
                    if not arg_nf:
                        for e_ind in arg:
                            if xk[i, e_ind] < self.bd_tol * 2:
                                continue
                            x_sum_nm[e_ind] += xk[i, e_ind]"""
            for v_ind in coms:
                for i in coms[v_ind]:
                    arg = list(set(delta_p_act[i][v_ind]) - set(argmin[i][v_ind]))
                    arg_nf = list(set(arg) - set(bounds_f[v_ind][i]))
                    if not arg_nf:
                        # In diesem Fall sind alle aktiven nichtminimalen Kanten bereits fixiert. Damit können keine Fortschritte mehr gemacht werden -> relaxiere bds
                        for e_ind in arg:
                            if xk[i, e_ind] <= self.bd_tol * coms_lens[v_ind]:
                                continue
                            w_ind = self.V.index(self.E[e_ind][1])
                            if self.q[e_ind] < self.eps:
                                # Falls q = 0 -> g_e nicht-negativ -> kleinster a-Wert kann nicht erreicht werden -> wähle kleinstmögliche untere Schranke: 0
                                bounds[v_ind][i][e_ind] = (0, xk[i, e_ind])
                            else:
                                # relaxiere untere Schranke, sodass der minimale a-Wert erreicht werden kann
                                bounds[v_ind][i][e_ind] = (np.max([0, xk[i, e_ind] - (x_sum[e_ind] + x_sum_fix[e_ind] - (1 + a[i, v_ind] - a[i, w_ind]) * self.nu[e_ind])]), xk[i, e_ind])
                            bounds_f[v_ind][i].remove(e_ind)
                            delta_p_act_nf[i][v_ind].append(e_ind)
                            b_nf[i, v_ind] += xk[i, e_ind]
                    if not argmin_nf[i][v_ind]:
                        # Dagegen sind in diesem Fall alle minimalen Kanten fixiert. Damit können ebenfalls keine Forschritte gemacht werden -> relaxiere bds
                        for e_ind in argmin[i][v_ind]:
                            if xk[i, e_ind] > self.b[i, v_ind] - self.bd_tol * coms_lens[v_ind]:
                                continue
                            w_ind = self.V.index(self.E[e_ind][1])
                            xe = x_sum[e_ind] + x_sum_fix[e_ind]
                            x2 = (a_min2[i, v_ind] - a[i, w_ind]) * self.nu[e_ind] - self.g(e_ind, xe, max_ld[v_ind])
                            if self.q[e_ind] < self.eps and xe < self.nu[e_ind] - self.eps:
                                x2 += self.nu[e_ind] - xe
                            bounds_f[v_ind][i].remove(e_ind)
                            delta_p_act_nf[i][v_ind].append(e_ind)
                            # relaxierte obere Schranke mit x2: Dies ist der Flusswert, sodass der zweitkleinste a-Wert (= 'a_min2[i, v_ind]') auch für Kante 'e_ind'
                            # erreicht werden kann.
                            bounds[v_ind][i][e_ind] = (xk[i, e_ind], np.min([x2 + xk[i, e_ind], self.b[i, v_ind]]))
                            b_nf[i, v_ind] += xk[i, e_ind]
            return

        def fix_nodes():
            def dfs_step(visited, v_ind, i):
                for e_ind in delta_p_act[i][v_ind]:
                    w_ind = self.V.index(self.E[e_ind][1])
                    if w_ind not in visited:
                        dfs_step(visited, w_ind, i)
                visited.append(v_ind)

            def fix_check(v_ind):
                # Prüfe, ob alle aktiven Kanten den gleichen a-Wert liefern
                # ??? Eigentlich nur für Güter mit nichtdisjunkten Mengen delta_p_act_nf[i][v_ind] ?
                if np.all([v_ind not in argmin_nf[i] or len(argmin_nf[i][v_ind]) == len(delta_p_act_nf[i][v_ind]) for i in coms[v_ind]]):
                    # alle Kanten, für die ein Flusswert >0 fixiert wurde, müssen ebenfalls in 'argmin[i][v_ind]' liegen
                    if np.any([[e_ind for e_ind in bounds_f[v_ind][i] if xk[i, e_ind] > self.bd_tol * coms_lens[v_ind] and e_ind not in argmin[i][v_ind]] for i in coms[v_ind]]):
                        # sonst -> Abbruch
                        return False
                    for i in coms[v_ind]:
                        # ist ein Knoten der zwischen 'v_ind' und 'ti' liegt noch nicht fixiert -> Abbruch
                        visited = []
                        dfs_step(visited, v_ind, i)
                        visited.remove(v_ind)
                        if np.any([w_ind in coms for w_ind in visited]):
                            return False
                    # sonst: 'v_ind' kann fixiert werden
                    return True
                return False

            ci = 0
            coms_keys = list(coms)
            cn = len(coms_keys)
            while ci < cn:
                v_ind = coms_keys[ci]
                if fix_check(v_ind):
                    for i in coms[v_ind]:
                        for e_ind in delta_p_act[i][v_ind]:
                            if xk[i, e_ind] < self.bd_tol * coms_lens[v_ind]:
                                xfp[i, e_ind] = 0
                            elif xk[i, e_ind] > self.b[i, v_ind] - self.bd_tol * coms_lens[v_ind]:
                                xfp[i, e_ind] = self.b[i, v_ind]
                                x_sum_fix[e_ind] += xfp[i, e_ind]
                            else:
                                xfp[i, e_ind] = xk[i, e_ind]
                                x_sum_fix[e_ind] += xfp[i, e_ind]
                            x_sum[e_ind] -= xk[i, e_ind]
                            xk[i, e_ind] = 0
                        delta_p_act_nf[i][v_ind] = []
                    del coms[v_ind]
                    coms_keys.remove(v_ind)
                    ci = 0
                    cn -= 1
                else:
                    ci += 1
            return

        # Initialisiere x0
        for i in range(self.I):
            nu_sum_act.append({})
            for v_ind in range(self.n):
                if self.b[i, v_ind] > 0:
                    coms_lens[v_ind] += 1
                    # self.flow_vol[-1][v_ind, 0] += self.b_ot[-1][i, v_ind]
                    if len(delta_p_act[i][v_ind]) == 1:
                        e_ind = delta_p_act[i][v_ind][0]
                        xfp[i, e_ind] = self.b[i, v_ind]
                        x_sum_fix[e_ind] += xfp[i, e_ind]
                        init_coms[v_ind].append(i)
                        continue
                    coms[v_ind].append(i)
                    init_coms[v_ind].append(i)
                    nu_sum_act[i][v_ind] = 0
                    for e_ind in delta_p_act[i][v_ind]:
                        nu_sum_act[i][v_ind] += self.nu[e_ind]
                    for e_ind in delta_p_act[i][v_ind]:
                        xk[i, e_ind] = self.b[i, v_ind] * self.nu[e_ind] / nu_sum_act[i][v_ind]
                        x_sum[e_ind] += xk[i, e_ind]

        # Teilmenge von 'delta_p_act', welche für jedes Gut nur die Kanten mit nicht bereits fixiertem Flusswert enthält
        delta_p_act_nf = copy.deepcopy(delta_p_act)

        bounds = {}  # enthält untere und obere Schranken für 'xfp[i, e_ind]' - Werte
        bounds_f = {}
        for v_ind in coms:
            bounds[v_ind] = {}
            bounds_f[v_ind] = {}
            for i in coms[v_ind]:
                bounds[v_ind][i] = {}
                bounds_f[v_ind][i] = []
        b_nf = copy.deepcopy(self.b)
        fp_comp = []

        while True:
            # nur im ersten Schritt leer, danach enthält 'fp_comp' Tupel der Form '(i, v_ind)', was impliziert, dass die Aufteilung von Gut 'i' im Knoten 'v_ind'
            # vollständig fixiert wurde (d.h. die Boundsintervalle aller aktiven ausgehenden Kanten sind hinreichend klein)
            for (i, v_ind) in fp_comp:
                for e_ind in bounds[v_ind][i]:
                    if xk[i, e_ind] < self.bd_tol:  # ???
                        xfp[i, e_ind] = 0
                    elif xk[i, e_ind] > self.b[i, v_ind] - self.bd_tol:
                        xfp[i, e_ind] = self.b[i, v_ind]
                        x_sum_fix[e_ind] += xfp[i, e_ind]
                    else:
                        xfp[i, e_ind] = xk[i, e_ind]
                        x_sum_fix[e_ind] += xfp[i, e_ind]
                del bounds[v_ind][i]
                del bounds_f[v_ind][i]
                coms[v_ind].remove(i)
                if not coms[v_ind]:
                    del coms[v_ind]

            if fp_comp:  # nur im 1. Schritt nicht erfüllt. In diesem wird x0 verwendet
                if not coms:  # Abbruchbedingung: keine Güter mehr übrig, die noch auf Kanten aufgeteilt werden müssen
                    return xfp, x_sum_fix, a, argmin, init_coms, coms_lens
                xk, x_sum = calc_flow_by_bounds(xk_old=xk)

            a, a_min2, argmin, argmin_nf = calc_a()

            relax_bounds()

            fix_nodes()

            if not coms:  # Abbruchbedingung: keine Güter mehr übrig, die noch auf Kanten aufgeteilt werden müssen
                return xfp, x_sum_fix, a, argmin, init_coms, coms_lens

            while True:
                forced_zeros = self.gamma(x_sum + x_sum_fix, a, argmin, delta_p_act_nf, coms, coms_lens)  # gamma ???
                fp_comp = []
                new_f = []
                for v_ind in coms:
                    # coms_len = len(coms[v_ind])  # beachtet Löschen der bereits fixierten Güter (anders als 'coms_lens[v_ind]')
                    for i in coms[v_ind]:
                        dpanf_start = len(delta_p_act_nf[i][v_ind])
                        fz_nf = list(set(forced_zeros[i][v_ind]) - set(bounds_f[v_ind][i]))
                        for e_ind in fz_nf:
                            if e_ind in bounds[v_ind][i]:
                                bounds[v_ind][i][e_ind] = (bounds[v_ind][i][e_ind][0], xk[i, e_ind])
                            else:
                                # initiale untere Schranke: so, dass durch Verringerung von xk[i, e_ind] der für 'v_ind' minimale a_i - Wert erreicht werden kann
                                bounds[v_ind][i][e_ind] = (np.max([0, xk[i, e_ind] + (a[i, v_ind] - a[i, self.V.index(self.E[e_ind][1])]) * self.nu[e_ind] - self.g(e_ind, x_sum[e_ind] + x_sum_fix[e_ind], max_ld[v_ind])]), xk[i, e_ind])
                            bd_diff = abs(bounds[v_ind][i][e_ind][1] - bounds[v_ind][i][e_ind][0])
                            if bd_diff < self.bd_tol:
                                if bd_diff == 0 and xk[i, e_ind] > self.eps:
                                    # bd wird relaxiert
                                    bounds[v_ind][i][e_ind] = (0.5 * xk[i, e_ind], xk[i, e_ind])
                                else:
                                    delta_p_act_nf[i][v_ind].remove(e_ind)
                                    new_f.append((i, e_ind))
                        if v_ind not in argmin_nf[i] or len(argmin_nf[i][v_ind]) == dpanf_start:
                            # Anpassung der bounds nicht notwendig
                            continue
                        for e_ind in argmin_nf[i][v_ind]:
                            if e_ind in bounds[v_ind][i]:
                                bounds[v_ind][i][e_ind] = (xk[i, e_ind], bounds[v_ind][i][e_ind][1])
                            else:
                                w_ind = self.V.index(self.E[e_ind][1])
                                xe = x_sum[e_ind] + x_sum_fix[e_ind]
                                x2 = (a_min2[i, v_ind] - a[i, w_ind]) * self.nu[e_ind] - self.g(e_ind, xe, max_ld[v_ind])
                                if self.q[e_ind] < self.eps and xe < self.nu[e_ind] - self.eps:
                                    x2 += self.nu[e_ind] - xe
                                bounds[v_ind][i][e_ind] = (xk[i, e_ind], np.min([x2 + xk[i, e_ind], self.b[i, v_ind]]))
                            bd_diff = abs(bounds[v_ind][i][e_ind][1] - bounds[v_ind][i][e_ind][0])
                            if bd_diff < self.bd_tol:
                                if bd_diff == 0 and xk[i, e_ind] < self.b[i, v_ind] - self.eps:
                                    # bd wird relaxiert
                                    bounds[v_ind][i][e_ind] = (xk[i, e_ind], np.min([2 * xk[i, e_ind], self.b[i, v_ind]]))
                                else:
                                    delta_p_act_nf[i][v_ind].remove(e_ind)
                                    new_f.append((i, e_ind))

                xk, x_sum = calc_flow_by_bounds(xk_old=xk)
                a, a_min2, argmin, argmin_nf = calc_a()

                for (i, e_ind) in new_f:
                    v_ind = self.V.index(self.E[e_ind][0])
                    bounds_f[v_ind][i].append(e_ind)
                    b_nf[i, v_ind] -= xk[i, e_ind]
                    if len(delta_p_act_nf[i][v_ind]) == 1:
                        e_last = delta_p_act_nf[i][v_ind][0]
                        delta_p_act_nf[i][v_ind] = []
                        bounds_f[v_ind][i].append(e_last)
                        b_nf[i, v_ind] -= xk[i, e_last]
                    if b_nf[i, v_ind] < self.bd_tol * coms_lens[v_ind]:
                        if np.all([e in argmin[i][v_ind] for e in bounds_f[v_ind][i] if xk[i, e] > self.bd_tol * coms_lens[v_ind]]):
                            if (i, v_ind) not in fp_comp:  # ???
                                fp_comp.append((i, v_ind))
                        else:
                            for e in bounds_f[v_ind][i]:
                                delta_p_act_nf[i][v_ind].append(e)
                                if e in argmin[i][v_ind]:
                                    argmin_nf[i][v_ind].append(e)
                                    w_ind = self.V.index(self.E[e][1])
                                    xe = x_sum[e] + x_sum_fix[e]
                                    x2 = (a_min2[i, v_ind] - a[i, w_ind]) * self.nu[e] - self.g(e, xe, max_ld[v_ind])
                                    if self.q[e] < self.eps and xe < self.nu[e] - self.eps:
                                        x2 += self.nu[e] - xe
                                    bounds[v_ind][i][e] = (xk[i, e], np.min([x2 + xk[i, e], self.b[i, v_ind]]))
                                else:
                                    bounds[v_ind][i][e] = (np.max([0, xk[i, e] + (a[i, v_ind] - a[i, self.V.index(self.E[e][1])]) * self.nu[e] - self.g(e, x_sum[e] + x_sum_fix[e], max_ld[v_ind])]), xk[i, e])
                            bounds_f[v_ind][i] = []
                            b_nf[i, v_ind] = self.b[i, v_ind]

                if fp_comp:
                    # es existieren Knoten in 'fp_comp', welche fixiert werden können
                    break

                relax_bounds()

                fix_nodes()

                if not coms:
                    # Abbruchbedingung: keine Güter mehr übrig, die noch auf Kanten aufgeteilt werden müssen
                    return xfp, x_sum_fix, a, argmin, init_coms, coms_lens

    def main(self):
        """
        Hauptteil: Konstruiert schrittweise einen kontinuierlichen IDE-Fluss zu den Eingabedaten aus 'Math.cont_data'.
        :return: 0
        """
        theta = 0
        theta_ind = -1
        # Obergrenze für theta
        T = 250
        label_difs_act = scipy.sparse.lil_matrix((self.I, self.n))
        label_dif_sum = np.zeros(self.n)
        max_ld = np.zeros(self.n)
        max_diff = 0
        inflow_lens = np.zeros(self.m)
        flow_tol = np.ones(self.n)
        delta_p_act = []
        delta_p_inact = []
        top_ords = []
        # berechne 'delta_p_act', delta_p_inact', 'top_ords' und 'self.b' für 'theta' = 0
        for i in range(self.I):
            delta_p_act.append({})
            delta_p_inact.append({})
            for v_ind in range(self.n):
                delta_p_act[i][v_ind] = self.get_outgoing_active_edges(i, v_ind)
                delta_p_inact[i][v_ind] = [e for e in self.delta_p[v_ind] if e not in delta_p_act[i][v_ind]]
                self.b[i, v_ind] = self.calc_b(i, v_ind, theta)
            # berechne top. Sortierung für aktive Kanten von Gut 'i'
            top_ords.append(self.topologicalSort(i, delta_p_act[i]))
        skip_count = 0
        while theta < T:
            s_err = [[], []]
            s_err_rel = [[], []]
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
            skip_calcs = False
            if theta > 0:
                old_b = copy.deepcopy(self.b)
                self.b = scipy.sparse.lil_matrix((self.I, self.n))
                for i in range(self.I):
                    for v_ind in range(self.n):
                        self.b[i, v_ind] = self.calc_b(i, v_ind, theta)
                # prüfe, ob 'self.b' unverändert -> es besteht Möglichkeit, dass vorige Flusswerte immer noch korrekt
                if not (self.b - old_b).count_nonzero():
                    # berechne 'a' für vorige Flusswerte 'x_sum', mit neuen 'delta_p_act' und 'top_ords'
                    a = scipy.sparse.lil_matrix((self.I, self.n))
                    sparse_items = []
                    for i in range(self.I):
                        t_ind = self.V.index('t{}'.format(i+1))
                        sparse_items.append((i, t_ind))
                    argmin = []
                    for i in range(self.I):
                        argmin.append({})
                        for v_ind in top_ords[i]:
                            for e_ind in delta_p_act[i][v_ind]:
                                w_ind = self.V.index(self.E[e_ind][1])
                                a_e = self.g(e_ind, x_sum[e_ind], max_ld[v_ind]) / self.nu[e_ind] + a[i, w_ind]
                                a_e_tol = (self.bd_tol * coms_lens[v_ind]) / self.nu[e_ind]
                                if (i, v_ind) not in sparse_items:
                                    if abs(a_e) > a_e_tol + self.eps:
                                        a[i, v_ind] = a_e
                                    sparse_items.append((i, v_ind))
                                    argmin[i][v_ind] = [e_ind]
                                elif a_e < a[i, v_ind] + a_e_tol:
                                    if a_e < a[i, v_ind] - a_e_tol:
                                        if abs(a_e) > a_e_tol + self.eps:
                                            a[i, v_ind] = a_e
                                        else:
                                            a[i, v_ind] = 0
                                        argmin[i][v_ind] = [e_ind]
                                    elif a_e < a[i, v_ind]:
                                        if abs(a_e) > a_e_tol + self.eps:
                                            a[i, v_ind] = a_e
                                        else:
                                            a[i, v_ind] = 0
                                        argminiv = copy.deepcopy(argmin[i][v_ind])
                                        for e_arg in argminiv:
                                            a_e_arg = self.g(e_arg, x_sum[e_arg], max_ld[v_ind]) / self.nu[e_arg] + a[i, self.V.index(self.E[e_arg][1])]
                                            if a_e_arg > a[i, v_ind] + a_e_tol:
                                                argmin[i][v_ind].remove(e_arg)
                                        argmin[i][v_ind].append(e_ind)
                                    else:
                                        argmin[i][v_ind].append(e_ind)
                    dpact_pos_flow = {}
                    no_skip = False
                    for v_ind in init_coms:
                        dpact_pos_flow[v_ind] = {}
                        for i in init_coms[v_ind]:
                            dpact_pos_flow[v_ind][i] = [e_ind for e_ind in delta_p_act[i][v_ind] if x_total[i, e_ind] > self.bd_tol * coms_lens[v_ind]]
                            if np.any([e_ind not in dpact_pos_flow[v_ind][i] for e_ind in dp_act_old[i][v_ind] if x_total[i, e_ind] > self.bd_tol * coms_lens[v_ind]]):
                                no_skip = True
                                break
                        if no_skip:
                            break
                    # prüfe, ob vorige Flussaufteilung immer noch korrekt
                    if not no_skip and np.all([np.all([v_ind not in argmin[i] or np.all([e_ind in argmin[i][v_ind] for e_ind in dpact_pos_flow[v_ind][i]]) for i in init_coms[v_ind]]) for v_ind in init_coms]):
                        # verwende gleiche Flussaufteilung nochmal, 'a', 'argmin', 'top_ords', 'delta_p_act', 'delta_p_inact' bereits aktualisiert
                        skip_calcs = True
                        skip_count += 1
                        print("SKIP", theta_ind, theta)
                    else:
                        x_total, x_sum, a, argmin, init_coms, coms_lens = self.fp_approx(theta_ind, max_ld, top_ords, delta_p_act)
                else:
                    x_total, x_sum, a, argmin, init_coms, coms_lens = self.fp_approx(theta_ind, max_ld, top_ords, delta_p_act)
            else:
                 x_total, x_sum, a, argmin, init_coms, coms_lens = self.fp_approx(theta_ind, max_ld, top_ords, delta_p_act)

            delta_c = []
            for v_ind in range(self.n):
                for e_ind in self.delta_p[v_ind]:
                    delta_c.append(self.change_of_cost(e_ind, x_sum[e_ind], max_ld[v_ind]))

            for ti in range(self.I):
                for v_ind in top_ords[ti]:
                    # betrachte aktive Kanten
                    for e_ind in delta_p_act[ti][v_ind]:
                        if not skip_calcs:  # 'fp' bleibt gleich
                            if theta == 0:
                                if x_total[ti, e_ind] > self.bd_tol:
                                    self.fp[ti][e_ind][0] = (0, x_total[ti, e_ind])
                                    self.fp_ind[ti][e_ind].append(0)
                            # falls sich f^+ -Wert in dieser Phase ändert, aktualisiere 'self.fp'
                            elif abs(self.fp[ti][e_ind][-1][1] - x_total[ti, e_ind]) > 2 * self.bd_tol * coms_lens[v_ind]:
                                self.fp[ti][e_ind].append((theta, x_total[ti, e_ind]))
                                self.fp_ind[ti][e_ind].append(theta_ind)
                                if self.q_ind[e_ind][-1] != theta_ind and (self.q[e_ind] > self.eps or self.q_global[e_ind][-1] > self.eps or x_sum[e_ind] > self.nu[e_ind]):
                                    self.q_ind[e_ind].append(theta_ind)
                        if x_sum[e_ind] > 0:
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
                        if abs(self.fm[ti][e_ind][fm_ind][1] - outflow) > self.bd_tol * 2 * coms_lens[v_ind]:
                            # ???? braucht man das überhaupt noch?
                            if abs(self.fm[ti][e_ind][fm_ind][0] - outflow_time) < self.I * self.bd_tol:
                                del self.fm[ti][e_ind][fm_ind]
                            else:
                                self.fm[ti][e_ind].append((outflow_time, outflow))
            for ti in range(self.I):
                # überspringen jeweilige Senke, da diese hier uninteressant
                for v_ind in top_ords[ti][1:]:
                    for e_ind in delta_p_act[ti][v_ind]:
                        last_fm_ind = self.last_fm_change(ti, e_ind, theta)
                        # bestimme, ob sich vor dem Ende der aktuellen Phase ein f^- -Wert ändert -> verkürze Phase
                        if len(self.fm[ti][e_ind]) > last_fm_ind + 1 and self.I * self.bd_tol < self.fm[ti][e_ind][last_fm_ind + 1][0] - theta < next_phase[0] + max_diff:
                            next_fm = self.fm[ti][e_ind][last_fm_ind + 1][0] - theta
                            if next_fm < next_phase[0]:
                                nph_rev = copy.deepcopy(next_phase)
                                nph_rev.reverse()
                                for p in nph_rev:
                                    if p > next_fm + max_diff:
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
                        if not skip_calcs:
                            # falls f^+ -Wert vorher > 0 war, so wird dieser hier auf 0 gesetzt, da Kante inaktiv
                            if abs(self.fp[ti][e_ind][-1][1]) > self.eps:
                                self.fp[ti][e_ind].append((theta, 0))
                                self.fp_ind[ti][e_ind].append(theta_ind)
                                if self.q_ind[e_ind][-1] != theta_ind and (self.q[e_ind] > self.eps or self.q_global[e_ind][-1] > self.eps or x_sum[e_ind] > self.nu[e_ind] + self.eps):
                                    self.q_ind[e_ind].append(theta_ind)
                                if self.q[e_ind] > self.eps:
                                    outflow_time = theta + self.q[e_ind]/self.nu[e_ind] + self.r[e_ind]
                                else:
                                    outflow_time = theta + self.r[e_ind]
                                fm_ind = self.last_fm_change(ti, e_ind, outflow_time)
                                # setze 'self.fm'- Wert auf 0, falls dies noch nicht geschehen ist
                                if abs(self.fm[ti][e_ind][fm_ind][1]) > self.I * self.bd_tol:
                                    self.fm[ti][e_ind].append((outflow_time, 0))

                        # bestimme, ob sich vor dem Ende der aktuellen Phase ein f^- -Wert ändert -> verkürze Phase
                        last_fm_ind = self.last_fm_change(ti, e_ind, theta)
                        if len(self.fm[ti][e_ind]) > last_fm_ind + 1 and self.I * self.bd_tol < self.fm[ti][e_ind][last_fm_ind + 1][0] - theta < next_phase[0] + max_diff:
                            next_fm = self.fm[ti][e_ind][last_fm_ind + 1][0] - theta
                            if next_fm < next_phase[0]:
                                nph_rev = copy.deepcopy(next_phase)
                                nph_rev.reverse()
                                for p in nph_rev:
                                    if p > next_fm + max_diff:
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
                        if x_total[ti, active_ind] > 0:
                            break
                    if len_act:
                        for e_ind in delta_p_inact[ti][v_ind]:
                            # prüfe, wann inaktive Kanten unter momentanem Einfluss aktiv werden
                            tar_ind = self.V.index(self.E[e_ind][1])
                            act_ind = self.V.index(self.E[active_ind][1])
                            np_tol = np.max([self.bd_tol, self.bd_tol * coms_lens[v_ind] * next_phase[0]])

                            if self.labels[ti][tar_ind] + self.q[e_ind]/self.nu[e_ind] + self.r[e_ind] + \
                                    (delta_c[e_ind] + a[ti, tar_ind]) * (next_phase[0] + np_tol) < \
                                    self.labels[ti][act_ind] + self.q[active_ind]/self.nu[active_ind] + \
                                    self.r[active_ind] + (delta_c[active_ind] + a[ti, act_ind]) * (next_phase[0] + np_tol) and \
                                    abs(delta_c[e_ind] + a[ti, tar_ind] - delta_c[active_ind] - a[ti, act_ind]) \
                                    > self.eps:
                                time_ub = (self.labels[ti][act_ind] + self.q[active_ind]/self.nu[active_ind] + self.r[active_ind] - \
                                          self.labels[ti][tar_ind] - self.q[e_ind]/self.nu[e_ind] - self.r[e_ind]) \
                                          / (delta_c[e_ind] + a[ti, tar_ind] - delta_c[active_ind] - a[ti, act_ind])

                                if time_ub < next_phase[0]:
                                    nph_rev = copy.deepcopy(next_phase)
                                    nph_rev.reverse()
                                    for p in nph_rev:
                                        if p > time_ub + np.max([self.bd_tol * coms_lens[v_ind] * time_ub, max_diff, self.bd_tol]):
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
            for v_ind in range(self.n):
                for e_ind in self.delta_p[v_ind]:
                    change_of_q = delta_c[e_ind] * self.nu[e_ind]
                    # Falls die Warteschlange von 'e' unter aktuellem Fluss abgebaut wird, bestimme Zeitpunkt, zu dem diese vollständig abgebaut ist (bei gleich bleibendem
                    # Fluss)
                    if change_of_q < -self.eps:
                        # 'phase_length': Dauer bis Warteschlangenlänge gleich 0
                        phase_length = - self.q[e_ind] / change_of_q
                        if phase_length < next_phase[0] + self.bd_tol * coms_lens[v_ind] * next_phase[0] + max_diff / (-change_of_q + self.bd_tol * coms_lens[v_ind]):
                            if phase_length < next_phase[0]:
                                nph_rev = copy.deepcopy(next_phase)
                                nph_rev.reverse()
                                for p in nph_rev:
                                    if p > phase_length + self.bd_tol * coms_lens[v_ind] * phase_length + max_diff / (-change_of_q + self.bd_tol * coms_lens[v_ind]):
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

            if not skip_calcs and theta_ind > 0:
                for e_ind in range(self.m):
                    if np.any([theta_ind in self.fp_ind[i][e_ind] for i in range(self.I)]) and (self.q[e_ind] > self.eps or self.q_global[e_ind][-1] > self.eps or x_sum[e_ind] > self.nu[e_ind]):
                        self.q_global[e_ind].append(self.q[e_ind])

            for v_ind in range(self.n):
                for e_ind in self.delta_p[v_ind]:
                    if np.any([self.fp[i][e_ind][-1][1] > 0 for i in range(self.I)]):
                        inflow_lens[e_ind] += next_phase[-1]
                        if inflow_lens[e_ind] > flow_tol[v_ind]:
                            flow_tol[v_ind] = inflow_lens[e_ind]

            theta += next_phase[-1]
            if next_phase[-1] != T:
                # aktualisiere Warteschlangenlängen und Kosten
                new_q = []
                for v_ind in range(self.n):
                    q_tol = np.max([self.bd_tol, self.bd_tol * coms_lens[v_ind] * next_phase[-1]])
                    for e_ind in self.delta_p[v_ind]:
                        next_q_len = self.q[e_ind] + self.g(e_ind, x_sum[e_ind], max_ld[v_ind]) * next_phase[-1]
                        if next_q_len < q_tol:
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
                label_difs = scipy.sparse.lil_matrix((self.I, self.m))
                max_diff = max([max_ld[v] for v in range(self.n)])
                for i in range(self.I):
                    for v_ind in range(self.n):
                        self.labels[i][v_ind] = self.labels[i][v_ind] + next_phase[-1] * a[i, v_ind]
                    # etwas effizienter als 'v_ind in range(self.n)', da t_i übersprungen wird. Aber: 'top_ords[i]' hier keine aktuelle top. Sortierung mehr
                    for v_ind in top_ords[i][1:]:
                        if np.isinf(self.labels[i][v_ind]):
                            continue
                        v = self.V[v_ind]
                        outneighbors = self.G[v].keys()
                        #dv_act = []
                        for w in outneighbors:
                            w_ind = self.V.index(w)
                            edge = self.E.index((v, w))
                            label_difs[i, edge] = self.labels[i][v_ind] - self.labels[i][w_ind] - self.c[edge]
                            if label_difs[i, edge] > 0 or abs(label_difs[i, edge]) < max_diff + self.bd_tol or abs(label_difs[i, edge]) < 2 * max_diff + 2 * self.I * self.bd_tol * next_phase[-1]:
                                # 'label_dif' > 0: geschieht nur aufgrund von Approximationsfehlern, für Kanten die aktiv sein sollten
                                self.E_active[i, edge] = 1
                                #dv_act.append(edge)
                        """if len(dv_act) == 0:
                            print()"""

                """for i in range(self.I):
                    err_i = 0
                    err_i_rel = 0
                    for v_ind in range(self.n):
                        if self.b[i, v_ind] > 0:
                            #delta_p_v_act = self.get_outgoing_active_edges(i, v_ind)
                            if len(delta_p_act[i][v_ind]) == 1:
                                continue
                            vmax = -np.Inf
                            vmin = np.Inf
                            for e_ind in delta_p_act[i][v_ind]:
                                if self.fp[i][e_ind][-1][1] > 0:
                                    w_ind = self.V.index(self.E[e_ind][1])
                                    val = self.labels[i][w_ind] + self.c[e_ind]
                                    if val > vmax:
                                        vmax = val
                                    if val < vmin:
                                        vmin = val
                            err_i += vmax - vmin
                            err_i_rel += (vmax - vmin) / self.b[i, v_ind]
                    s_err[i].append((theta, err_i))
                    s_err_rel[i].append((theta, err_i_rel))
                with open('output_examples/no_term_errors250-6.txt', 'ab') as f:
                    pickle.dump(s_err, f)
                    pickle.dump(s_err_rel, f)
                    f.close()"""

                dp_act_old = copy.deepcopy(delta_p_act)
                for i in range(self.I):
                    for v_ind in range(self.n):
                        delta_p_act[i][v_ind] = self.get_outgoing_active_edges(i, v_ind)
                        delta_p_inact[i][v_ind] = [e for e in self.delta_p[v_ind] if e not in delta_p_act[i][v_ind]]
                    # berechne top. Sortierung zur aktualisierten Menge aktiver Kanten für Gut 'i'
                    top_ords[i] = self.topologicalSort(i, delta_p_act[i])
                    for v_ind in top_ords[i][1:]:
                        """len_dpaiv = len(delta_p_act[i][v_ind])
                        if len_dpaiv == 1:
                            # labels werden angepasst, damit diese durch Approximationsfehler nicht zu stark auseinanderdriften.
                            e = self.E[delta_p_act[i][v_ind][0]]
                            w_ind = self.V.index(e[1])
                            #labelvor = self.labels[i][v_ind]
                            self.labels[i][w_ind] = self.labels[i][v_ind] - self.c[delta_p_act[i][v_ind][0]]
                            if self.labels[i][w_ind] < self.bd_tol * flow_tol[v_ind]:
                                self.labels[i][w_ind] = 0
                            #self.labels[i][v_ind] = self.labels[i][w_ind] + self.c[delta_p_act[i][v_ind][0]]
                            #labelnach = self.labels[i][v_ind]
                            #qi[dv_act[v_ind][0]][i] = self.q[dv_act[v_ind][0]]
                        elif len_dpaiv > 1:
                            #l_sum = 0"""

                        for e_ind in delta_p_act[i][v_ind]:
                            w_ind = self.V.index(self.E[e_ind][1])
                            label_dif = self.labels[i][v_ind] - self.labels[i][w_ind] - self.c[e_ind]
                            if abs(label_dif) > label_difs_act[i, v_ind]:
                                label_dif_sum[v_ind] += abs(label_dif) - label_difs_act[i, v_ind]
                                label_difs_act[i, v_ind] = abs(label_dif)
                                if label_difs_act[i, v_ind] > max_ld[v_ind]:
                                    max_ld[v_ind] = label_difs_act[i, v_ind]

                        for e_ind in delta_p_act[i][v_ind]:
                            w_ind = self.V.index(self.E[e_ind][1])
                            self.labels[i][w_ind] = self.labels[i][v_ind] - self.c[e_ind]
                            if self.labels[i][w_ind] < label_difs_act[i, v_ind]:
                                self.labels[i][w_ind] = 0

                        for e_ind in delta_p_act[i][v_ind]:
                            w_ind = self.V.index(self.E[e_ind][1])
                            # 'qi[e_ind][i]': Warteschlangenlönge von Kante 'e_ind', passend zu den labels von Gut 'i'
                            qi[e_ind][i] = (self.labels[i][v_ind] - self.labels[i][w_ind] - self.r[e_ind]) * self.nu[e_ind]
                            if qi[e_ind][i] < label_difs_act[i, v_ind] * self.nu[e_ind] + self.eps:
                                qi[e_ind][i] = 0

                for e_ind in range(self.m):
                    len_qie = len(qi[e_ind])
                    if len_qie > 0:
                        mean_qi = 0
                        for i in qi[e_ind]:
                            mean_qi += qi[e_ind][i]
                        mean_qi /= len_qie
                        self.q[e_ind] = mean_qi
                        self.c[e_ind] = self.q[e_ind] / self.nu[e_ind] + self.r[e_ind]

                for p in fm_to_zero:
                    for e_ind in fm_to_zero[p]:
                        for i in range(self.I):
                            if self.fm[i][e_ind][-1][1] > self.eps:
                                self.fm[i][e_ind].append((theta + self.r[e_ind], 0))

        for i in range(self.I):
            # am Ende sind alle f^+ -Werte 0
            for e in range(self.m):
                if abs(self.fp[i][e][-1][1]) > self.eps:
                    self.fp[i][e].append((theta - next_phase[-1], 0))
                    self.fp_ind[i][e].append(theta_ind)

        """s_times0 = [t for (t, v) in s_err[0]]
        s_vals0 = [v for (t, v) in s_err[0]]
        #plt.plot(s_times0, s_vals0, 'r')
        s_times1 = [t for (t, v) in s_err[1]]
        s_vals1 = [v for (t, v) in s_err[1]]
        #plt.plot(s_times1, s_vals1, 'b')
        s_times2 = [t for (t, v) in s_err[2]]
        s_vals2 = [v for (t, v) in s_err[2]]
        #plt.plot(s_times2, s_vals2, 'g')
        s_vals_rel0 = [v for (t, v) in s_err_rel[0]]
        s_vals_rel1 = [v for (t, v) in s_err_rel[1]]
        s_vals_rel2 = [v for (t, v) in s_err_rel[2]]
        total_err = []
        rel_err = []
        for t in range(len(s_times0)):
            total_err.append(s_vals0[t] + s_vals1[t] + s_vals2[t])
            rel_err.append(s_vals_rel0[t] + s_vals_rel1[t] + s_vals_rel2[t])
        #plt.plot(s_times0, total_err, 'k')
        fig, axs = plt.subplots(2)
        # fig.suptitle('Vertically stacked subplots')
        axs[0].plot(s_times0, s_vals0, 'r')
        axs[0].plot(s_times0, s_vals1, 'b')
        axs[0].plot(s_times0, s_vals2, 'g')
        axs[1].plot(s_times0, rel_err, 'c')
        axs[1].plot(s_times0, total_err, 'k')
        axs[1].legend(['relativ', 'absolut'])
        axs[0].set(ylabel='Fehler', xlim=(0, 7), ylim=0)
        axs[1].set(xlabel='Zeit', ylabel='Fehler', xlim=(0, 7), ylim=0)
        plt.show()"""

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

        end_time = time.time()
        timediff1 = end_time - self.time_vor_main
        timediff2 = end_time - self.start_time
        print("times")
        print(timediff1)
        print(timediff2)
        print(timediff2)

        """with open('output_examples/holzkirchen_komplett-6.txt', 'ab') as f:
            pickle.dump(self.fp, f)
            pickle.dump(self.fm, f)
            pickle.dump(self.q_global, f)
            pickle.dump(self.q_ind, f)
            pickle.dump(self.global_phase, f)
            f.close()"""

        with open('output_examples/no_term250-6.txt', 'ab') as f:
            pickle.dump(self.fp, f)
            pickle.dump(self.fm, f)
            pickle.dump(self.q_global, f)
            pickle.dump(self.q_ind, f)
            pickle.dump(self.global_phase, f)
            f.close()

        """with open('output_examples/no_term200-6.txt', 'wb') as f:
            pickle.dump(self.fp, f)
            f.close()
        with open('output_examples/no_term200-6.txt', 'ab') as f:
            pickle.dump(self.fm, f)
            pickle.dump(self.q_global, f)
            pickle.dump(self.q_ind, f)
            pickle.dump(self.global_phase, f)
            f.close()"""

        '''output_json = open("output_examples/output-flow2.json", "w")
        output_json.write('{"network": {\n "nodes": [')
        output_json.close()
        output_json = open("output_examples/output-flow2.json", "a")
        """for v_ind in range(self.n):
           output_json.write(' {{"id": {0}, "x": {1}, "y": 0.0}},'.format(v_ind, 0.0 + v_ind))"""
        output_json.write(' {"id": 0, "x": 0, "y": 0},')
        output_json.write(' {"id": 1, "x": 1.5, "y": -1.5},')
        output_json.write(' {"id": 2, "x": 1.5, "y": 0},')
        output_json.write(' {"id": 3, "x": 1.5, "y": 1.5},')
        output_json.write(' {"id": 4, "x": 3, "y": -2.5},')
        output_json.write(' {"id": 5, "x": 3, "y": -1},')
        output_json.write(' {"id": 6, "x": 3, "y": 0.5},')
        output_json.write(' {"id": 7, "x": 3, "y": 2},')
        output_json.write(' {"id": 8, "x": 4.5, "y": -1},')
        output_json.write(' {"id": 9, "x": 4.5, "y": 1.25},')
        output_json.write(' {"id": 10, "x": 6, "y": -2},')
        output_json.write(' {"id": 11, "x": 6, "y": 0},')
        output_json.write(' {"id": 12, "x": 8, "y": 1.25},')
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
        output_json.close()'''

        nachwrite = time.time()
        writetime = nachwrite - end_time
        print("writetime", writetime)
        print("phases", len(self.global_phase))
        print("davon geskipt:", skip_count)
        return 0
