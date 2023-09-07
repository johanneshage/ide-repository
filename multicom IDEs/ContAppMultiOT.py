import numpy as np
import scipy.sparse
import copy
from collections import defaultdict
import time
import pickle
from OutputTableMulti import OutputTableMulti
# import os


class ContAppMultiOT:
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
        self.b = [scipy.sparse.lil_matrix((self.I, self.n))]
        self.items = G.items()
        self.keys = G.keys()
        self.eps = 10**(-13)  # Rechengenauigkeit
        self.bd_tol = 10**(-5)  # Approximationsgenauigkeit
        self.flow_vol = []  # merke Flusswerte in den einzelnen Knoten für OutputTable # NUR FÜR OutputTable BENÖTIGT UND KOSTET SPEICHER
        self.delta_p = {}

        for delta in self.items:
            for w in list(delta[1].keys()):
                self.E.append((delta[0], w))  # speichere alle Kanten in E
                self.r.append(list(self.items)[list(self.keys).index(delta[0])][1][w][1])  # speichere Reisezeiten in r
                self.nu.append(list(self.items)[list(self.keys).index(delta[0])][1][w][0])  # Kapazitäten in nu
        self.m = len(self.E)  # Anz. Kanten
        self.c = []
        self.global_phase = [0]  # speichere Startzeitpunkte der "globalen" Phasen, also alle Zeitpunkte, zu denen sich mindestens eine lokale Phase ändert

        self.fp = [[] for _ in range(self.I)]  # f^+    [[],[]]
        self.fp_ind = [[] for _ in range(self.I)]  # Liste der Indizes der Phasen
        self.fm = [[] for _ in range(self.I)]  # f^-
        self.q_global = [[0] for _ in range(self.m)]  # enthält für jede Kante alle 'piecewise_linear' Eckpunkte der Warteschlangen
        self.q_ind = [[0] for _ in range(self.m)]  # enthält die zu diesen Eckpunkten zugehörigen Phasen
        self.q = np.zeros(self.m)  # enthält alle aktuellen Warteschalngenlängen
        for i in range(self.I):
            for k in range(self.m):
                self.fp[i].append([(0, 0)])
                self.fp_ind[i].append([])
                self.fm[i].append([(0, 0)])

        self.c.append(np.zeros(self.m))
        for e_ind in range(self.m):  # Initialisierung "self.c" (Kosten)
            self.c[0][e_ind] = self.r[e_ind]

        # Zeitpunkte, zu denen sich der externe Zufluss in mindestens einem Quellknoten ändert
        self.u_start = set()
        for com in self.u:
            for s_list in com:
                for t in s_list:
                    self.u_start.add(t[0])

        self.graphReversed = self.reverse_graph(G)

        self.E_active = np.ones((self.I, self.m))
        for i in range(self.I):
            self.labels.append([self.dijkstra(self.graphReversed, 't{}'.format(i+1), visited=[], distances={})])

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

        """items = loadall('output_examples/no_term_dijkstra.txt')
        for item in items:
            self.labels.append(item.copy())"""

        """items = loadall('output_examples/new_no_term40-5.txt')

        for (no, item) in enumerate(items):
            if no == 0:
                self.fp = item.copy()
            elif no == 1:
                self.fm = item.copy()
            elif no == 2:
                self.q_global = item.copy()
            elif no == 3:
                self.q_ind = item.copy()
            elif no == 4:
                self.global_phase = item.copy()
            elif no == 5:
                self.E_active = item.copy()
            elif no == 6:
                theta = item.copy()
            elif no == 7:
                self.b = item.copy()
            elif no == 8:
                self.labels = item.copy()
            elif no == 9:
                self.q = item.copy()
            elif no == 10:
                self.c = item.copy()
            elif no == 11:
                x_total = item.copy()
            elif no == 12:
                x_sum = item.copy()
            elif no == 13:
                coms_lens = item.copy()
            elif no == 14:
                init_coms = item.copy()
            elif no == 15:
                dp_act_old = item.copy()
            elif no == 16:
                max_dif = item
            elif no == 17:
                delta_p_act = item.copy()
            elif no == 18:
                delta_p_inact = item.copy()
            elif no == 19:
                top_ords = item.copy()
            else:
                nu_min = item.copy()"""

        self.E_active = scipy.sparse.lil_matrix((self.I, self.m))
        for (v_ind, v) in enumerate(self.V):
            outneighbors = self.G[v].keys()
            self.delta_p[v_ind] = [self.E.index((v,u)) for u in outneighbors]
            for e_ind in self.delta_p[v_ind]:
                w_ind = self.V.index(self.E[e_ind][1])
                for i in range(self.I):
                    if abs(self.labels[i][-1][v_ind] - self.labels[i][-1][w_ind] - self.c[-1][e_ind]) < self.eps:
                        self.E_active[i, e_ind] = 1

        self.time_vor_main = time.time()
        self.main()
        # self.main(theta, x_total, x_sum, coms_lens, init_coms, dp_act_old, max_dif, delta_p_act, delta_p_inact, top_ords, nu_min)

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
            if np.isinf(self.labels[i][-1][k]):
                continue
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
                    new_distance = distances[src] + self.c[-1][self.E.index((neighbor,src))]
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

    def get_outgoing_active_edges(self, i, v_ind):
        """
        bestimmt alle aus 'v' ausgehende, momentan für Gut 'i' aktive Kanten
        :param i: Index des Guts
        :param v_ind: Knoten
        :return: Liste der Indizes der aktiven Kanten
        """
        return [e for e in self.delta_p[v_ind] if self.E_active[i, e]]

    def change_of_cost(self, e_ind, z):
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

    def g(self, e_ind, x):
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

    def gamma(self, argmin, delta_p_act, coms):
        """
        Erzeugt Liste 'forced_zeros': 'forced_zeroes[i][v_ind]' enthält alle für 'i' aktiven, von 'v_ind' ausgehenden, Kanten, deren Flusswert von 'i' 0 sein muss.
        :param argmin:
        :param delta_p_act:
        :param coms:
        :return:
        """
        forced_zeros = []
        for i in range(self.I):
            forced_zeros.append({})
        for v_ind in coms:
            for i in coms[v_ind]:
                forced_zeros[i][v_ind] = []
                for e_ind in delta_p_act[i][v_ind]:
                    # w_ind = self.V.index(self.E[e_ind][1])
                    if e_ind not in argmin[i][v_ind]:
                        forced_zeros[i][v_ind].append(e_ind)
        return forced_zeros

    '''def get_items(self, s):
        """
        Bestimmt Indizes der gesetzen Werte der (dünnbesetzten) Matrix 's'
        :param s: Matrix vom Typ 'scipy.sparse.lil_matrix'
        :return: Menge der 2-Tupel mit gesuchten Indizes
        """
        s_coo = s.tocoo()
        return set(zip(s_coo.row, s_coo.col))'''

    def fp_approx(self, theta_ind, top_ords, delta_p_act, nu_min):
        theta = self.global_phase[-1]
        xk = scipy.sparse.lil_matrix((self.I, self.m))
        xfp = scipy.sparse.lil_matrix((self.I, self.m))
        x_sum = np.zeros(self.m)
        x_sum_fix = np.zeros(self.m)
        nu_sum_act = []
        coms = defaultdict(list)
        init_coms = defaultdict(list)
        coms_lens = np.zeros(self.n)  # 'coms_lens[v_ind]' entspricht der Länge von 'init_coms[v_ind]'
        self.flow_vol.append(scipy.sparse.lil_matrix((self.n, 1)))
        print(theta_ind, theta)

        def calc_flow_by_bounds(xk_old, subset=None):
            """
            :param xk_old:
            :param subset: Teilmenge von 'self.I' x 'self.V': falls angegeben, so enthält 'subset' alle (i, v)-Paare, für die 'bounds' relaxiert wurden. Dann werden nur für diese
                           Paare neue Flusswerte berechnet (passend zu den neuen 'bounds')
            :return:
            """

            xk = scipy.sparse.lil_matrix((self.I, self.m))
            x_sum = np.zeros(self.m)

            if subset is None:
                calc_v = coms.copy()
            else:
                calc_v = defaultdict(list)
                for (i, v_ind) in subset:
                    calc_v[v_ind].append(i)
                for v_ind in coms:
                    for i in coms[v_ind]:
                        if v_ind not in calc_v or i not in calc_v[v_ind]:
                            for e_ind in delta_p_act[i][v_ind]:
                                xk[i, e_ind] = xk_old[i, e_ind]
                                x_sum[e_ind] += xk[i, e_ind]

            for v_ind in calc_v:
                # enthält für jedes Gut in 'v_ind' Menge des bereits einer Kante zugewiesenen Flussvolumens
                flow_sent = {}
                for i in calc_v[v_ind]:
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
                        xk[i, e_ind] = xk_old[i, e_ind]
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
                    a_e_min = 1
                    for e_ind in delta_p_act[i][v_ind]:
                        w_ind = self.V.index(self.E[e_ind][1])
                        a_e = self.g(e_ind, x_sum[e_ind] + x_sum_fix[e_ind]) / self.nu[e_ind] + a[i, w_ind]
                        a_e_tol = self.bd_tol / self.nu[e_ind] + self.bd_tol / a_e_min
                        if (i, v_ind) not in sparse_items:
                            if abs(a_e) > a_e_tol + self.eps:
                                a[i, v_ind] = a_e
                            sparse_items.append((i, v_ind))
                            argmin[i][v_ind] = [e_ind]
                            a_e_min = self.nu[e_ind]
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
                                a_e_min = self.nu[e_ind]
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
                                a_e_min = self.nu[e_ind]
                                argminiv = copy.deepcopy(argmin[i][v_ind])
                                for e_arg in argminiv:
                                    a_e_arg = self.g(e_arg, x_sum[e_arg] + x_sum_fix[e_arg]) / self.nu[e_arg] + a[i, self.V.index(self.E[e_arg][1])]
                                    a_e_tol = self.bd_tol / self.nu[e_arg] + self.bd_tol / a_e_min
                                    if a_e_arg > a[i, v_ind] + a_e_tol:
                                        argmin[i][v_ind].remove(e_arg)
                                        if e_arg not in bounds_f[v_ind][i] and (i, e_arg) not in new_f:
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
            relaxed = []
            for v_ind in coms:
                for i in coms[v_ind]:
                    arg = list(set(delta_p_act[i][v_ind]) - set(argmin[i][v_ind]))
                    arg_nf = list(set(arg) - set(bounds_f[v_ind][i]))
                    if not arg_nf:
                        # In diesem Fall sind alle aktiven nichtminimalen Kanten bereits fixiert. Damit können keine Fortschritte mehr gemacht werden -> relaxiere bds
                        relaxed.append((i, v_ind))
                        for e_ind in arg:
                            if xk[i, e_ind] < 2 * self.bd_tol / nu_min[i, v_ind]:
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
                        relaxed.append((i, v_ind))
                        for e_ind in argmin[i][v_ind]:
                            if xk[i, e_ind] > self.b[-1][i, v_ind] - 2 * self.bd_tol / nu_min[i, v_ind]:
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
                            bounds[v_ind][i][e_ind] = (xk[i, e_ind], np.min([x2 + xk[i, e_ind], self.b[-1][i, v_ind]]))
                            b_nf[i, v_ind] += xk[i, e_ind]
            return relaxed

        def fix_nodes():
            def search_path(v_ind, future_w, i):
                visit = 0
                len_visited = 1
                visited = [v_ind]
                while visit < len_visited:
                    node = visited[visit]
                    for e_ind in delta_p_act[i][node]:
                        successor = self.V.index(self.E[e_ind][1])
                        if successor not in visited:
                            if successor in future_w:
                                return True
                            visited.append(successor)
                            len_visited += 1
                    visit += 1
                return False

            def fix_check(v_ind):
                delta_p_act_nf_pos_v = {}
                for i in coms[v_ind]:
                    delta_p_act_nf_pos_v[i] = [e_ind for e_ind in delta_p_act_nf[i][v_ind] if xk[i, e_ind] > 2 * self.bd_tol / nu_min[i, v_ind]]
                # Prüfe, ob alle aktiven Kanten mit positivem Flusswert den gleichen a-Wert liefern
                if np.all([v_ind not in argmin_nf[i] or np.all([e_ind in argmin_nf[i][v_ind] for e_ind in delta_p_act_nf_pos_v[i]]) for i in coms[v_ind]]):
                    # alle Kanten, für die ein Flusswert >0 fixiert wurde, müssen ebenfalls in 'argmin[i][v_ind]' liegen
                    if np.any([[e_ind for e_ind in bounds_f[v_ind][i] if xk[i, e_ind] > 2 * self.bd_tol / nu_min[i, v_ind] and e_ind not in argmin[i][v_ind]] for i in coms[v_ind]]):
                        # sonst -> Abbruch
                        return False
                    for i in coms[v_ind]:
                        to_ind = top_ords[i].index(v_ind)
                        future_w = [w_ind for w_ind in top_ords[i][1:to_ind] if w_ind in coms]
                        if search_path(v_ind, future_w, i):
                            # ist ein Knoten der zwischen 'v_ind' und 'ti' liegt noch nicht approximiert -> Abbruch
                            return False
                    # sonst: 'v_ind' kann approximiert werden
                    return True
                return False

            ci = 0
            coms_keys = list(coms)
            cn = len(coms_keys)
            while ci < cn:
                v_ind = coms_keys[ci]
                if fix_check(v_ind):
                    fc_err = 0
                    fc = True
                    for i in coms[v_ind]:
                        for e_ind in delta_p_act[i][v_ind]:
                            if xk[i, e_ind] < 2 * self.bd_tol / nu_min[i, v_ind]:
                                xfp[i, e_ind] = 0
                                fc_err += xk[i, e_ind]
                            elif xk[i, e_ind] > self.b[-1][i, v_ind] - 2 * self.bd_tol / nu_min[i, v_ind]:
                                xfp[i, e_ind] = self.b[-1][i, v_ind]
                                x_sum_fix[e_ind] += xfp[i, e_ind]
                                fc = False
                            else:
                                xfp[i, e_ind] = xk[i, e_ind]
                                x_sum_fix[e_ind] += xfp[i, e_ind]
                            x_sum[e_ind] -= xk[i, e_ind]
                            xk[i, e_ind] = 0
                        if fc_err and fc:
                            # Stelle Flusserhaltung wieder her, falls diese durch 0 - Setzung obiger Flusswerte verletzt
                            for e_ind in delta_p_act[i][v_ind]:
                                if xfp[i, e_ind]:
                                    add_flow = fc_err * xfp[i, e_ind] / (self.b[-1][i, v_ind] - fc_err)
                                    xfp[i, e_ind] += add_flow
                                    x_sum_fix[e_ind] += add_flow
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
                if self.b[-1][i, v_ind] > 0:
                    coms_lens[v_ind] += 1
                    self.flow_vol[-1][v_ind, 0] += self.b[-1][i, v_ind]
                    if len(delta_p_act[i][v_ind]) == 1:
                        e_ind = delta_p_act[i][v_ind][0]
                        xfp[i, e_ind] = self.b[-1][i, v_ind]
                        x_sum_fix[e_ind] += xfp[i, e_ind]
                        init_coms[v_ind].append(i)
                        continue
                    coms[v_ind].append(i)
                    init_coms[v_ind].append(i)
                    nu_sum_act[i][v_ind] = 0
                    for e_ind in delta_p_act[i][v_ind]:
                        nu_sum_act[i][v_ind] += self.nu[e_ind]
                    for e_ind in delta_p_act[i][v_ind]:
                        xk[i, e_ind] = self.b[-1][i, v_ind] * self.nu[e_ind] / nu_sum_act[i][v_ind]
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
        b_nf = copy.deepcopy(self.b[-1])
        fp_comp = []

        while True:
            # nur im ersten Schritt leer, danach enthält 'fp_comp' Tupel der Form '(i, v_ind)', was impliziert, dass die Aufteilung von Gut 'i' im Knoten 'v_ind'
            # vollständig fixiert wurde (d.h. die Boundsintervalle aller aktiven ausgehenden Kanten sind hinreichend klein)
            for (i, v_ind) in fp_comp:
                fc_err = 0
                fc = True
                for e_ind in bounds[v_ind][i]:
                    if xk[i, e_ind] < self.bd_tol * 1 / coms_lens[v_ind]:  # ???
                        xfp[i, e_ind] = 0
                        fc_err += xk[i, e_ind]
                    elif xk[i, e_ind] > self.b[-1][i, v_ind] - self.bd_tol * 1 / coms_lens[v_ind]:
                        xfp[i, e_ind] = self.b[-1][i, v_ind]
                        x_sum_fix[e_ind] += xfp[i, e_ind]
                        fc = False
                    else:
                        xfp[i, e_ind] = xk[i, e_ind]
                        x_sum_fix[e_ind] += xfp[i, e_ind]
                    x_sum[e_ind] -= xk[i, e_ind]
                if fc_err and fc:
                    # Stelle Flusserhaltung wieder her, falls diese durch 0 - Setzung obiger Flusswerte verletzt
                    for e_ind in bounds[v_ind][i]:
                        if xfp[i, e_ind]:
                            add_flow = fc_err * xfp[i, e_ind] / (self.b[-1][i, v_ind] - fc_err)
                            xfp[i, e_ind] += add_flow
                            x_sum_fix[e_ind] += add_flow
                del bounds[v_ind][i]
                del bounds_f[v_ind][i]
                coms[v_ind].remove(i)
                if not coms[v_ind]:
                    del coms[v_ind]

            if fp_comp:  # nur im 1. Schritt nicht erfüllt. In diesem wird x0 verwendet
                if not coms:  # Abbruchbedingung: keine Güter mehr übrig, die noch auf Kanten aufgeteilt werden müssen
                    # verbessere potentielle Approximationsfehler in 'a' durch erneuten Aufruf von 'calc_a'. Dies ist möglich, da eventuell in vorheriger 'for' - Schleife
                    # Approximationsfehler in 'xk' verbessert wurden
                    a, a_min2, argmin, argmin_nf = calc_a()
                    return xfp, x_sum_fix, a, argmin, init_coms
                xk, x_sum = calc_flow_by_bounds(xk_old=xk)

            a, a_min2, argmin, argmin_nf = calc_a()

            relaxed = relax_bounds()

            if relaxed:
                xk, x_sum = calc_flow_by_bounds(xk_old=xk, subset=relaxed)

            fix_nodes()

            if not coms:  # Abbruchbedingung: keine Güter mehr übrig, die noch auf Kanten aufgeteilt werden müssen
                # verbessere potentielle Approximationsfehler in 'a' durch erneuten Aufruf von 'calc_a'. Dies ist möglich, da eventuell in vorheriger 'for' - Schleife
                # Approximationsfehler in 'xk' verbessert wurden
                a, a_min2, argmin, argmin_nf = calc_a()
                return xfp, x_sum_fix, a, argmin, init_coms

            while True:
                forced_zeros = self.gamma(argmin, delta_p_act_nf, coms)
                fp_comp = []
                new_f = []
                for v_ind in coms:
                    for i in coms[v_ind]:
                        dpanf_start = len(delta_p_act_nf[i][v_ind])
                        fz_nf = list(set(forced_zeros[i][v_ind]) - set(bounds_f[v_ind][i]))
                        for e_ind in fz_nf:
                            if e_ind in bounds[v_ind][i]:
                                bounds[v_ind][i][e_ind] = (bounds[v_ind][i][e_ind][0], xk[i, e_ind])
                            else:
                                # initiale untere Schranke: so, dass durch Verringerung von xk[i, e_ind] der für 'v_ind' minimale a_i - Wert erreicht werden kann
                                bounds[v_ind][i][e_ind] = (np.max([0, xk[i, e_ind] + (a[i, v_ind] - a[i, self.V.index(self.E[e_ind][1])]) * self.nu[e_ind] - self.g(e_ind, x_sum[e_ind] + x_sum_fix[e_ind])]), xk[i, e_ind])
                            bd_dif = abs(bounds[v_ind][i][e_ind][1] - bounds[v_ind][i][e_ind][0])
                            if bd_dif < self.bd_tol * 1 / coms_lens[v_ind]:
                                if bd_dif == 0 and xk[i, e_ind] > self.eps:
                                    # bd wird relaxiert
                                    w_ind = self.V.index(self.E[e_ind][1])
                                    if self.q[e_ind] < self.eps:
                                        # Falls q = 0 -> g_e nicht-negativ -> kleinster a-Wert kann nicht erreicht werden -> wähle kleinstmögliche untere Schranke: 0
                                        bounds[v_ind][i][e_ind] = (0, xk[i, e_ind])
                                    else:
                                        # relaxiere untere Schranke, sodass der minimale a-Wert erreicht werden kann
                                        bounds[v_ind][i][e_ind] = (np.max([0, xk[i, e_ind] - (x_sum[e_ind] + x_sum_fix[e_ind] - (1 + a[i, v_ind] - a[i, w_ind]) * self.nu[e_ind])]), xk[i, e_ind])
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
                                x2 = (a_min2[i, v_ind] - a[i, w_ind]) * self.nu[e_ind] - self.g(e_ind, xe)
                                if self.q[e_ind] < self.eps and xe < self.nu[e_ind] - self.eps:
                                    x2 += self.nu[e_ind] - xe
                                bounds[v_ind][i][e_ind] = (xk[i, e_ind], np.min([x2 + xk[i, e_ind], self.b[-1][i, v_ind]]))
                            bd_dif = abs(bounds[v_ind][i][e_ind][1] - bounds[v_ind][i][e_ind][0])
                            if bd_dif < self.bd_tol * 1 / coms_lens[v_ind]:
                                if bd_dif == 0 and xk[i, e_ind] < self.b[-1][i, v_ind] - self.eps:
                                    # bd wird relaxiert
                                    w_ind = self.V.index(self.E[e_ind][1])
                                    xe = x_sum[e_ind] + x_sum_fix[e_ind]
                                    x2 = (a_min2[i, v_ind] - a[i, w_ind]) * self.nu[e_ind] - self.g(e_ind, xe)
                                    if self.q[e_ind] < self.eps and xe < self.nu[e_ind] - self.eps:
                                        x2 += self.nu[e_ind] - xe
                                    # relaxierte obere Schranke mit x2: Dies ist der Flusswert, sodass der zweitkleinste a-Wert (= 'a_min2[i, v_ind]') auch für Kante 'e_ind'
                                    # erreicht werden kann.
                                    bounds[v_ind][i][e_ind] = (xk[i, e_ind], np.min([x2 + xk[i, e_ind], self.b[-1][i, v_ind]]))
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
                    if b_nf[i, v_ind] < 2 * self.bd_tol / nu_min[i, v_ind]:
                        if np.all([e in argmin[i][v_ind] for e in bounds_f[v_ind][i] if xk[i, e] > 2 * self.bd_tol / nu_min[i, v_ind]]):
                            if (i, v_ind) not in fp_comp:  # ???
                                fp_comp.append((i, v_ind))
                        else:
                            for e in bounds_f[v_ind][i]:
                                delta_p_act_nf[i][v_ind].append(e)
                                if e in argmin[i][v_ind]:
                                    if e not in argmin_nf[i][v_ind]:
                                        argmin_nf[i][v_ind].append(e)
                                    w_ind = self.V.index(self.E[e][1])
                                    xe = x_sum[e] + x_sum_fix[e]
                                    x2 = (a_min2[i, v_ind] - a[i, w_ind]) * self.nu[e] - self.g(e, xe)
                                    if self.q[e] < self.eps and xe < self.nu[e] - self.eps:
                                        x2 += self.nu[e] - xe
                                    bounds[v_ind][i][e] = (xk[i, e], np.min([x2 + xk[i, e], self.b[-1][i, v_ind]]))
                                else:
                                    if self.q[e] > self.eps:
                                        bounds[v_ind][i][e] = (np.max([0, xk[i, e] + (a[i, v_ind] - a[i, self.V.index(self.E[e][1])]) * self.nu[e] - self.g(e, x_sum[e] + x_sum_fix[e])]), xk[i, e])
                                    else:
                                        bounds[v_ind][i][e] = (0, xk[i, e])
                            bounds_f[v_ind][i] = []
                            b_nf[i, v_ind] = self.b[-1][i, v_ind]

                if fp_comp:
                    # es existieren Knoten in 'fp_comp', welche fixiert werden können
                    break

                relaxed = relax_bounds()

                if relaxed:
                    xk, x_sum = calc_flow_by_bounds(xk_old=xk, subset=relaxed)

                fix_nodes()

                if not coms: # Abbruchbedingung: keine Güter mehr übrig, die noch auf Kanten aufgeteilt werden müssen
                    # verbessere potentielle Approximationsfehler in 'a' durch erneuten Aufruf von 'calc_a'. Dies ist möglich, da eventuell in vorheriger 'for' - Schleife
                    # Approximationsfehler in 'xk' verbessert wurden
                    a, a_min2, argmin, argmin_nf = calc_a()
                    return xfp, x_sum_fix, a, argmin, init_coms

    def main(self, theta=None, x_total=None, x_sum=None, init_coms=None, dp_act_old=None, max_dif=None, delta_p_act=None, delta_p_inact=None, top_ords=None, nu_min=None):
        """
        Hauptteil: Konstruiert schrittweise einen kontinuierlichen IDE-Fluss zu den Eingabedaten aus 'Math.cont_data'.
        :return: 0
        """
        if theta is None:
            theta = 0
            theta_ind = -1
            max_dif = self.bd_tol
            delta_p_act = []
            delta_p_inact = []
            top_ords = []
            nu_min = scipy.sparse.lil_matrix((self.I, self.n))
        else:
            theta_ind = len(self.global_phase) - 2
        # Obergrenze für theta
        T = 1000
        if theta == 0:
            # berechne 'delta_p_act', delta_p_inact', 'top_ords' und 'self.b' für 'theta' = 0
            for i in range(self.I):
                delta_p_act.append({})
                delta_p_inact.append({})
                for v_ind in range(self.n):
                    delta_p_act[i][v_ind] = self.get_outgoing_active_edges(i, v_ind)
                    delta_p_inact[i][v_ind] = [e for e in self.delta_p[v_ind] if e not in delta_p_act[i][v_ind]]
                    if delta_p_act[i][v_ind]:
                        nu_min[i, v_ind] = np.min([self.nu[e] for e in delta_p_act[i][v_ind]])
                    self.b[-1][i, v_ind] = self.calc_b(i, v_ind, theta)
                # berechne top. Sortierung für aktive Kanten von Gut 'i'
                top_ords.append(self.topologicalSort(i, delta_p_act[i]))
        skip_count = 0
        while theta < T:
            # in der Zukunft liegende Zeitpunkte aus der Liste 'self.u_start'
            start_points = [t for t in self.u_start if t > theta]
            theta_ind += 1
            # Liste aller Kanten, deren Warteschlange in der aktuellen Phase 0 wird, und deren f^- -Werte auf 0
            # gesetzt werden müssen
            fm_to_zero = {}
            # 'next_phase[-1]' bestimmt am Ende der Schleife, wie lange aktuelle Phase andauert
            if start_points:
                next_sp = np.min(start_points)
                next_phase = [next_sp - theta]
            else:
                next_phase = [T]
                next_sp = np.Inf
            skip_calcs = False
            if theta > 0:
                old_b = copy.deepcopy(self.b[-1])
                self.b.append(scipy.sparse.lil_matrix((self.I, self.n)))
                for i in range(self.I):
                    for v_ind in range(self.n):
                        self.b[-1][i, v_ind] = self.calc_b(i, v_ind, theta)
                # prüfe, ob 'self.b' unverändert -> es besteht Möglichkeit, dass vorige Flusswerte immer noch korrekt
                if not (self.b[-1] - old_b).count_nonzero():
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
                            a_e_min = 1
                            for e_ind in delta_p_act[i][v_ind]:
                                w_ind = self.V.index(self.E[e_ind][1])
                                a_e = self.g(e_ind, x_sum[e_ind]) / self.nu[e_ind] + a[i, w_ind]
                                a_e_tol = self.bd_tol / self.nu[e_ind] + self.bd_tol / a_e_min
                                if (i, v_ind) not in sparse_items:
                                    if abs(a_e) > a_e_tol + self.eps:
                                        a[i, v_ind] = a_e
                                    sparse_items.append((i, v_ind))
                                    argmin[i][v_ind] = [e_ind]
                                    a_e_min = self.nu[e_ind]
                                elif a_e < a[i, v_ind] + a_e_tol:
                                    if a_e < a[i, v_ind] - a_e_tol:
                                        if abs(a_e) > a_e_tol + self.eps:
                                            a[i, v_ind] = a_e
                                        else:
                                            a[i, v_ind] = 0
                                        argmin[i][v_ind] = [e_ind]
                                        a_e_min = self.nu[e_ind]
                                    elif a_e < a[i, v_ind]:
                                        if abs(a_e) > a_e_tol + self.eps:
                                            a[i, v_ind] = a_e
                                        else:
                                            a[i, v_ind] = 0
                                        a_e_min = self.nu[e_ind]
                                        argminiv = copy.deepcopy(argmin[i][v_ind])
                                        for e_arg in argminiv:
                                            a_e_arg = self.g(e_arg, x_sum[e_arg]) / self.nu[e_arg] + a[i, self.V.index(self.E[e_arg][1])]
                                            a_e_tol = self.bd_tol / self.nu[e_arg] + self.bd_tol / a_e_min
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
                            dpact_pos_flow[v_ind][i] = [e_ind for e_ind in delta_p_act[i][v_ind] if x_total[i, e_ind] > self.eps]
                            if np.any([e_ind not in dpact_pos_flow[v_ind][i] for e_ind in dp_act_old[i][v_ind] if x_total[i, e_ind] > self.eps]):
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
                        x_total, x_sum, a, argmin, init_coms = self.fp_approx(theta_ind, top_ords, delta_p_act, nu_min)
                else:
                    x_total, x_sum, a, argmin, init_coms = self.fp_approx(theta_ind, top_ords, delta_p_act, nu_min)
            else:
                x_total, x_sum, a, argmin, init_coms = self.fp_approx(theta_ind, top_ords, delta_p_act, nu_min)

            delta_c = []
            for v_ind in range(self.n):
                for e_ind in self.delta_p[v_ind]:
                    delta_c.append(self.change_of_cost(e_ind, x_sum[e_ind]))

            for ti in range(self.I):
                for v_ind in top_ords[ti]:
                    if self.delta_p[v_ind]:
                        dp_min_nu = np.min([self.nu[e] for e in self.delta_p[v_ind]])
                    # betrachte aktive Kanten
                    for e_ind in delta_p_act[ti][v_ind]:
                        if not skip_calcs:  # 'fp' bleibt gleich
                            if theta == 0:
                                if x_total[ti, e_ind] > self.eps:
                                    self.fp[ti][e_ind][0] = (0, x_total[ti, e_ind])
                                    self.fp_ind[ti][e_ind].append(0)
                            # falls sich f^+ -Wert in dieser Phase ändert, aktualisiere 'self.fp'
                            # Vergleichen hier zwei unterschiedliche Zeitpunkte
                            elif abs(self.fp[ti][e_ind][-1][1] - x_total[ti, e_ind]) > 2 * (self.bd_tol / self.nu[e_ind] + self.bd_tol / dp_min_nu):
                                self.fp[ti][e_ind].append((theta, x_total[ti, e_ind]))
                                self.fp_ind[ti][e_ind].append(theta_ind)
                                if self.q_ind[e_ind][-1] != theta_ind and (self.q[e_ind] > self.eps or self.q_global[e_ind][-1] > self.eps or x_sum[e_ind] > self.nu[e_ind] + self.eps):
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
                        # Vergleichen hier zwei unterschiedliche Zeitpunkte
                        if abs(self.fm[ti][e_ind][fm_ind][1] - outflow) > 2 * (self.bd_tol / self.nu[e_ind] + self.bd_tol / dp_min_nu):
                            if abs(self.fm[ti][e_ind][fm_ind][0] - outflow_time) < self.eps:
                                # zuvor berechneter 'outflow' kann sich durch Änderung der Flusswerte verzögern
                                del self.fm[ti][e_ind][fm_ind]
                            else:
                                self.fm[ti][e_ind].append((outflow_time, outflow))
            for ti in range(self.I):
                # überspringen jeweilige Senke, da diese hier uninteressant
                for v_ind in top_ords[ti][1:]:
                    for e_ind in delta_p_act[ti][v_ind]:
                        last_fm_ind = self.last_fm_change(ti, e_ind, theta)
                        max_err = self.bd_tol / self.nu[e_ind] + self.bd_tol / self.nu[argmin[ti][v_ind][0]]
                        # bestimme, ob sich vor dem Ende der aktuellen Phase ein f^- -Wert ändert -> verkürze Phase
                        if len(self.fm[ti][e_ind]) > last_fm_ind + 1 and self.eps < self.fm[ti][e_ind][last_fm_ind + 1][0] - theta < next_phase[0] + max_dif + max_err * next_phase[0] / self.nu[e_ind]:
                            next_fm = self.fm[ti][e_ind][last_fm_ind + 1][0] - theta
                            if next_fm < next_phase[0]:
                                nph_rev = copy.deepcopy(next_phase)
                                nph_rev.reverse()
                                for p in nph_rev:
                                    if p > next_fm + max_dif + max_err * next_fm / self.nu[e_ind]:
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
                                    if e_ind == 29:
                                        print()
                                    self.q_ind[e_ind].append(theta_ind)
                                if self.q[e_ind] > self.eps:
                                    outflow_time = theta + self.q[e_ind]/self.nu[e_ind] + self.r[e_ind]
                                else:
                                    outflow_time = theta + self.r[e_ind]
                                fm_ind = self.last_fm_change(ti, e_ind, outflow_time)
                                # setze 'self.fm'- Wert auf 0, falls dies noch nicht geschehen ist
                                if abs(self.fm[ti][e_ind][fm_ind][1]) > 0:
                                    self.fm[ti][e_ind].append((outflow_time, 0))

                        # bestimme, ob sich vor dem Ende der aktuellen Phase ein f^- -Wert ändert -> verkürze Phase
                        last_fm_ind = self.last_fm_change(ti, e_ind, theta)
                        max_err = self.bd_tol / self.nu[e_ind] + self.bd_tol / self.nu[argmin[ti][v_ind][0]]
                        if len(self.fm[ti][e_ind]) > last_fm_ind + 1 and self.eps < self.fm[ti][e_ind][last_fm_ind + 1][0] - theta < next_phase[0] + max_dif + max_err * next_phase[0] / self.nu[e_ind]:
                            next_fm = self.fm[ti][e_ind][last_fm_ind + 1][0] - theta
                            if next_fm < next_phase[0]:
                                nph_rev = copy.deepcopy(next_phase)
                                nph_rev.reverse()
                                for p in nph_rev:
                                    if p > next_fm + max_dif + max_err * next_fm / self.nu[e_ind]:
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
                            max_err = self.bd_tol / self.nu[e_ind] + self.bd_tol / self.nu[argmin[ti][v_ind][0]]
                            dec_rate = self.g(active_ind, x_sum[active_ind]) / self.nu[active_ind] + a[ti, act_ind] - self.g(e_ind, x_sum[e_ind]) / self.nu[e_ind] - a[ti, tar_ind]
                            if dec_rate < 2 * max_err / nu_min[ti, v_ind] + self.eps:
                                # Kante wird sicher nicht aktiv
                                continue
                            np_tol = max_dif + (max_err * next_phase[0] / self.nu[active_ind] + max_err * next_phase[0] / self.nu[e_ind]) / dec_rate

                            if self.labels[ti][-1][tar_ind] + self.q[e_ind]/self.nu[e_ind] + self.r[e_ind] + \
                                    (delta_c[e_ind] + a[ti, tar_ind]) * (next_phase[0] + np_tol) < \
                                    self.labels[ti][-1][act_ind] + self.q[active_ind]/self.nu[active_ind] + \
                                    self.r[active_ind] + (delta_c[active_ind] + a[ti, act_ind]) * (next_phase[0] + np_tol) and \
                                    abs(delta_c[e_ind] + a[ti, tar_ind] - delta_c[active_ind] - a[ti, act_ind]) \
                                    > self.eps:
                                time_ub = (self.labels[ti][-1][act_ind] + self.q[active_ind]/self.nu[active_ind] + self.r[active_ind] - \
                                           self.labels[ti][-1][tar_ind] - self.q[e_ind]/self.nu[e_ind] - self.r[e_ind]) \
                                          / (delta_c[e_ind] + a[ti, tar_ind] - delta_c[active_ind] - a[ti, act_ind])

                                if time_ub < next_phase[0]:
                                    nph_rev = copy.deepcopy(next_phase)
                                    nph_rev.reverse()
                                    np_tol = max_dif + (max_err * time_ub / self.nu[active_ind] + max_err * time_ub / self.nu[e_ind]) / dec_rate
                                    for p in nph_rev:
                                        if p > time_ub + np_tol:
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
                if self.delta_p[v_ind]:
                    max_err = 2 * self.bd_tol / np.min([self.nu[e] for e in self.delta_p[v_ind]])
                for e_ind in self.delta_p[v_ind]:
                    change_of_q = delta_c[e_ind] * self.nu[e_ind]
                    # Falls die Warteschlange von 'e' unter aktuellem Fluss abgebaut wird, bestimme Zeitpunkt, zu dem diese vollständig abgebaut ist (bei gleich bleibendem
                    # Fluss)
                    if change_of_q < -self.eps:
                        # 'phase_length': Dauer bis Warteschlangenlänge gleich 0
                        phase_length = - self.q[e_ind] / change_of_q
                        q_tol = max_dif + max_err * next_phase[0] / - change_of_q
                        if phase_length < next_phase[0] + q_tol:
                            if phase_length < next_phase[0]:
                                nph_rev = copy.deepcopy(next_phase)
                                nph_rev.reverse()
                                for p in nph_rev:
                                    q_tol = max_dif + max_err * phase_length / - change_of_q
                                    if p > phase_length + q_tol:
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
                    if np.any([theta_ind in self.fp_ind[i][e_ind] for i in range(self.I)]) and (self.q[e_ind] > self.eps or self.q_global[e_ind][-1] > self.eps or x_sum[e_ind] > self.nu[e_ind] + self.eps):
                        self.q_global[e_ind].append(self.q[e_ind])

            alpha = (next_phase[-1] + next_phase[0]) / 2
            if next_sp - theta - alpha < alpha - next_phase[0] + self.eps:
                # Ausnahme bei der Schrittweitenwahl: Zeitpunkte mit Änderungen der externen Einflussraten dürfen nicht durch einen früheren Zeitpunkt approximiert werden, da sie
                # sonst übersprungen werden
                alpha = next_sp - theta
            theta += alpha
            if alpha != T:
                # aktualisiere Warteschlangenlängen und Kosten
                new_q = []
                self.c.append(np.zeros(self.m))
                for v_ind in range(self.n):
                    if self.delta_p[v_ind]:
                        q_len_tol = np.max([self.bd_tol, 2 * self.bd_tol / np.min([self.nu[e] for e in self.delta_p[v_ind]]) * alpha])
                    for e_ind in self.delta_p[v_ind]:
                        next_q_len = self.q[e_ind] + self.g(e_ind, x_sum[e_ind]) * alpha
                        if next_q_len < q_len_tol:
                            next_q_len = 0
                            if self.q[e_ind] > self.eps:
                                self.q_global[e_ind].append(0)
                                self.q_ind[e_ind].append(theta_ind + 1)
                        new_q.append(next_q_len)
                        self.c[-1][e_ind] = new_q[e_ind] / self.nu[e_ind] + self.r[e_ind]

                # speichere aktuelle Warteschlangenlängen
                self.q = new_q

                # speichere Phase
                self.global_phase.append(theta)
                self.E_active = scipy.sparse.lil_matrix((self.I, self.m))
                qi = []
                for e_ind in range(self.m):
                    qi.append({})
                for i in range(self.I):
                    self.labels[i].append(np.zeros(self.n))
                    for v_ind in range(self.n):
                        self.labels[i][-1][v_ind] = self.labels[i][-2][v_ind] + alpha * a[i, v_ind]
                    # etwas effizienter als 'v_ind in range(self.n)', da t_i übersprungen wird. Aber: 'top_ords[i]' hier keine aktuelle top. Sortierung mehr
                    for v_ind in top_ords[i][1:]:
                        for edge in self.delta_p[v_ind]:
                            w_ind = self.V.index(self.E[edge][1])
                            label_dif = self.labels[i][-1][v_ind] - self.labels[i][-1][w_ind] - self.c[-1][edge]
                            max_err = self.bd_tol / self.nu[edge] + self.bd_tol / self.nu[argmin[i][v_ind][0]]
                            if label_dif > -(max_dif + 2 * max_err / nu_min[i, v_ind] * alpha):
                                # 'label_difs[i, edge]' > 0: geschieht nur aufgrund von Approximationsfehlern, für Kanten die aktiv sein sollten
                                self.E_active[i, edge] = 1

                dp_act_old = copy.deepcopy(delta_p_act)
                for i in range(self.I):
                    for v_ind in range(self.n):
                        delta_p_act[i][v_ind] = self.get_outgoing_active_edges(i, v_ind)
                        delta_p_inact[i][v_ind] = [e for e in self.delta_p[v_ind] if e not in delta_p_act[i][v_ind]]
                    # berechne top. Sortierung zur aktualisierten Menge aktiver Kanten für Gut 'i'
                    top_ords[i] = self.topologicalSort(i, delta_p_act[i])

                for i in range(self.I):
                    for v_ind in range(self.n):
                        for e_ind in delta_p_act[i][v_ind]:
                            w_ind = self.V.index(self.E[e_ind][1])
                            label_dif = abs(self.labels[i][-1][v_ind] - self.labels[i][-1][w_ind] - self.c[-1][e_ind])
                            if label_dif > max_dif:
                                max_dif = label_dif

                for i in range(self.I):
                    for v_ind in top_ords[i][1:]:
                        ell_sum = 0
                        for e_ind in delta_p_act[i][v_ind]:
                            w_ind = self.V.index(self.E[e_ind][1])
                            ell_sum += self.labels[i][-1][w_ind] + self.c[-1][e_ind]
                        self.labels[i][-1][v_ind] = ell_sum / len(delta_p_act[i][v_ind])
                        if self.labels[i][-1][v_ind] < max_dif:
                            self.labels[i][-1][v_ind] = 0
                        for e_ind in delta_p_act[i][v_ind]:
                            w_ind = self.V.index(self.E[e_ind][1])
                            # 'qi[e_ind][i]': Warteschlangenlänge von Kante 'e_ind', passend zu den labels von Gut 'i'
                            qi[e_ind][i] = (self.labels[i][-1][v_ind] - self.labels[i][-1][w_ind] - self.r[e_ind]) * self.nu[e_ind]
                            if qi[e_ind][i] < max_dif * self.nu[e_ind] + self.eps:
                                qi[e_ind][i] = 0

                for e_ind in range(self.m):
                    len_qie = len(qi[e_ind])
                    if len_qie > 0:
                        mean_qi = 0
                        for i in qi[e_ind]:
                            mean_qi += qi[e_ind][i]
                        mean_qi /= len_qie
                        self.q[e_ind] = mean_qi
                        self.c[-1][e_ind] = self.q[e_ind] / self.nu[e_ind] + self.r[e_ind]

                for p in fm_to_zero:
                    for e_ind in fm_to_zero[p]:
                        for i in range(self.I):
                            if self.fm[i][e_ind][-1][1] > self.eps:
                                self.fm[i][e_ind].append((theta + self.r[e_ind], 0))

        print("phases", len(self.global_phase))
        print("davon geskipt:", skip_count)

        # erzeuge Ausgabe
        OutputTableMulti(self.G, self.V, self.E, self.I, self.r, self.nu, self.fp, self.fp_ind, self.fm, self.q_global, self.q_ind, self.global_phase, self.labels, self.flow_vol, self.bd_tol)
        return 0
