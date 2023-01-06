import numpy as np
from OutputTableMulti import OutputTableMulti
from scipy.sparse.linalg import spsolve
from scipy.optimize import linprog
import scipy
import copy
import itertools
from collections import defaultdict
import time


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
        self.c = []  # Kantenkosten für alle Kanten zu Beginn jeder Phase (Reisezeit + Warteschlangenlänge / Kapazität)
        self.G = G
        self.u = u
        self.I = len(self.u)
        self.labels = []  # Knotenlabels
        self.V = list(G.keys())  # Liste der Knoten
        self.n = len(self.V)  # Anzahl Knoten
        self.b = []
        self.items = G.items()
        self.keys = G.keys()
        # self.eps = 10**(-12)  # Für Rundungsfehler
        self.eps = 10**(-6)
        self.bd_tol = 10**(-6)
        self.flow_vol = []  # merke Flusswerte in den einzelnen Knoten für OutputTable
        self.delta_p = []

        for delta in self.items:
            for w in list(delta[1].keys()):
                self.E.append((delta[0], w))  # speichere alle Kanten in E
                self.r.append(list(self.items)[list(self.keys).index(delta[0])][1][w][1])  # speichere Reisezeiten in r
                self.nu.append(list(self.items)[list(self.keys).index(delta[0])][1][w][0])  # Kapazitäten in nu
        self.m = len(self.E)  # Anz. Kanten
        self.global_phase = [0]  # speichere Startzeitpunkte der "globalen" Phasen, also alle Zeitpunkte, zu denen sich mindestens eine lokale Phase ändert
        self.q_global = [self.m * [0]]  # speichere für jede globale Phase alle momentanen Warteschlangenlängen zu Beginn der Phase

        self.fp = [[] for _ in range(self.I)]  # f^+    [[],[]]
        self.fp_ind = [[] for _ in range(self.I)]  # Liste der Indizes der Phasen
        self.fm = [[] for _ in range(self.I)]  # f^-
        for i in range(self.I):
            for k in range(self.m):
                self.fp[i].append([(0, 0)])
                self.fp_ind[i].append([])
                self.fm[i].append([(0, 0)])

        self.c.append(np.zeros(self.m))
        for e_ind in range(self.m):  # Initialisierung "self.c" (Kosten)
            self.c[0][e_ind] = self.r[e_ind]

        # Zeitpunkte, zu denen sich der Zufluss in mindestens einem Quellknoten ändert
        self.u_start = set()
        for com in self.u:
            for s_list in com:
                for t in s_list:
                    self.u_start.add(t[0])

        self.graphReversed = self.reverse_graph(G)

        self.E_active = np.ones((self.I, self.m))
        for i in range(self.I):
            self.labels.append([self.dijkstra(self.graphReversed, "t{}".format(i+1), 0, visited=[], distances={})])

        self.E_active = scipy.sparse.lil_matrix((self.I, self.m))
        for (v_ind, v) in enumerate(self.V):
            outneighbors = self.G[v].keys()
            self.delta_p.append([])
            self.delta_p[v_ind] = outneighbors
            for w in outneighbors:
                w_ind = self.V.index(w)
                edge = self.E.index((v, w))
                for i in range(self.I):
                    if abs(self.labels[i][0][v_ind] - self.labels[i][0][w_ind] - self.c[0][edge]) < self.eps:
                        self.E_active[i, edge] = 1

        self.time_vor_main = time.time()
        self.main()

    # Quelle:
    # https://www.geeksforgeeks.org/topological-sorting/#:~:text=Topological%20sorting%20for%20Directed%20Acyclic,4%202%203%201%200%E2%80%9D
    # A recursive function used by topologicalSort
    def topologicalSortUtil(self, v, i, visited, stack):

        # Mark the current node as visited.
        visited[self.V.index(v)] = True

        # Recur for all the vertices adjacent to this vertex
        dp_act = self.get_outgoing_active_edges(i, v)
        for e_ind in dp_act:
            w = self.E[e_ind][1]
            k = self.V.index(w)
            if not visited[k]:
                self.topologicalSortUtil(w, i, visited, stack)

        # Push current vertex to stack which stores result
        stack.append(v)

    # The function to do Topological Sort. It uses recursive
    # topologicalSortUtil()
    def topologicalSort(self, i):
        # Mark all the vertices as not visited
        ti = 't{}'.format(i+1)
        visited = [False]*self.n
        visited[self.V.index(ti)] = True
        stack = [ti]

        # Call the recursive helper function to store Topological
        # Sort starting from all vertices one by one
        for k in range(self.n):
            if not visited[k]:
                self.topologicalSortUtil(self.V[k], i, visited, stack)
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
        while True:
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
            src=min(unvisited, key=unvisited.get)
        # return self.dijkstra(graph, x, phase_ind, visited, distances, predecessors)

    # Quelle: http://www.gilles-bertrand.com/2014/03/dijkstra-algorithm-python-example-source-code-shortest-path.html
    def dijkstraalt(self, graph, src, phase_ind, visited=[], distances={}, predecessors={}):
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
            alpha[i] += self.nu[e_ind] * flow_ratio
            if self.q_global[-1][e_ind] == 0:
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
        v = self.V[v_ind]
        if v == 't{}'.format(i+1):
            return 0
        in_paths = self.get_ingoing_edges(v)
        for e_ind in in_paths:
            ind = self.last_fm_change(i, e_ind, phase)
            """try:
                if abs(self.fm[ind + 1][i][e_ind][0] - phase) < self.eps:
                    ind += 1
            except (IndexError, KeyError):
                continue"""
            b += self.fm[i][e_ind][ind][1]
        u_v = self.u[i][v_ind]
        u_v_len = len(u_v)
        for tuple_ind in range(u_v_len - 1, -1, -1):
            if u_v[tuple_ind][0] <= phase:
                b += u_v[tuple_ind][1]
                break
        return b

    def get_ingoing_edges(self, v):
        """
        bestimmt alle in 'v' eingehenden Kanten
        :param v: Knoten
        :return: Liste der Indizes der Kanten
        """
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

    def get_outgoing_edges(self, v):
        """
        bestimmt alle aus 'v' ausgehenden Kanten
        :param v: Knoten
        :return: Liste der Indizes der Kanten
        """
        gv_keys = self.G[v].keys()
        return [self.E.index((v,u)) for u in gv_keys]

    def get_outgoing_active_edges(self, i, v):
        """
        bestimmt alle aus 'v' ausgehende, momentan für Gut 'i' aktive Kanten
        :param i: Index des Guts
        :param v: Knoten
        :return: Liste der Indizes der aktiven Kanten
        """
        outneighbors = self.G[v].keys()
        delta_p = [self.E.index((v,u)) for u in outneighbors]
        return [e for e in delta_p if self.E_active[i, e]]

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

    def last_fm_change(self, i, e_ind, theta):
        """
        Bestimmt den größten Zeitpunkt <= 'theta', zu dem sich der f_i^- -Wert von Kante 'e_ind' ändert
        :param i: Index des Guts
        :param e_ind: Kantenindex
        :param theta: aktueller Zeitpunkt
        :return: Index des gesuchten Zeitpunkts
        """
        fm_len = len(self.fm[i][e_ind])
        for t in range(fm_len - 1, 0, -1):
            if self.fm[i][e_ind][t][0] < theta + self.eps:
                return t
        return 0

    '''def find_cover_simplex(self, theta_ind):
        simplex = [[], [], []]  # CCS - Format für Eckpunkte des Simplex
        # Dimension: Anzahl Einträge, die ungleich 0 sein können -> len(S[0])
        for e_ind in range(self.m):
            e_start_ind = self.V.index(self.E[e_ind][0])
            simplex[2].append(len(simplex[0]))
            for i in range(self.I):
                if b[i][e_start_ind] > self.eps and self.E_active[i][theta_ind][e_ind]:
                    simplex[0].append(b[i][e_start_ind])  # Wert
                    simplex[1].append(i)  # Zeilenindex

                    """if len(simplex[2]) - 1 < e_ind:
                        dif = e_ind - len(simplex[2]) + 1
                        for c in range(dif):
                            simplex[2].append(len(simplex[0]) - 1)  # Spalten - pointer"""
        simplex[2].append(len(simplex[0]))
        return simplex

    def bary_coords(self, v, x):
        n = len(x)
        if n == 1:
            return [1]
        q, r = np.linalg.qr(v)
        z = q.T.dot(x)
        b = np.zeros(n)
        b[n - 1] = z[n - 1] / r[n - 1, n - 1]
        if b[n - 1] < self.eps:
            b[n - 1] = 0
        for ind in range(n - 2, -1, -1):
            brsum = 0
            for j in range(n - 1, ind, -1):
                brsum += b[j] * r[ind, j]
            b[ind] = (z[ind] - brsum) / r[ind, ind]
            if b[ind] < self.eps:
                b[ind] = 0
        # b = np.linalg.solve(v, x)  # anschauen: multigrid verfahren
        bsum = sum(b)
        for i in range(n):
            b[i] *= 1/bsum
        return b'''

    '''def barycentric_trafo_sp(self, simplex, x):
        n = len(x)
        s1 =
        s2 = np.ones(n)
        b = spsolve(A,x)'''

    '''def find_ancestor_simplex(self, k, p):
        n = len(p)
        q = k * p
        if n == 1:
            return q, [1]
        x = np.zeros(n)
        lam = np.zeros(n-1)
        # perturbiere q, damit q nicht am Rand liegt!
        q_zeros = np.where(q < self.eps)[0]
        epsilon = 0.5
        for i in range(n):
            if q[i] < epsilon and i not in q_zeros:
                epsilon = q[i]
        epsilon *= 0.5
        # epsilon = random() * 0.5
        v_eps = np.zeros(n)
        for i in range(0, n, 2):
            v_eps[i] += epsilon * (10**(-(i % 6) - 1))
        for i in range(1, n, 2):
            if i in q_zeros:
                v_eps[i] += epsilon * (10**(-(i % 6) - 1))
            else:
                v_eps[i] -= epsilon * (10**(-(i % 6) - 1))
        v_sum = sum(v_eps)
        for i in range(n):
            # einfacher?
            if i not in q_zeros:
                v_eps[i] -= v_sum
                break
        q += v_eps

        x[0] = np.floor(q[0] + 1)
        lam[0] = round(x[0] - q[0], 9)
        for i in range(1, n-1):
            x[i] = np.floor(q[i] - lam[i-1] +1)
            lam[i] = round(x[i] - q[i] + lam[i-1], 9)
        x[n-1] = q[n-1] - lam[n-2]

        """
        # if any(x < -self.eps):
        x[0] = np.ceil(q[0])
        lam[0] = round(x[0] - q[0], 9)
        for i in range(1, n-1):
            x[i] = np.ceil(q[i] - lam[i-1])
            lam[i] = round(x[i] - q[i] + lam[i-1], 9)
        x[n-1] = q[n-1] - lam[n-2]"""

        pi = range(1, n)
        sort = sorted(zip(lam, pi),  key=lambda y: (y[0], -y[1]), reverse=True)
        pi = [element for _, element in sort]
        return x, pi'''

    '''def in_k(self, x, con, all_iv, k):
        xk = k * x
        A, b = con
        all_iv_keys = all_iv.keys()
        for (i, v) in all_iv_keys:
            iv_sum = 0
            for d in all_iv[(i, v)]:
                iv_sum += xk[d]
            if iv_sum > self.eps:
                for d in all_iv[(i, v)]:
                    xk[d] *= 1.0 / iv_sum
        if np.linalg.norm(A.dot(xk) - b) < self.eps and min(xk) > -self.eps:  # ?
            return True
        return False'''

    def g(self, e_ind, theta_ind, x):
        dif = x - self.nu[e_ind]
        if abs(dif) < self.bd_tol:
            return 0
        if self.q_global[theta_ind][e_ind] > self.eps:
            return dif
        return max([dif, 0])

    '''def find_fv(self, node, theta_ind, all_iv, simplex):
        rows, col_pointer = simplex[1:]
        dim = len(rows)
        fv = -np.ones(dim)
        g_e = np.zeros(self.m)
        for e_ind in range(self.m):
            if col_pointer[e_ind] - col_pointer[e_ind + 1]:
                v_ind = self.V.index(self.E[e_ind][0])
                xsum = 0
                for ro in range(col_pointer[e_ind], col_pointer[e_ind+1]):
                    # Summe effizienter möglich?
                    if node[ro] > self.eps:
                        xsum += b[rows[ro]][v_ind] * node[ro] / sum([node[d] for d in all_iv[(rows[ro], v_ind)]])
                g_e[e_ind] = self.g(e_ind, theta_ind, xsum)
            else:
                g_e[e_ind] = self.g(e_ind, theta_ind, 0)
        for i in range(self.I):
            act_edges = [e_ind for e_ind in range(self.m) if self.E_active[i][theta_ind][e_ind]]
            top_ord = self.topologicalSort(i)[1:]
            a = np.zeros(self.n)
            for v in top_ord:
                v_ind = self.V.index(v)
                delta_p_act = self.get_outgoing_active_edges(i, v)
                a[v_ind] = np.Inf
                for e_ind in delta_p_act:
                    w_ind = self.V.index(self.E[e_ind][1])
                    if a[v_ind] > g_e[e_ind] + a[w_ind]:
                        a[v_ind] = g_e[e_ind] + a[w_ind]
                if not a[v_ind] < np.Inf:
                    a[v_ind] = 0

            for e_ind in act_edges:
                v, w = self.E[e_ind]
                v_ind = self.V.index(v)
                w_ind = self.V.index(w)
                if b[i, v_ind] < self.eps:
                    continue
                if g_e[e_ind]/self.nu[e_ind] + a[w_ind] - a[v_ind] > self.eps:
                    col_start = col_pointer[e_ind]
                    col_end = col_pointer[e_ind+1]
                    col_skip = len([r for r in rows[col_start:col_end] if r < i])
                    fv[col_start + col_skip] = 0
        co = 0
        for d in range(dim):
            if fv[d] > -self.eps:
                continue
            try:
                while col_pointer[co+1] <= d:
                    co += 1
            except IndexError:
                pass
            if not self.E_active[rows[d]][-1][co]:
                fv[d] = 0
                continue
            v_ind = self.V.index(self.E[co][0])
            all_iv_keys = all_iv.keys()
            v_count = len([v for (i, v) in all_iv_keys if i == rows[d]])
            fv[d] = 1.0 / (self.I * v_count)
            delta_p = self.get_outgoing_active_edges(rows[d], self.V[v_ind])
            delta_p.remove(co)
            for e_ind in delta_p:
                try:
                    elt_ind = rows[col_pointer[e_ind]:col_pointer[e_ind + 1]].index(rows[d]) + col_pointer[e_ind]
                    fv[elt_ind] = 0
                except ValueError:
                    continue
        return fv

    def fk_on_tk(self, x0, pi, con, simplex, all_iv, theta_ind, c, k, part=None):

        dim = len(simplex[0])
        if dim == 1:
            return x0, [1]
        if part is None:
            part = range(dim)
        U1 = np.eye(dim-1, k=-1) - np.eye(dim-1)
        U2 = np.zeros(dim-1)
        U2[-1] = 1
        U = np.concatenate((U1, U2.reshape((1,dim-1))))
        X = np.zeros((dim, dim))
        X[:, 0] = 1.0 / k * x0
        # X[:, 0] = x0
        f_X = -np.ones((dim, dim))
        for j in part:
            if j > 0:
                if len(part) < dim:
                    X[:, j] = X[:, 0]
                    for uj in range(j):
                        X[:, j] += (1.0 / k) * U[:, pi[uj] - 1]
                else:
                    X[:, j] = X[:, j-1] + (1.0 / k) * U[:, pi[j-1] - 1]
            if self.in_k(X[:, j], con, all_iv, k):
                f_X[:, j] = self.find_fv(X[:, j], theta_ind, all_iv, simplex)
            else:
                f_X[:, j] = c
        if len(part) == 1:
            return [xi[0] for xi in X[:, part]], [fxi[0] for fxi in f_X[:, part]]
        return X[:, part], f_X[:, part]

    def compute_fp(self, simplex, con, theta_ind):

        erglistx0 = []
        erglistpi = []
        [vals, rows, col_pointer] = simplex
        dim = len(vals)
        if dim == 0:
            return 0, {}
        cs_bz = np.zeros(dim)
        co = 0
        all_iv = {}
        for d in range(dim):
            ro = rows[d]
            try:
                while d >= col_pointer[co + 1]:
                    co += 1
            except IndexError:
                pass
            v_ind = self.V.index(self.E[co][0])
            len_dp_act = len(self.get_outgoing_active_edges(rows[d], self.V[v_ind], theta_ind=theta_ind))
            cs_bz[d] = 1 / len_dp_act
            if (ro, v_ind) not in all_iv.keys():
                all_iv[(ro, v_ind)] = [d]
            else:
                all_iv[(ro, v_ind)].append(d)

        k = 3
        # transformiere 'cs_bz' auf Standardsimplex
        cs_bz *= 1.0 / len(all_iv.keys())
        epsilon = 0.5
        for i in range(dim):
            if self.eps < cs_bz[i] < epsilon:
                epsilon = cs_bz[i]
        epsilon = 1/(2*k+2)
        epsilon *= 0.5
        #epsilon = random() * 0.5
        v_eps = np.zeros(dim)
        for i in range(0, dim, 2):
            v_eps[i] += epsilon**(i+1)
        for i in range(1, dim, 2):
            v_eps[i] -= epsilon**(i+1)
        v_eps_sum = sum(v_eps)
        for i in range(dim):
            if cs_bz[i] + v_eps[i] > v_eps_sum:
                v_eps[i] -= v_eps_sum
                break
        #cs_bz += v_eps
        fpk = []
        while True:
            if len(fpk) > 2 and np.linalg.norm(fpk[-2] - fpk[-1]) < self.eps:
                return fpk[-1], all_iv
            if k == 3:
                start_node = cs_bz
            else:
                start_node = fpk[-1]
            #start_node = np.zeros(dim)
            #start_node[-1] = 1
            xc0 , pi_c = self.find_ancestor_simplex(k, start_node)
            Xc, f_Xc = self.fk_on_tk(xc0, pi_c, con, simplex, all_iv, theta_ind, cs_bz, k)

            # bestimme baryzentrische Koordinaten von cs_bz bzgl. Xc[:, 0], ..., Xc[:, dim]
            bz = bary_coords(Xc, start_node)
            """
            while np.any(bz < self.eps):
                bz0 = np.where(bz < self.eps)[0][0]
                cs_bz[bz0] += 1/(k+3) 
                for i in range(dim):
                    if cs_bz[i] > 1/(k+3) + self.eps and i != bz0:
                        cs_bz[i] -= 1/(k+3)
                        break
                bz = bary_coords(Xc, cs_bz)"""
            # berechne f_k(c) mit den bary. Koordinaten von c bzgl. X und deren Funktionswerten f_k(X)
            fk_c = np.zeros(dim)
            if dim == 1:
                fk_c[0] = f_Xc[0]
            else:
                for i in range(dim):
                    fk_c += bz[i] * f_Xc[:, i]

            r = fk_c - start_node  # cs = c + kappa * r, mit kappa max., s.d. cs >= 0.
            if np.max(r) < self.eps:
                fpk.append(start_node)
                k += 2
                continue
            epsilon = 0.5
            for i in range(dim):
                if self.eps < abs(r[i]) < epsilon:
                    epsilon = abs(r[i])
            epsilon *= 0.5
            v_eps = np.zeros(dim)
            """
            if dim % 2:
                # wähle 'sign' so, dass v_eps[dim] > 0. (Weil sonst r[dim] < 0 und somit später kappa = 0)
                sign = -1
            else:
                sign = 1
            for i in range(0, dim+1, 2):
                v_eps[i] = sign * epsilon/(i+1)
            for i in range(1, dim+1, 2):
                v_eps[i] = -sign * epsilon/i
            # falls 'dim' gerade, so wird zweite 'for'-Schleife einmal weniger durchlaufen als die erste
            if not dim % 2:
                v_eps[dim-1] -= epsilon / (dim+1)
            """
            start_node_zeros = np.where(start_node < self.eps)[0]

            for i in range(0, dim, 2):
                v_eps[i] = epsilon * (10**(-(i % 6) - 1))
            for i in range(1, dim, 2):
                if i in start_node_zeros:
                    v_eps[i] = epsilon * (10**(-(i % 6) - 1))
                else:
                    v_eps[i] = -epsilon * (10**(-(i % 6) - 1))
            v_eps_sum = sum(v_eps)
            if v_eps_sum > self.eps:
                # einfacher?
                for i in range(dim):
                    if i not in start_node_zeros:
                        v_eps[i] -= v_eps_sum
                        break
            """
            for i in range(0, dim, 2):
                v_eps[i] = epsilon**(i+1)
            for i in range(1, dim, 2):
                if i in cs_bz_zeros:
                    v_eps[i] = epsilon**i
                else:
                    v_eps[i] = -epsilon**i
            v_eps_sum = sum(v_eps)
            if v_eps_sum > self.eps:
                # einfacher?
                for i in range(dim):
                    if i not in cs_bz_zeros:
                        v_eps[i] -= v_eps_sum
                        break"""
            # Perturbation von 'r'
            r += v_eps
            # perturbation von r so nicht korrekt!? da epsilon so klein wird r nicht ausreichend perturbiert!?
            #r = np.zeros(dim)
            kappa = 1/self.eps
            for i in range(dim):
                if start_node[i] + kappa * r[i] < -self.eps:
                    kappa = -start_node[i] / r[i]
            cs = start_node + kappa * r
            x0, pi = self.find_ancestor_simplex(k, cs)
            X, f_X = self.fk_on_tk(x0, pi, con, simplex, all_iv, theta_ind, cs_bz, k)

            L = np.zeros((dim+1, dim))
            for d in range(dim):
                L[:dim, d] = f_X[:, d] - X[:, d]
                L[dim, d] = 1

            tf = np.concatenate((np.zeros(dim+1), np.ones(dim+1)))
            r = np.append(r, 1).reshape((dim+1, 1))
            A_total = np.concatenate((L, r, np.eye(dim+1)), axis=1)
            b_total = np.zeros(dim+1)
            b_total[-1] = 1
            basic_sol = linprog(tf, A_eq=A_total, b_eq=b_total, method='simplex') #, options={'bland': True})
            y = basic_sol.x[:dim]
            z = basic_sol.x[dim]
            fun = basic_sol.fun

            # aktueller Simplex hat keine zulässige Lsg. Siehe CEavesAnOddThm: S0 als (n-1)-Simplex enthält bfs. Aber
            # wie
            # ist das hier zu interpretieren, wo ich doch die dim schon um 1 verringert habe? -> 0 wieder zu
            # überdeckungssimplex hinzu nehmen!? dann v0 aus paper gleich 0: f(v0) - v0 = f(v0) = csbz ...
            # v1, ..., vn sind dann die knoten dich ich hier schon habe
            if fun > self.eps:
                raise TypeError('Basic feasible solution can\'t be found')
            if z < self.eps:
                xk = np.zeros(dim)
                for d in range(dim):
                    xk += y[d] * X[:, d]
                fpk.append(xk)
                k += 2
                continue

            s = np.where(y < self.eps)[0][0]
            while True:
                if s == dim-1:
                    L_s = L[:, :s]
                elif s == 0:
                    L_s = L[:, 1:]
                else:
                    L_s = np.concatenate((L[:, :s], L[:, s+1:]), axis=1)
                L_s = np.concatenate((L_s, r), axis=1)
                q, rtriangle = np.linalg.qr(L_s)
                qLs = q.T.dot(L[:, s])
                mu = np.zeros(dim)
                mu[dim-1] = qLs[dim-1] / rtriangle[dim-1, dim-1]
                if abs(mu[dim-1]) < self.eps:
                    mu[dim-1] = 0
                for ind in range(dim - 2, -1, -1):
                    if abs(rtriangle[ind, ind]) < self.eps:
                        # ???
                        # mu[ind] = -1
                        continue
                    mursum = 0
                    for j in range(dim - 1, ind, -1):
                        mursum += mu[j] * rtriangle[ind, j]
                    mu[ind] = (qLs[ind] - mursum) / rtriangle[ind, ind]
                    if abs(mu[ind]) < self.eps:
                        mu[ind] = 0
                #lam = np.linalg.solve(L_s[:-1, :], L[:, s])  # Darstellung von 'L[:, s]' als Linearkomb. der Spalten von 'L_s'
                #lam = np.linalg.solve(L_s, L[:, s])  # Darstellung von 'L[:, s]' als Linearkomb. der Spalten von 'L_s'

                delt = 1/self.eps
                t = 1/self.eps
                y_range = list(set(range(dim)) - {s})
                for col in y_range:
                    if col > s:
                        mucol = col - 1
                    else:
                        mucol = col
                    if y[col] - mu[mucol] * delt < -self.eps:
                        delt = y[col]/mu[mucol]
                        t = col
                if z - mu[dim-1] * delt < -self.eps:
                    delt = z/mu[dim-1]

                delta_y = np.concatenate((list(-mu[:s] * delt), [delt], list(-mu[s:-1] * delt)))
                delta_z = -mu[-1] * delt
                yss = y + delta_y
                zss = z + delta_z

                if zss < self.eps:
                    xk = np.zeros(dim)
                    for d in range(dim):
                        xk += yss[d] * X[:, d]
                    fpk.append(xk)
                    erglistx0.append(x0)
                    erglistpi.append(pi)
                    break
                    # return x0, pi, yss

                L_rows = L.shape[0]
                if t == 0:
                    x0[pi[0]-1] -= 1
                    x0[pi[0]] += 1
                    j1 = pi.pop(0)
                    pi.append(j1)
                    yss = np.delete(yss, 0)
                    yss = np.append(yss, 0)
                    xt, fxt = self.fk_on_tk(x0, pi, con, simplex, all_iv, theta_ind, cs_bz, k, [dim-1])
                    lt = fxt - xt
                    lt = np.append(lt, 1)
                    L = np.concatenate((L[:, 1:], lt.reshape((L_rows, 1))), axis=1)
                    s = dim - 1
                elif t == dim - 1:
                    x0[pi[-1]-1] += 1
                    x0[pi[-1]] -= 1
                    jn = pi.pop(-1)
                    pi.insert(0, jn)
                    yss = np.delete(yss, -1)
                    yss = np.insert(yss, 0, 0)
                    xt, fxt = self.fk_on_tk(x0, pi, con, simplex, all_iv, theta_ind, cs_bz, k, [0])
                    lt = fxt - xt
                    lt = np.append(lt, 1)
                    L = np.concatenate((lt.reshape((L_rows, 1)), L[:, :-1]), axis=1)
                    s = 0
                else:
                    """
                    # table 1 falsch!? eigener Versuch: Vertausche Spalten in L:
                    L_old = L.copy()
                    if t < dim - 2:
                        L = np.concatenate((L_old[:, :t], L_old[:, t+1].reshape((L_rows, 1)), L_old[:, t].reshape((L_rows, 1)), L_old[:, t+2:]), axis=1)
                    else:
                        # t = dim - 2
                        L = np.concatenate((L_old[:, :t], L_old[:, t+1].reshape((L_rows, 1)), L_old[:, t].reshape((L_rows, 1))), axis=1)
                    pi_old = pi.copy()
                    pi[t-1] = pi[t]
                    pi[t] = pi_old[t-1]
                    xt, fxt = self.fk_on_tk(x0, pi, con, simplex, all_iv, theta_ind, cs_bz, k, [t])
                    lt = fxt - xt
                    lt = np.append(lt, 1)
                    L[:, t] = lt
                    s = t
                y = yss
                z = zss

            k += 2

    def lt_algo(self, simplex, con, theta_ind):
        [vals, rows, col_pointer] = simplex
        dim = len(vals)
        if dim == 0:
            return 0, {}
        cs_bz = np.zeros(dim)
        co = 0
        all_iv = {}
        for d in range(dim):
            ro = rows[d]
            try:
                while d >= col_pointer[co + 1]:
                    co += 1
            except IndexError:
                pass
            v_ind = self.V.index(self.E[co][0])
            len_dp_act = len(self.get_outgoing_active_edges(rows[d], self.V[v_ind], theta_ind=theta_ind))
            cs_bz[d] = 1 / len_dp_act
            if (ro, v_ind) not in all_iv.keys():
                all_iv[(ro, v_ind)] = [d]
            else:
                all_iv[(ro, v_ind)].append(d)

        k = 3
        # transformiere 'cs_bz' auf Standardsimplex
        cs_bz *= 1.0 / len(all_iv.keys())
        clss = []

        while True:
            if len(clss) > 2 and np.linalg.norm(clss[-1][-1] - clss[-2][-1]) < self.eps:
                return clss[-1]
            v0, pi = self.find_ancestor_simplex(k, cs_bz)
            # Im Fall 'part=[0]' ist 'pi' irrelevant
            x, fv0 = self.fk_on_tk(v0, [], con, simplex, all_iv, theta_ind, cs_bz, k, part=[0])
            lv0 = fv0 - v0 / k + np.ones(dim)
            base = np.eye(dim)
            base_v = np.zeros((dim, dim))
            base_y = np.ones(dim)
            wt_next = v0 / k
            T = []
            gamma = []
            R = np.zeros(dim)
            L = [lv0]
            tauw = [v0 / k]
            t = 0
            mu_no = dim
            while True:
                if mu_no == 0:
                    clss.append(tauw)
                    k += 1
                    break
                # Im Fall 'part=[0]' ist 'pi' irrelevant
                fvt_next = self.fk_on_tk(k * wt_next, [], con, simplex, all_iv, theta_ind, cs_bz, k, part=[0])[1]
                lvt_next = fvt_next - wt_next + np.ones(dim)
                # pivotiere lvt_next in base
                q, r = np.linalg.qr(base)
                z = q.T.dot(lvt_next)
                lam = np.zeros(dim)
                lam[dim-1] = z[dim-1] / r[dim-1, dim-1]
                for ind in range(dim-2, -1, -1):
                    lamrsum = 0
                    for j in range(dim-1, ind, -1):
                        lamrsum += lam[j] * r[ind, j]
                    lam[ind] = (z[ind] - lamrsum) / r[ind, ind]
                delt = 1/self.eps
                out = 1/self.eps
                for j in range(dim):
                    if base_y[j] - lam[j] * delt < -self.eps:
                        delt = base_y[j] / lam[j]
                        out = j
                for j in range(dim):
                    base_y[j] -= delt * lam[j]
                base_y[out] = delt
                if len(np.where(base[:, out] > self.eps)[0]) > 1:
                    wt_old = wt_next.copy()
                    wt_next = base_v[:, out]
                    m = 0
                    for wi in range(t):
                        if np.linalg.norm(wt_next - tauw[wi]) < self.eps:
                            m = wi
                            break
                    base[:, out] = lvt_next.reshape(dim)
                    base_v[:, out] = wt_old
                else:
                    mu_no -= 1
                    qi = np.zeros(dim)
                    i = np.where(base[:, out] > self.eps)[0][0]
                    qi[i] = -1
                    qi[i+1] = 1
                    tauw.append(tauw[-1] + qi / k)
                    base[:, out] = lvt_next.reshape(dim)
                    base_v[:, out] = wt_next
                    wt_next = tauw[-1]
                    T.append(i)
                    gamma.append(i)
                    t += 1
                    continue
                # if lvt_next in L:
                    # s = L.index(lvt_next)
                while True:
                    if m == 0:
                        qgamma1 = np.zeros(dim)
                        qgamma1[gamma[0]] = -1 / k
                        qgamma1[gamma[0]+1] = 1 / k
                        tauw[0] += qgamma1
                        g0 = gamma[0]
                        gamma = np.delete(gamma, 0)
                        gamma = np.append(gamma, g0)
                        R[g0] += 1
                        qj = np.zeros(dim)
                        qj[gamma[0]] = -1 / k
                        qj[gamma[0]+1] = 1 / k
                        tauw.append(tauw[-1] + qj)
                        wt_next = tauw[-1]
                        del tauw[0]
                        del L[0]
                        L.append(lvt_next)
                        break
                    elif m == t:
                        gt = gamma[-1]
                        if R[gt] < self.eps:  # 'R[gt]' wird negativ
                            R[gt] = 0
                            gamma = np.delete(gamma, -1)
                            T.pop()
                            lam = np.zeros(dim)
                            q, r = np.linalg.qr(base)
                            z = q.T[:, gt]
                            lam[dim-1] = z[dim-1] / r[dim-1, dim-1]
                            for ind in range(dim-2, -1, -1):
                                lamrsum = 0
                                for j in range(dim-1, ind, -1):
                                    lamrsum += lam[j] * r[ind, j]
                                lam[ind] = (z[ind] - lamrsum) / r[ind, ind]
                            delt = 1/self.eps
                            out = 1/self.eps
                            for j in range(dim):
                                if base_y[j] * base[j] - lam[j] * delt < -self.eps:
                                    delt = base_y[j] / lam[j]
                                    out = j
                            for j in range(dim):
                                base_y[j] -= delt * lam[j]
                            base_y[out] = delt
                            if len(np.where(base[:, out] > self.eps)[0]) > 1:
                                wt_next = base_v[:, out]
                                base[:, out] = np.zeros(dim)
                                base[gt, out] = 1
                            else:
                                mu_no -= 1
                                qi = np.zeros(dim)
                                i = np.where(base[:, out] > self.eps)[0][0]
                                base[:, out] = np.zeros(dim)
                                base[gt, out] = 1
                                qi[i] = -1
                                qi[i+1] = 1
                                tauw.append(tauw[-1] + qi / k)
                                wt_next = tauw[-1]
                                T.append(i)
                                gamma.append(i)
                                t += 1
                                break
                        else:
                            qgammat = np.zeros(dim)
                            gamma = np.delete(gamma, -1)
                            gamma = np.insert(gamma, 0, gt)
                            qgammat[gt] = -1 / k
                            qgammat[gt + 1] = 1 / k
                            tauw[0] -= qgammat
                            R[gt] -= 1
                            del tauw[-1]
                            tauw.append(tauw[0])
                            del L[-1]
                            L.append(lvt_next)
                            break
                    else:
                        gs = gamma[m]
                        gamma[m] = gamma[out-1]
                        gamma[m-1] = gs
                        qj = np.zeros(dim)
                        qj[gamma[m]] = -1 / k
                        qj[gamma[m]+1] = 1 / k
                        tauw.append(tauw[m] + qj)
                        del tauw[m]
                        del L[m]
                        L.append(lvt_next)
                        break
                # else:
                #    tauw.append(wt_next)
                #    L.append(lvt_next)
                #    T.append(lvt_next)
                #    gamma.append(lvt_next)'''

    def gamma(self, x, x_sum, a, delta_p_act, coms, theta_ind):
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
        gx = scipy.sparse.lil_matrix((self.I, self.m))
        forced_zeros = []
        for i in range(self.I):
            forced_zeros.append({})
        for v_ind in coms:
            for i in coms[v_ind]:
                forced_zeros[i][v_ind] = []
                x_sum_nz = 0
                cor = 0
                for e_ind in delta_p_act[i][v_ind]:
                    w_ind = self.V.index(self.E[e_ind][1])
                    if abs(a[i, v_ind] - self.g(e_ind, theta_ind, x_sum[e_ind]) / self.nu[e_ind] - a[i, w_ind]) > self.eps:
                        # 'gx[i, e_ind]' = 0 muss gelten
                        forced_zeros[i][v_ind].append(e_ind)
                        # merke Unterschied zwischen 'x[i, e_ind]' und 'gx[i, e_ind]' um später Flusserhaltung wiederherzustellen
                        cor += x[i, e_ind]
                    else:
                        # merke bereits festgelegte Flussmenge
                        x_sum_nz += x[i, e_ind]
                dp_nz = list(set(delta_p_act[i][v_ind]) - set(forced_zeros[i][v_ind]))
                for e_ind in dp_nz:
                    # Verteile übriges Flussvolumen ('cor' viel) auf aktive, nicht-0-Kanten im Verhältnis der bereits zugeteilten Flusswerte
                    gx[i, e_ind] = x[i, e_ind] * (1 + cor / x_sum_nz)
        return gx, forced_zeros

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
        nu_sum_act = []
        coms = defaultdict(list)
        self.b.append(scipy.sparse.lil_matrix((self.I, self.n)))
        self.flow_vol.append(scipy.sparse.lil_matrix((self.n, 1)))

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
                    cor = 0
                    e_diff = {}
                    e_diff_sum = 0
                    e_cor = []
                    if flow_sent[i] > b_nf[i, v_ind] + self.eps:
                        for e_ind in e_bds_nf:
                            x_sum[e_ind] -= xk[i, e_ind]
                            xk[i, e_ind] *= b_nf[i, v_ind] / flow_sent[i]
                            # befindet sich nun ein 'xk[i, e_ind]' - Wert außerhalb seines 'bounds' - Intervalls, so wird dieser Wert an die Grenze des Intervalls gesetzt
                            # und diese Korrektur später durch die Kanten in 'e_cor' ausgeglichen, um Flusserhaltung beizubehalten
                            if xk[i, e_ind] < bounds[v_ind][i][e_ind][0]:
                                cor += bounds[v_ind][i][e_ind][0] - xk[i, e_ind]
                                xk[i, e_ind] = bounds[v_ind][i][e_ind][0]
                            else:
                                e_cor.append(e_ind)
                            x_sum[e_ind] += xk[i, e_ind]

                        if len(e_cor) < len(e_bds_nf):
                            # mindestens eine bound wurde verletzt
                            for e_ind in e_cor:
                                e_diff[e_ind] = xk[i, e_ind] - bounds[v_ind][i][e_ind][0]
                                e_diff_sum += e_diff[e_ind]
                            for e_ind in e_cor:
                                x_cor = cor * e_diff[e_ind] / e_diff_sum
                                xk[i, e_ind] -= x_cor
                                x_sum[e_ind] -= x_cor
                    elif flow_sent[i] < b_nf[i, v_ind] - self.eps:
                        for e_ind in e_bds_nf:
                            x_sum[e_ind] -= xk[i, e_ind]
                            xk[i, e_ind] *= b_nf[i, v_ind] / flow_sent[i]
                            # befindet sich nun ein 'xk[i, e_ind]' - Wert außerhalb seines 'bounds' - Intervalls, so wird dieser Wert an die Grenze des Intervalls gesetzt
                            # und diese Korrektur später durch die Kanten in 'e_cor' ausgeglichen, um Flusserhaltung beizubehalten
                            if xk[i, e_ind] > bounds[v_ind][i][e_ind][1]:
                                cor += xk[i, e_ind] - bounds[v_ind][i][e_ind][1]
                                xk[i, e_ind] = bounds[v_ind][i][e_ind][1]
                            else:
                                e_cor.append(e_ind)
                            x_sum[e_ind] += xk[i, e_ind]

                        if len(e_cor) < len(e_bds_nf):
                            # mindestens eine bound wurde verletzt
                            for e_ind in e_cor:
                                e_diff[e_ind] = bounds[v_ind][i][e_ind][1] - xk[i, e_ind]
                                e_diff_sum += e_diff[e_ind]
                            for e_ind in e_cor:
                                x_cor = cor * e_diff[e_ind] / e_diff_sum
                                xk[i, e_ind] += x_cor
                                x_sum[e_ind] += x_cor
            return xk, x_sum

        # Initialisiere x0
        for i in range(self.I):
            delta_p_act.append({})
            nu_sum_act.append({})
            for v_ind in range(self.n):
                v = self.V[v_ind]
                delta_p_act[i][v_ind] = self.get_outgoing_active_edges(i, v)
                self.b[-1][i, v_ind] = self.calc_b(i, v_ind, theta)
                if self.b[-1][i, v_ind] > self.bd_tol:
                    self.flow_vol[-1][v_ind, 0] += self.b[-1][i, v_ind]
                    if len(delta_p_act[i][v_ind]) == 1:
                        e_ind = delta_p_act[i][v_ind][0]
                        xfp[i, e_ind] = self.b[-1][i, v_ind]
                        x_sum_fix[e_ind] += xfp[i, e_ind]
                        continue
                    coms[v_ind].append(i)
                    nu_sum_act[i][v_ind] = 0
                    for e_ind in delta_p_act[i][v_ind]:
                        nu_sum_act[i][v_ind] += self.nu[e_ind]
                    for e_ind in delta_p_act[i][v_ind]:
                        xk[i, e_ind] = self.b[-1][i, v_ind] * self.nu[e_ind] / nu_sum_act[i][v_ind]
                        x_sum[e_ind] += xk[i, e_ind]

        # Teilmenge von 'delta_p_act', welche für jedes Gut nur die Kanten mit nicht bereits fixiertem Flusswert enthält
        delta_p_act_nf = copy.deepcopy(delta_p_act)

        # Bestimme a_i,v - Werte zu x0
        a = scipy.sparse.lil_matrix((self.I, self.n))
        top_ords = []
        sparse_items = []
        for i in range(self.I):
            t_ind = self.V.index('t{}'.format(i+1))
            a[i, t_ind] = 10**(-14)  # < 'self.eps', wird also als 0 behandelt
            sparse_items.append((i, t_ind))
        argmin = []
        for i in range(self.I):
            argmin.append({})
            top_ords.append(self.topologicalSort(i))
            for v in top_ords[i]:
                v_ind = self.V.index(v)
                a_min2 = {}
                for e_ind in delta_p_act[i][v_ind]:
                    w_ind = self.V.index(self.E[e_ind][1])
                    a_e = self.g(e_ind, theta_ind, x_sum[e_ind] + x_sum_fix[e_ind]) / self.nu[e_ind] + a[i, w_ind]
                    if (i, v_ind) not in sparse_items:
                        a[i, v_ind] = a_e
                        sparse_items.append((i, v_ind))
                        argmin[i][v_ind] = [e_ind]
                    # elif a_e < a[i, v_ind] + self.eps:
                    # ???
                    elif a_e < a[i, v_ind] + self.bd_tol:
                        if a_e < a[i, v_ind] - self.bd_tol:
                            a_min2[e_ind] = a[i, v_ind]
                            a[i, v_ind] = a_e
                            argmin[i][v_ind] = [e_ind]
                        else:
                            argmin[i][v_ind].append(e_ind)

        coms_keys = list(coms)
        for v_ind in coms_keys:
            # für x0 ist 'argmin' = 'argmin_nf' und 'delta_p_act' = 'delta_p_act_nf'
            if np.all([len(argmin[i][v_ind]) == len(delta_p_act_nf[i][v_ind]) for i in coms[v_ind]]):
                outneighbors = list(set([self.V.index(self.E[e_ind][1]) for e_ind in argmin[i][v_ind] for i in coms[v_ind]]))
                # sind zusätzlich für alle Endknoten 'w_ind' die Flussaufteilungen bereits fixiert, kann auch die Aufteilung in Knoten 'v_ind' fixiert werden
                if np.all([w_ind not in coms for w_ind in outneighbors]):
                    for i in coms[v_ind]:
                        for e_ind in delta_p_act_nf[i][v_ind]:
                            if xk[i, e_ind] < self.bd_tol:
                                xfp[i, e_ind] = 0
                            else:
                                xfp[i, e_ind] = xk[i, e_ind]
                                x_sum_fix[e_ind] += xfp[i, e_ind]
                        delta_p_act_nf[i][v_ind] = []
                    del coms[v_ind]
                    if not coms:
                        # Abbruchbedingung: keine Güter mehr übrig, die noch auf Kanten aufgeteilt werden müssen
                        return xfp, x_sum_fix, a, top_ords

        # 'bounds' enthält untere und obere Schranken für 'xfp[i, e_ind]' - Werte
        bounds = {}
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
            # print("XK", xk)
            # nur im ersten Schritt nicht erfüllt, danach enthält 'fp_comp' Tupel der Form '(i, v_ind)', was impliziert, dass die Aufteilung von Gut 'i' im Knoten 'v_ind'
            # vollständig fixiert wurde (d.h. die Boundsintervalle aller aktiven ausgehenden Kanten sind hinreichend klein)
            for (i, v_ind) in fp_comp:
                for e_ind in bounds[v_ind][i]:
                    xfp[i, e_ind] = (bounds[v_ind][i][e_ind][0] + bounds[v_ind][i][e_ind][1]) * 0.5
                    if xfp[i, e_ind] < self.bd_tol:
                        xfp[i, e_ind] = 0
                    else:
                        x_sum_fix[e_ind] += xfp[i, e_ind]
                    del bounds[v_ind]
                    del bounds_f[v_ind]
                    coms[v_ind].remove(i)
                    if not coms[v_ind]:
                        del coms[v_ind]
            """
            if fp_comp:
                for (i, e_ind) in fp_comp:
                    v_ind = self.V.index(self.E[e_ind][0])
                    xfp[i, e_ind] = (bounds[v_ind][i][e_ind][0] + bounds[v_ind][i][e_ind][1]) * 0.5
                    if xfp[i, e_ind] < self.bd_tol:
                        xfp[i, e_ind] = 0
                    x_sum_fix[e_ind] += xfp[i, e_ind]
                    delta_p_act_nf[i][v_ind].remove(e_ind)
                    nu_sum_act[i][v_ind] -= self.nu[e_ind]
                    b[i, v_ind] -= xfp[i, e_ind]
                    # ist nur noch eine aktive, nicht-fixierte Kante übrig, so erhält diese den Rest des Flusses von Gut 'i'
                    if len(delta_p_act_nf[i][v_ind]) == 1:
                        e_last = delta_p_act_nf[i][v_ind][0]
                        if (i, e_last) not in fp_comp:
                            if b[i, v_ind] > self.bd_tol:
                                xfp[i, e_last] = b[i, v_ind]
                                x_sum_fix[e_last] += xfp[i, e_last]
                            b[i, v_ind] = 0
                    del bounds[v_ind][i][e_ind]
                    if b[i, v_ind] < self.bd_tol:
                        coms[v_ind].remove(i)
                        if not coms[v_ind]:
                            del coms[v_ind]
            """
            """
                for e_ind in fp_comp:
                    v_ind = self.V.index(self.E[e_ind][0])
                    approx = (bounds[v_ind][e_ind][0] + bounds[v_ind][e_ind][1]) * 0.5
                    vol_per_cap_sum = 0
                    for i in coms[v_ind]:
                        if self.E_active[i][-1][e_ind]:
                            vol_per_cap_sum += b[i, v_ind] / nu_sum_act[i][v_ind]
                    for i in coms[v_ind]:
                        if self.E_active[i][-1][e_ind]:
                            xfp[i, e_ind] = approx * (b[i, v_ind] / nu_sum_act[i][v_ind]) / vol_per_cap_sum
                            x_sum_fix[e_ind] += xfp[i, e_ind]
                            delta_p_act_nf[i][v_ind].remove(e_ind)
                            nu_sum_act[i][v_ind] -= self.nu[e_ind]
                            # ist nur noch eine aktive, nicht-fixierte Kante übrig, so erhält diese den Rest des Flusses von Gut 'i'
                            if len(delta_p_act_nf[i][v_ind]) == 1:
                                e_last = delta_p_act_nf[i][v_ind][0]
                                if e_last not in fp_comp:  # ist es überhaupt möglich, dass dies 'False' ist?
                                    xfp[i, e_last] = b[i, v_ind] - xfp[i, e_ind]
                                    flow_vol_nf[v_ind] -= xfp[i, e_last]
                                    x_sum_fix[e_last] += xfp[i, e_last]
                                    b[i, v_ind] = xfp[i, e_ind]
                    del bounds[v_ind][e_ind]
                for e_ind in fp_comp:
                    v_ind = self.V.index(self.E[e_ind][0])
                    coms_vind = copy.deepcopy(coms[v_ind])
                    for i in coms_vind:
                        b[i, v_ind] -= xfp[i, e_ind]
                        flow_vol_nf[v_ind] -= xfp[i, e_ind]
                        if b[i, v_ind] < self.eps:
                            coms[v_ind].remove(i)
                            if not coms[v_ind]:
                                del coms[v_ind]
            """

            if not coms:  # Abbruchbedingung: keine Güter mehr übrig, die noch auf Kanten aufgeteilt werden müssen
                return xfp, x_sum_fix, a, top_ords

            all_bds = np.concatenate([[bounds[v_ind][i] for i in coms[v_ind]] for v_ind in coms])
            if np.any([len(d) for d in all_bds]):
                xk, x_sum = calc_flow_by_bounds()
            else:
                xk = scipy.sparse.lil_matrix((self.I, self.m))
                x_sum = np.zeros(self.m)
                # setze 'xk'-Werte der nicht-fixierten Kanten
                for v_ind in coms:
                    for i in coms[v_ind]:
                        for e_ind in delta_p_act[i][v_ind]:
                            # ??? nu_sum_act oder nu_sum_act_nf einführen ???
                            xk[i, e_ind] = b_nf[i, v_ind] * self.nu[e_ind] / nu_sum_act[i][v_ind]
                            x_sum[e_ind] += xk[i, e_ind]

            # berechne zugehörige a-Werte
            a = scipy.sparse.lil_matrix((self.I, self.n))
            sparse_items = []
            for i in range(self.I):
                t_ind = self.V.index('t{}'.format(i+1))
                a[i, t_ind] = 10**(-14)  # < 'self.eps', wird also als 0 behandelt
                sparse_items.append((i, t_ind))
            argmin = []
            argmin_nf = []
            for i in range(self.I):
                argmin.append({})
                argmin_nf.append({})
                for v in top_ords[i]:
                    v_ind = self.V.index(v)
                    a_min2 = {}
                    for e_ind in delta_p_act[i][v_ind]:
                        w_ind = self.V.index(self.E[e_ind][1])
                        a_e = self.g(e_ind, theta_ind, x_sum[e_ind] + x_sum_fix[e_ind]) / self.nu[e_ind] + a[i, w_ind]
                        if (i, v_ind) not in sparse_items:
                            a[i, v_ind] = a_e
                            sparse_items.append((i, v_ind))
                            argmin[i][v_ind] = [e_ind]
                            if e_ind in delta_p_act_nf[i][v_ind]:
                                argmin_nf[i][v_ind] = [e_ind]
                            else:
                                # min wird in bereits fixierter Kante angenommen. 'argmin_nf' enthält jedoch nur nicht-fixierte, minimale Kanten
                                argmin_nf[i][v_ind] = []
                        elif a_e < a[i, v_ind] + self.bd_tol:
                            if a_e < a[i, v_ind] - self.bd_tol:
                                a_min2[e_ind] = a[i, v_ind]
                                a[i, v_ind] = a_e
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
                    if v_ind in coms and i in coms[v_ind]:
                        arg = list(set(delta_p_act[i][v_ind]) - set(argmin[i][v_ind]))
                        arg_nf = list(set(arg) - set(bounds_f[v_ind][i]))
                        if not arg_nf:
                            # In diesem Fall sind alle aktiven nichtminimalen Kanten bereits fixiert. Damit können keine Fortschritte mehr gemacht werden -> relaxiere bds
                            for e_ind in arg:
                                if xk[i, e_ind] < self.bd_tol:
                                    continue
                                w_ind = self.V.index(self.E[e_ind][1])
                                if self.q_global[-1][e_ind] < self.eps:
                                    # Falls q = 0 -> g_e nicht-negativ -> kleinster a-Wert kann nicht erreicht werden -> wähle kleinstmögliche untere Schranke: 0
                                    bounds[v_ind][i][e_ind] = (0, xk[i, e_ind])
                                else:
                                    # relaxiere untere Schranke, sodass der minimale a-Wert erreicht werden kann
                                    bounds[v_ind][i][e_ind] = ((1 + a[i, v_ind] - a[i, w_ind]) * self.nu[e_ind], xk[i, e_ind])
                                bounds_f[v_ind][i].remove(e_ind)
                                delta_p_act_nf[i][v_ind].append(e_ind)
                                b_nf[i, v_ind] += xk[i, e_ind]
                        if not argmin_nf[i][v_ind]:
                            # Dagegen sind in diesem Fall alle minimalen Kanten fixiert. Damit können ebenfalls keine Forschritte gemacht werden -> relaxiere bds
                            for e_ind in argmin[i][v_ind]:
                                if xk[i, e_ind] > self.b[-1][i, v_ind] - self.bd_tol:
                                    continue
                                w_ind = self.V.index(self.E[e_ind][1])
                                xe = x_sum[e_ind] + x_sum_fix[e_ind]
                                x2 = (a_min2[e_ind] - a[i, w_ind]) * self.nu[e_ind] - self.g(e_ind, theta_ind, xe)
                                if self.q_global[-1][e_ind] < self.eps and xe < self.nu[e_ind] - self.eps:
                                    x2 += self.nu[e_ind] - xe
                                bounds_f[v_ind][i].remove(e_ind)
                                delta_p_act_nf[i][v_ind].append(e_ind)
                                # relaxierte obere Schranke mit x2: Dies ist der Flusswert, sodass der zweitkleinste a-Wert (= 'a_min2[e_ind]') auch für Kante 'e_ind' erreicht
                                # werden kann.
                                bounds[v_ind][i][e_ind] = (xk[i, e_ind], np.min([x2, self.b[-1][i, v_ind]]))
                                b_nf[i, v_ind] += xk[i, e_ind]

            coms_keys = list(coms)
            for v_ind in coms_keys:
                # Prüfe, ob alle aktiven Kanten den gleichen a-Wert liefern
                # ??? Eigentlich nur für Güter mit nichtdisjunkten Mengen delta_p_act_nf[i][v_ind] ?
                if np.all([v_ind not in argmin_nf[i] or len(argmin_nf[i][v_ind]) == len(delta_p_act_nf[i][v_ind]) for i in coms[v_ind]]):
                    outneighbors = list(set([self.V.index(self.E[e_ind][1]) for e_ind in argmin_nf[i][v_ind] for i in coms[v_ind]]))
                    # sind zusätzlich für alle Endknoten 'w_ind' die Flussaufteilungen bereits fixiert, kann auch die Aufteilung in Knoten 'v_ind' fixiert werden
                    if np.all([w_ind not in coms for w_ind in outneighbors]):
                        for i in coms[v_ind]:
                            for e_ind in delta_p_act_nf[i][v_ind]:
                                if xk[i, e_ind] < self.bd_tol:
                                    xfp[i, e_ind] = 0
                                else:
                                    xfp[i, e_ind] = xk[i, e_ind]
                                    x_sum_fix[e_ind] += xfp[i, e_ind]
                                x_sum[e_ind] -= xk[i, e_ind]
                                xk[i, e_ind] = 0
                            delta_p_act_nf[i][v_ind] = []
                        del coms[v_ind]
                        if not coms:  # Abbruchbedingung: keine Güter mehr übrig, die noch auf Kanten aufgeteilt werden müssen
                            return xfp, x_sum_fix, a, top_ords

            while True:
                print("theta_ind", theta_ind, "thta", theta)
                gx, forced_zeros = self.gamma(xk + xfp, x_sum + x_sum_fix, a, delta_p_act_nf, coms, theta_ind)
                # lil_matrix wird beim Rechnen automatisch in csr Matrix umgewandelt!?
                # if scipy.sparse.linalg.norm(gx.tocsr() - xk.tocsr()) < self.bd_tol:
                if scipy.sparse.linalg.norm(gx - xk) < self.bd_tol:
                    sol = xfp + xk
                    sparse_items = self.get_items(xk)
                    for (i, e_ind) in sparse_items:
                        if xk[i, e_ind] < self.bd_tol:
                            x_sum[e_ind] -= xk[i, e_ind]
                            sol[i, e_ind] -= xk[i, e_ind]
                    return sol, x_sum + x_sum_fix, a, top_ords

                fp_comp = []
                for v_ind in coms:
                    for i in coms[v_ind]:
                        fz_nf = list(set(forced_zeros[i][v_ind]) - set(bounds_f[v_ind][i]))
                        am_nf = list(set(argmin[i][v_ind]) - set(bounds_f[v_ind][i]))
                        for e_ind in fz_nf:
                            if e_ind in bounds[v_ind][i]:
                                bounds[v_ind][i][e_ind] = (bounds[v_ind][i][e_ind][0], xk[i, e_ind])
                            else:
                                bounds[v_ind][i][e_ind] = (0, xk[i, e_ind])
                            bd_diff = abs(bounds[v_ind][i][e_ind][0] - bounds[v_ind][i][e_ind][1])
                            if bd_diff < self.bd_tol:
                                if bd_diff == 0 and xk[i, e_ind] > 0:
                                    # bd wird relaxiert
                                    bounds[v_ind][i][e_ind] = (0, xk[i, e_ind])
                                else:
                                    bounds_f[v_ind][i].append(e_ind)
                                    b_nf[i, v_ind] -= xk[i, e_ind]
                                    delta_p_act_nf[i][v_ind].remove(e_ind)
                                    if b_nf[i, v_ind] < self.eps:
                                        fp_comp.append((i, v_ind))
                        if v_ind not in argmin_nf[i] or len(argmin_nf[i][v_ind]) == len(delta_p_act_nf[i][v_ind]):
                            continue
                        for e_ind in am_nf:
                            if e_ind in bounds[v_ind][i]:
                                bounds[v_ind][i][e_ind] = (xk[i, e_ind], bounds[v_ind][i][e_ind][1])
                            else:
                                bounds[v_ind][i][e_ind] = (xk[i, e_ind], self.b[-1][i, v_ind])
                            bd_diff = abs(bounds[v_ind][i][e_ind][0] - bounds[v_ind][i][e_ind][1])
                            if bd_diff < self.bd_tol:
                                if bd_diff == 0 and xk[i, e_ind] < self.b[-1][i, v_ind]:
                                    # bd wird relaxiert
                                    bounds[v_ind][i][e_ind] = (xk[i, e_ind], self.b[-1][i, v_ind])
                                else:
                                    bounds_f[v_ind][i].append(e_ind)
                                    b_nf[i, v_ind] -= xk[i, e_ind]
                                    delta_p_act_nf[i][v_ind].remove(e_ind)
                                    if b_nf[i, v_ind] < self.eps:
                                        fp_comp.append((i, v_ind))

                '''for v_ind in coms:
                    # wird kein Gut in Knoten 'v_ind' auf mehrere Kanten aufgeteilt, so können 'bounds' angepasst werden
                    if np.all([v_ind not in argmin[i] or len(argmin[i][v_ind]) <= 1 for i in coms[v_ind]]):
                        """coms_skip = []
                        coms_left = copy.deepcopy(coms[v_ind])"""
                        for i in coms[v_ind]:
                            """if i in coms_skip:
                                continue"""
                            try:
                                e_ind = argmin[i][v_ind][0]
                            except KeyError:
                                continue
                            if e_ind in bounds[v_ind][i]:
                                bd = bounds[v_ind][i][e_ind]
                                if bd[1] - xk[i, e_ind] - xfp[i, e_ind] > self.eps:
                                    # 'e_ind' ist einzige Kante in 'argmin[i][v_ind]', also kann Gesamtflusswert als untere Schranke verwendet werden
                                    bounds[v_ind][i][e_ind] = (xk[i, e_ind] + xfp[i, e_ind], bd[1])
                                    if abs(bounds[v_ind][i][e_ind][0] - bounds[v_ind][i][e_ind][1]) < self.bd_tol:
                                        # liegen Schranken nahe genug beieinander, besteht die Möglichkeit, dass der Fluss für Kante 'e_ind' fixiert werden kann. Andererseits
                                        # können die berechneten Schranken auch falsch sein, da sie nur für approximierte a-Werte berechnet wurden (nicht für die a-Werte
                                        # korrespondierend zu einer optimalen Flussaufteilung). Welcher der beiden Fälle hier vorliegt, wird im Folgenden geprüft.
                                        mean_a = 0
                                        for f_ind in delta_p_act_nf[i][v_ind]:
                                            w_ind = self.V.index(self.E[f_ind][1])
                                            mean_a += self.g(f_ind, theta_ind, x_sum[f_ind] + x_sum_fix[f_ind]) / self.nu[f_ind] + a[i, w_ind]
                                        mean_a /= len(delta_p_act_nf[i][v_ind])
                                        e2 = self.V.index(self.E[e_ind][1])
                                        a_e = self.g(e_ind, theta_ind, x_sum[e_ind] + x_sum_fix[e_ind]) / self.nu[e_ind] + a[i, e2]
                                        if a_e > mean_a + self.eps:
                                            # a-Wert unrealistisch -> verwerfe Schranken
                                            del bounds[v_ind][i][e_ind]
                                            break
                                        else:
                                            # a-Wert realistisch -> fixiere Gut 'i' für 'e_ind'
                                            fp_comp.append((i, e_ind))
                                        """e_min = []
                                        i_ind = coms_left.index(i)
                                        for j in coms_left[i_ind:]:
                                            if v_ind in argmin[j] and e_ind == argmin[j][v_ind][0]:
                                                e_min.append(j)
                                        for j in e_min:
                                            mean_a = 0
                                            for f_ind in delta_p_act_nf[j][v_ind]:
                                                w_ind = self.V.index(self.E[f_ind][1])
                                                mean_a += self.g(f_ind, theta_ind, x_sum[f_ind] + x_sum_fix[f_ind]) / self.nu[f_ind] + a[j, w_ind]
                                            mean_a /= len(delta_p_act_nf[j][v_ind])
                                            e2 = self.V.index(self.E[e_ind][1])
                                            a_e = self.g(e_ind, theta_ind, x_sum[e_ind] + x_sum_fix[e_ind]) / self.nu[e_ind] + a[j, e2]
                                            if a_e > mean_a + self.eps:
                                                # a-Wert unrealistisch -> verwerfe Schranken
                                                del bounds[v_ind][e_ind]
                                                break
                                            elif j == e_min[-1]:
                                                # a-Wert realistisch -> fixiere 'e_ind'
                                                fp_comp.append(e_ind)
                                        for j in e_min:
                                            coms_skip.append(j)
                                            coms_left.remove(j)"""

                            else:
                                bounds[v_ind][i][e_ind] = (xk[i, e_ind] + xfp[i, e_ind], b[i, v_ind])
                                if abs(bounds[v_ind][i][e_ind][0] - bounds[v_ind][i][e_ind][1]) < self.bd_tol:
                                    # liegen Schranken nahe genug beieinander, besteht die Möglichkeit, dass der Fluss für Kante 'e_ind' fixiert werden kann. Andererseits
                                    # können die berechneten Schranken auch falsch sein, da sie nur für approximierte a-Werte berechnet wurden (nicht für die a-Werte
                                    # korrespondierend zu einer optimalen Flussaufteilung). Welcher der beiden Fälle hier vorliegt, wird im Folgenden geprüft.
                                    mean_a = 0
                                    for f_ind in delta_p_act_nf[i][v_ind]:
                                        w_ind = self.V.index(self.E[f_ind][1])
                                        mean_a += self.g(f_ind, theta_ind, x_sum[f_ind] + x_sum_fix[f_ind]) / self.nu[f_ind] + a[i, w_ind]
                                    mean_a /= len(delta_p_act_nf[i][v_ind])
                                    e2 = self.V.index(self.E[e_ind][1])
                                    a_e = self.g(e_ind, theta_ind, x_sum[e_ind] + x_sum_fix[e_ind]) / self.nu[e_ind] + a[i, e2]
                                    if a_e > mean_a + self.eps:
                                        # a-Wert unrealistisch -> verwerfe Schranken
                                        del bounds[v_ind][i][e_ind]
                                        break
                                    else:
                                        # a-Wert realistisch -> fixiere Gut 'i' für 'e_ind'
                                        fp_comp.append((i, e_ind))
                                    """e_min = []
                                    i_ind = coms_cpy.index(i)
                                    for j in coms_cpy[i_ind:]:
                                        if v_ind in argmin[j] and e_ind == argmin[j][v_ind][0]:
                                            e_min.append(j)
                                    for j in e_min:
                                        mean_a = 0
                                        for f_ind in delta_p_act_nf[j][v_ind]:
                                            w_ind = self.V.index(self.E[f_ind][1])
                                            mean_a += self.g(f_ind, theta_ind, x_sum[f_ind] + x_sum_fix[f_ind]) / self.nu[f_ind] + a[j, w_ind]
                                        mean_a /= len(delta_p_act_nf[j][v_ind])
                                        e2 = self.V.index(self.E[e_ind][1])
                                        a_e = self.g(e_ind, theta_ind, x_sum[e_ind] + x_sum_fix[e_ind]) / self.nu[e_ind] + a[j, e2]
                                        if a_e > mean_a + self.eps:
                                            # a-Wert unrealistisch -> verwerfe Schranken
                                            del bounds[v_ind][e_ind]
                                            break
                                        elif j == e_min[-1]:
                                            # a-Wert realistisch -> fixiere 'e_ind'
                                            fp_comp.append(e_ind)
                                    for j in e_min:
                                        coms_cpy.remove(j)"""
                    elif np.all(np.array([len(forced_zeros[i][v_ind]) <= 1 for i in coms[v_ind]])):
                        for i in coms[v_ind]:
                            try:
                                e_ind = forced_zeros[i][v_ind][0]
                            except IndexError:
                                continue
                            if e_ind in bounds[v_ind][i]:
                                bd = bounds[v_ind][i][e_ind]
                                if xk[i, e_ind] + xfp[i, e_ind] - bd[0] > self.eps:
                                    bounds[v_ind][i][e_ind] = (bd[0], xk[i, e_ind] + xfp[i, e_ind])
                                    if abs(bounds[v_ind][i][e_ind][0] - bounds[v_ind][i][e_ind][1]) < self.bd_tol:
                                        # liegen Schranken nahe genug beieinander, besteht die Möglichkeit, dass der Fluss für Kante 'e_ind' fixiert werden kann. Andererseits
                                        # können die berechneten Schranken auch falsch sein, da sie nur für approximierte a-Werte berechnet wurden (nicht für die a-Werte
                                        # korrespondierend zu einer optimalen Flussaufteilung). Welcher der beiden Fälle hier vorliegt, wird im Folgenden geprüft.
                                        mean_a = 0
                                        for f_ind in delta_p_act_nf[i][v_ind]:
                                            w_ind = self.V.index(self.E[f_ind][1])
                                            mean_a += self.g(f_ind, theta_ind, x_sum[f_ind] + x_sum_fix[f_ind]) / self.nu[f_ind] + a[i, w_ind]
                                        mean_a /= len(delta_p_act_nf[i][v_ind])
                                        e2 = self.V.index(self.E[e_ind][1])
                                        a_e = self.g(e_ind, theta_ind, x_sum[e_ind] + x_sum_fix[e_ind]) / self.nu[e_ind] + a[i, e2]
                                        if a_e > mean_a + self.eps:
                                            # a-Wert unrealistisch -> verwerfe Schranken
                                            del bounds[v_ind][i][e_ind]
                                            break
                                        else:
                                            # a-Wert realistisch -> fixiere Gut 'i' für 'e_ind'
                                            fp_comp.append((i, e_ind))
                                        """e_fc = []
                                        i_ind = coms_cpy.index(i)
                                        for j in coms_cpy[i_ind:]:
                                            if v_ind in forced_zeros[j] and e_ind == forced_zeros[j][v_ind][0]:
                                                e_fc.append(j)
                                        for j in e_fc:
                                            mean_a = 0
                                            for f_ind in delta_p_act_nf[j][v_ind]:
                                                w_ind = self.V.index(self.E[f_ind][1])
                                                mean_a += self.g(f_ind, theta_ind, x_sum[f_ind] + x_sum_fix[f_ind]) / self.nu[f_ind] + a[j, w_ind]
                                            mean_a /= len(delta_p_act_nf[j][v_ind])
                                            e2 = self.V.index(self.E[e_ind][1])
                                            a_e = self.g(e_ind, theta_ind, x_sum[e_ind] + x_sum_fix[e_ind]) / self.nu[e_ind] + a[j, e2]
                                            if a_e > mean_a + self.eps:
                                                # a-Wert unrealistisch -> verwerfe Schranken
                                                del bounds[v_ind][e_ind]
                                                break
                                            elif j == e_fc[-1]:
                                                # a-Wert realistisch -> fixiere 'e_ind'
                                                fp_comp.append(e_ind)
                                        for j in e_fc:
                                            coms_cpy.remove(j)"""
                            else:
                                bounds[v_ind][i][e_ind] = (0, xk[i, e_ind] + xfp[i, e_ind])
                                if abs(bounds[v_ind][i][e_ind][0] - bounds[v_ind][i][e_ind][1]) < self.bd_tol:
                                    # liegen Schranken nahe genug beieinander, besteht die Möglichkeit, dass der Fluss für Kante 'e_ind' fixiert werden kann. Andererseits
                                    # können die berechneten Schranken auch falsch sein, da sie nur für approximierte a-Werte berechnet wurden (nicht für die a-Werte
                                    # korrespondierend zu einer optimalen Flussaufteilung). Welcher der beiden Fälle hier vorliegt, wird im Folgenden geprüft.
                                    mean_a = 0
                                    for f_ind in delta_p_act_nf[i][v_ind]:
                                        w_ind = self.V.index(self.E[f_ind][1])
                                        mean_a += self.g(f_ind, theta_ind, x_sum[f_ind] + x_sum_fix[f_ind]) / self.nu[f_ind] + a[i, w_ind]
                                    mean_a /= len(delta_p_act_nf[i][v_ind])
                                    e2 = self.V.index(self.E[e_ind][1])
                                    a_e = self.g(e_ind, theta_ind, x_sum[e_ind] + x_sum_fix[e_ind]) / self.nu[e_ind] + a[i, e2]
                                    if a_e > mean_a + self.eps:
                                        # a-Wert unrealistisch -> verwerfe Schranken
                                        del bounds[v_ind][i][e_ind]
                                        break
                                    else:
                                        # a-Wert realistisch -> fixiere Gut 'i' für 'e_ind'
                                        fp_comp.append((i, e_ind))
                                    """e_fc = []
                                    i_ind = coms_cpy.index(i)
                                    for j in coms_cpy[i_ind:]:
                                        if v_ind in forced_zeros[j] and e_ind == forced_zeros[j][v_ind][0]:
                                            e_fc.append(j)
                                    for j in e_fc:
                                        mean_a = 0
                                        for f_ind in delta_p_act_nf[j][v_ind]:
                                            w_ind = self.V.index(self.E[f_ind][1])
                                            mean_a += self.g(f_ind, theta_ind, x_sum[f_ind] + x_sum_fix[f_ind]) / self.nu[f_ind] + a[j, w_ind]
                                        mean_a /= len(delta_p_act_nf[j][v_ind])
                                        e2 = self.V.index(self.E[e_ind][1])
                                        a_e = self.g(e_ind, theta_ind, x_sum[e_ind] + x_sum_fix[e_ind]) / self.nu[e_ind] + a[j, e2]
                                        if a_e > mean_a + self.eps:
                                            # a-Wert unrealistisch -> verwerfe Schranken
                                            del bounds[v_ind][e_ind]
                                            break
                                        elif j == e_fc[-1]:
                                            # a-Wert realistisch -> fixiere 'e_ind'
                                            fp_comp.append(e_ind)
                                    for j in e_fc:
                                        coms_cpy.remove(j)"""'''

                if fp_comp:
                    # es existieren Knoten in 'fp_comp', welche fixiert werden können
                    break

                xk, x_sum = calc_flow_by_bounds(xk_old=xk)

                # berechne zugehörige a-Werte
                a = scipy.sparse.lil_matrix((self.I, self.n))
                sparse_items = []
                for i in range(self.I):
                    t_ind = self.V.index('t{}'.format(i+1))
                    a[i, t_ind] = 10**(-14)  # < 'self.eps', wird also als 0 behandelt
                    sparse_items.append((i, t_ind))
                argmin = []
                argmin_nf = []
                for i in range(self.I):
                    argmin.append({})
                    argmin_nf.append({})
                    for v in top_ords[i]:
                        v_ind = self.V.index(v)
                        a_min2 = {}
                        for e_ind in delta_p_act[i][v_ind]:
                            w_ind = self.V.index(self.E[e_ind][1])
                            a_e = self.g(e_ind, theta_ind, x_sum[e_ind] + x_sum_fix[e_ind]) / self.nu[e_ind] + a[i, w_ind]
                            if (i, v_ind) not in sparse_items:
                                a[i, v_ind] = a_e
                                sparse_items.append((i, v_ind))
                                argmin[i][v_ind] = [e_ind]
                                if e_ind in delta_p_act_nf[i][v_ind]:
                                    argmin_nf[i][v_ind] = [e_ind]
                                else:
                                    # min wird in bereits fixierter Kante angenommen. 'argmin_nf' enthält jedoch nur nicht-fixierte, minimale Kanten
                                    argmin_nf[i][v_ind] = []
                            elif a_e < a[i, v_ind] + self.bd_tol:
                                if a_e < a[i, v_ind] - self.bd_tol:
                                    a_min2[e_ind] = a[i, v_ind]
                                    a[i, v_ind] = a_e
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
                        if v_ind in coms and i in coms[v_ind]:
                            arg = list(set(delta_p_act[i][v_ind]) - set(argmin[i][v_ind]))
                            arg_nf = list(set(arg) - set(bounds_f[v_ind][i]))
                            if not arg_nf:
                                # In diesem Fall sind alle aktiven nichtminimalen Kanten bereits fixiert. Damit können keine Fortschritte mehr gemacht werden -> relaxiere bds
                                for e_ind in arg:
                                    if xk[i, e_ind] < self.bd_tol:
                                        continue
                                    w_ind = self.V.index(self.E[e_ind][1])
                                    if self.q_global[-1][e_ind] < self.eps:
                                        # Falls q = 0 -> g_e nicht-negativ -> kleinster a-Wert kann nicht erreicht werden -> wähle kleinstmögliche untere Schranke: 0
                                        bounds[v_ind][i][e_ind] = (0, xk[i, e_ind])
                                    else:
                                        # relaxiere untere Schranke, sodass der minimale a-Wert erreicht werden kann
                                        bounds[v_ind][i][e_ind] = ((1 + a[i, v_ind] - a[i, w_ind]) * self.nu[e_ind], xk[i, e_ind])
                                    bounds_f[v_ind][i].remove(e_ind)
                                    delta_p_act_nf[i][v_ind].append(e_ind)
                                    b_nf[i, v_ind] += xk[i, e_ind]
                            if not argmin_nf[i][v_ind]:
                                # Dagegen sind in diesem Fall alle minimalen Kanten fixiert. Damit können ebenfalls keine Forschritte gemacht werden -> relaxiere bds
                                for e_ind in argmin[i][v_ind]:
                                    if xk[i, e_ind] > self.b[-1][i, v_ind] - self.bd_tol:
                                        continue
                                    w_ind = self.V.index(self.E[e_ind][1])
                                    xe = x_sum[e_ind] + x_sum_fix[e_ind]
                                    x2 = (a_min2[e_ind] - a[i, w_ind]) * self.nu[e_ind] - self.g(e_ind, theta_ind, xe)
                                    if self.q_global[-1][e_ind] < self.eps and xe < self.nu[e_ind] - self.eps:
                                        x2 += self.nu[e_ind] - xe
                                    bounds_f[v_ind][i].remove(e_ind)
                                    delta_p_act_nf[i][v_ind].append(e_ind)
                                    # relaxierte obere Schranke mit x2: Dies ist der Flusswert, sodass der zweitkleinste a-Wert (= 'a_min2[e_ind]') auch für Kante 'e_ind' erreicht
                                    # werden kann.
                                    bounds[v_ind][i][e_ind] = (xk[i, e_ind], np.min([x2, self.b[-1][i, v_ind]]))
                                    b_nf[i, v_ind] += xk[i, e_ind]

                """for v_ind in coms:
                    if not delta_z[v_ind]:
                        bounds_mid_sum = 0
                        for e_ind in e_diff[v_ind]:
                            bounds_mid_sum += 0.5 * (bounds[v_ind][e_ind][0] + bounds[v_ind][e_ind][1])
                        if bounds_mid_sum > flow_vol_nf[v_ind]:
                            max_e = -np.Inf
                            max_av_inds = []
                            for e_ind in e_diff[v_ind]:
                                w_ind = self.V.index(self.E[e_ind][1])
                                a_e = self.g(e_ind, theta_ind, x_sum[e_ind]) / self.nu[e_ind] + sum([a[i, w_ind] for i in coms[v_ind] if self.E_active[i][-1][e_ind]])
                                if a_e > max_e + self.eps:
                                    max_e = a_e
                                    max_av_inds = [e_ind]
                                elif a_e > max_e - self.eps:
                                    max_av_inds.append(e_ind)
                            for e_ind in max_av_inds:
                                lb = bounds[v_ind][e_ind][0]
                                bounds[v_ind][e_ind] = (lb, 0.5 * (bounds[v_ind][e_ind][0] + bounds[v_ind][e_ind][1]))
                                if abs(bounds[v_ind][e_ind][0] - bounds[v_ind][e_ind][1]) < self.bd_tol:
                                    fp_comp.append(e_ind)
                if fp_comp:
                    break"""
                """for v_ind in coms:
                    if not delta_z[v_ind]:
                        av_max = {}
                        av_min = {}
                        sparse_items = get_items(a)
                        for e_ind in e_diff[v_ind]:
                            w_ind = self.V.index(self.E[e_ind][1])
                            av_max[e_ind] = self.g(e_ind, theta_ind, x_sum[e_ind]) / self.nu[e_ind] + max([a[i, w_ind] for i in coms[v_ind]])
                            av_min[e_ind] = self.g(e_ind, theta_ind, x_sum[e_ind]) / self.nu[e_ind] + min([a[i, w_ind] for i in coms[v_ind] if (i, w_ind) in sparse_items])

                        # enthält Kanten, die für den maximalen a-Wert verantwortlich sind
                        max_av_inds = []
                        # enthält Kanten, die für den minimalen a-Wert verantwortlich sind
                        min_av_inds = []
                        max_e = max(av_max.values())
                        min_e = min(av_min.values())
                        for e_ind in e_diff[v_ind]:
                            if abs(max_e - av_max[e_ind]) < self.eps:
                                max_av_inds.append(e_ind)
                            elif abs(min_e - av_min[e_ind]) < self.eps:
                                min_av_inds.append(e_ind)
                        for e_ind in max_av_inds:
                            int_c = 0.5 * (bounds[v_ind][e_ind][0] + bounds[v_ind][e_ind][1])
                            # 'int_c': erster approximierter Gesamtflusswert für 'e_ind'
                            # 'x_sum[e_ind]': Gesamtflusswert nach Korrektur
                            if x_sum[e_ind] < int_c:
                                # 'int_c' war zu groß, da selbst für kleineren Wert 'x_sum[e_ind]' der maximale a-Wert erzeugt wird
                                bounds[v_ind][e_ind] = (bounds[v_ind][e_ind][0], int_c)
                            if abs(bounds[v_ind][e_ind][0] - bounds[v_ind][e_ind][1]) < self.bd_tol:
                                # liegen Schranken nahe genug beieinander, besteht die Möglichkeit, dass der Fluss für Kante 'e_ind' fixiert werden kann. Andererseits
                                # können die berechneten Schranken auch falsch sein, da sie nur für approximierte a-Werte berechnet wurden (nicht für die a-Werte
                                # korrespondierend zu einer optimalen Flussaufteilung). Welcher der beiden Fälle hier vorliegt, wird im Folgenden geprüft.
                                for j in coms[v_ind]:
                                    if e_ind not in delta_p_act_nf[j][v_ind]:
                                        continue
                                    mean_a = 0
                                    for f_ind in delta_p_act_nf[j][v_ind]:
                                        w_ind = self.V.index(self.E[f_ind][1])
                                        mean_a += self.g(f_ind, theta_ind, x_sum[f_ind] + x_sum_fix[f_ind]) / self.nu[f_ind] + a[j, w_ind]
                                    mean_a /= len(delta_p_act_nf[j][v_ind])
                                    e2 = self.V.index(self.E[e_ind][1])
                                    a_e = self.g(e_ind, theta_ind, x_sum[e_ind] + x_sum_fix[e_ind]) / self.nu[e_ind] + a[j, e2]
                                    if a_e > mean_a + self.eps:
                                        # a-Wert unrealistisch -> verwerfe Schranken
                                        del bounds[v_ind][e_ind]
                                        break
                                if e_ind in bounds[v_ind]:
                                    # a-Wert realistisch -> fixiere 'e_ind'
                                    fp_comp.append(e_ind)
                        for e_ind in min_av_inds:
                            int_c = 0.5 * (bounds[v_ind][e_ind][0] + bounds[v_ind][e_ind][1])
                            # 'int_c': erster approximierter Gesamtflusswert für 'e_ind'
                            # 'x_sum[e_ind]': Gesamtflusswert nach Korrektur
                            if x_sum[e_ind] > int_c:
                                # 'int_c' war zu klein, da selbst für größeren Wert 'x_sum[e_ind]' der minimale a-Wert erzeugt wird
                                bounds[v_ind][e_ind] = (int_c, bounds[v_ind][e_ind][1])
                            if abs(bounds[v_ind][e_ind][0] - bounds[v_ind][e_ind][1]) < self.bd_tol:
                                # liegen Schranken nahe genug beieinander, besteht die Möglichkeit, dass der Fluss für Kante 'e_ind' fixiert werden kann. Andererseits
                                # können die berechneten Schranken auch falsch sein, da sie nur für approximierte a-Werte berechnet wurden (nicht für die a-Werte
                                # korrespondierend zu einer optimalen Flussaufteilung). Welcher der beiden Fälle hier vorliegt, wird im Folgenden geprüft.
                                for j in coms[v_ind]:
                                    if e_ind not in delta_p_act_nf[j][v_ind]:
                                        continue
                                    mean_a = 0
                                    for f_ind in delta_p_act_nf[j][v_ind]:
                                        w_ind = self.V.index(self.E[f_ind][1])
                                        mean_a += self.g(f_ind, theta_ind, x_sum[f_ind] + x_sum_fix[f_ind]) / self.nu[f_ind] + a[j, w_ind]
                                    mean_a /= len(delta_p_act_nf[j][v_ind])
                                    e2 = self.V.index(self.E[e_ind][1])
                                    a_e = self.g(e_ind, theta_ind, x_sum[e_ind] + x_sum_fix[e_ind]) / self.nu[e_ind] + a[j, e2]
                                    if a_e > mean_a + self.eps:
                                        # a-Wert unrealistisch -> verwerfe Schranken
                                        del bounds[v_ind][e_ind]
                                        break
                                if e_ind in bounds[v_ind]:
                                    # a-Wert realistisch -> fixiere 'e_ind'
                                    fp_comp.append(e_ind)"""

                """coms_keys = list(coms)
                for v_ind in coms_keys:
                    # Prüfe, ob alle aktiven Kanten den gleichen a-Wert liefern
                    if np.all([v_ind not in argmin[i] or len(argmin[i][v_ind]) == len(delta_p_act_nf[i][v_ind]) for i in coms[v_ind]]):
                        outneighbors = list(set([self.V.index(self.E[e_ind][1]) for e_ind in argmin[i][v_ind] for i in coms[v_ind]]))
                        # sind zusätzlich für alle Endknoten 'w_ind' die Flussaufteilungen bereits fixiert, kann auch die Aufteilung in Knoten 'v_ind' fixiert werden
                        if np.all([w_ind not in coms for w_ind in outneighbors]):
                            for i in coms[v_ind]:
                                for e_ind in delta_p_act_nf[i][v_ind]:
                                    if xk[i, e_ind] < self.bd_tol:
                                        xfp[i, e_ind] = 0
                                    else:
                                        xfp[i, e_ind] = xk[i, e_ind]
                                    x_sum_fix[e_ind] += xfp[i, e_ind]
                                    x_sum[e_ind] -= xk[i, e_ind]
                                    xk[i, e_ind] = 0
                                delta_p_act_nf[i][v_ind] = []
                            del coms[v_ind]
                            if not coms:
                                # Abbruchbedingung: keine Güter mehr übrig, die noch auf Kanten aufgeteilt werden müssen
                                # return xfp, x_sum + x_sum_fix
                                return xfp, x_sum_fix, top_ords"""

                coms_keys = list(coms)
                for v_ind in coms_keys:
                    if np.all([v_ind not in argmin_nf[i] or len(argmin_nf[i][v_ind]) == len(delta_p_act_nf[i][v_ind]) for i in coms[v_ind]]):
                        outneighbors = list(set([self.V.index(self.E[e_ind][1]) for e_ind in argmin[i][v_ind] for i in coms[v_ind]]))
                        # sind zusätzlich für alle Endknoten 'w_ind' die Flussaufteilungen bereits fixiert, kann auch die Aufteilung in Knoten 'v_ind' fixiert werden
                        if np.all([w_ind not in coms for w_ind in outneighbors]):
                            for i in coms[v_ind]:
                                for e_ind in delta_p_act_nf[i][v_ind]:
                                    if xk[i, e_ind] < self.bd_tol:
                                        xfp[i, e_ind] = 0
                                    else:
                                        xfp[i, e_ind] = xk[i, e_ind]
                                        x_sum_fix[e_ind] += xfp[i, e_ind]
                                    x_sum[e_ind] -= xk[i, e_ind]
                                    xk[i, e_ind] = 0
                                delta_p_act_nf[i][v_ind] = []
                            del coms[v_ind]
                            if not coms:
                                # Abbruchbedingung: keine Güter mehr übrig, die noch auf Kanten aufgeteilt werden müssen
                                return xfp, x_sum_fix, a, top_ords

    def main(self):
        """
        Hauptteil: Konstruiert schrittweise einen kontinuierlichen IDE-Fluss zu den Eingabedaten aus 'Math.cont_data'. Erzeugt Tabelle der
        Klasse 'Graphics.OutputTable' mit den entsprechenden Daten des Flusses.
        :return: 0
        """
        theta = 0
        theta_ind = -1
        # Obergrenze für theta
        T = 14
        stop_outflow = []
        # simplex = [[], [], []]
        # all_iv = {}
        while theta < T:
            # in der Zukunft liegende Zeitpunkte aus der Liste 'self.u_start'
            start_points = [t for t in self.u_start if t > theta]
            # in der Zukunft liegende Zeitpunkte, zu denen der f^- -Wert von mindestens einer Kante auf 0 springt (wird
            # während der Laufzeit aktualisiert)
            stop_outflow = [t for t in stop_outflow if t > theta]
            theta_ind += 1
            # Liste aller Kanten, deren Warteschlange in der aktuellen Phase 0 wird, und deren f^- -Werte auf 0
            # gesetzt werden müssen
            fm_to_zero = []
            # 'next_phase' bestimmt am Ende der Schleife, wie lange aktuelle Phase andauert
            if start_points or stop_outflow:
                next_phase = np.min(start_points + stop_outflow) - theta
            else:
                next_phase = T
            x_total, x_sum, a, top_ords = self.fp_approx(theta_ind)
            a_items = self.get_items(a)
            for (i, v_ind) in a_items:
                if abs(a[i, v_ind]) < self.eps:
                    a[i, v_ind] = 0
            """
            self.flow_vol.append(np.zeros(self.n))
            # Anz. Güter je Knoten, die aktuell eine Entscheidung über Flussaufteilung treffen müssen
            com_no = np.zeros(self.n)
            for v_ind in range(self.n):
                v = self.V[v_ind]
                for i in range(self.I):
                    # speichere b - Wert
                    b[i, v_ind] = self.calc_b(i, v, theta)
                    if b[i, v_ind] > self.eps:
                        com_no[v_ind] += 1
                        self.flow_vol[-1][v_ind] += b[i, v_ind]"""
            """simplex = self.find_cover_simplex(theta_ind)
            var_no = len(simplex[0])
            rows, col_pointer = simplex[1:]
            A_1 = np.zeros((self.I, var_no))
            # A_1 = np.zeros((self.I * (self.n-1), self.I * self.m))
            # A_2 = np.zeros((self.I, self.I * self.m))

            b_1 = np.zeros(self.I * (self.n-1))

            co = 0
            for ro in range(var_no):
                i = rows[ro]
                while col_pointer[co + 1] <= ro:
                    co += 1
                e = self.E[co]
                v = e[0]
                v_ind = self.V.index(v)
                A_1[i, ro] = 1
            b_1 = np.ones(A_1.shape[0])"""

            '''
            for i in range(self.I):
                ti_ind = self.V.index('t{}'.format(i+1))
                for v_ind in range(self.n):
                    if v_ind < ti_ind:
                        vi_ind = v_ind
                    elif v_ind == ti_ind:
                        continue
                    else:
                        vi_ind = v_ind - 1
                    # b_1[i * (self.n-1) + vi_ind] = b[i, v_ind]
                    # nicht b_{i,v}, sondern seine bary. Koordinaten, da wir in Standardsimplex arbeiten
                    if b[i, v_ind] > self.eps:
                        b_1[i * (self.n-1) + vi_ind] = 1
                    else:
                        b_1[i * (self.n-1) + vi_ind] = 0
                    delta_v = self.get_outgoing_edges(self.V[v_ind])  # Liste der Kantenindizes
                    for e_ind in delta_v:
                        # A_1[i * (self.n-1) + vi_ind, i * self.m + e_ind] = 1
                        A_1[i * (self.n-1) + vi_ind, self.I * e_ind + i] = 1

                for v in self.V:
                    delta_v_inact = list(set(self.get_outgoing_edges(v)) - set(self.get_outgoing_active_edges(i, v)))
                    for e_ind in delta_v_inact:
                        # A_2[i, i * self.m + e_ind] = 1
                        A_2[i, self.I * e_ind + i] = 1'''

            '''col_changes = []
            zero_cols = []
            for i in range(self.I):
                for cp in range(1, self.m + 1):
                    if simplex[2][cp] != simplex[2][cp-1]:
                        col_changes.append(cp - 1 + i * self.m)
            for i in range(self.I):
                for v_ind in range(self.n):
                    if b[i, v_ind] < self.eps:
                        delta_p = self.get_outgoing_edges(self.V[v_ind])
                        for e_ind in delta_p:
                            zero_cols.append(i * self.m + e_ind)
            col_changes_cpy = col_changes.copy()
            col_changes = [c for c in col_changes_cpy if c not in zero_cols]
            A_1 = A_1[:, col_changes]
            A_1 = np.concatenate((A_1, np.zeros(self.I * (self.n - 1)).reshape(self.I * (self.n - 1), 1)), axis=1)
            A_2 = A_2[:, col_changes]
            A_2 = np.concatenate((A_2, np.zeros(self.I).reshape(self.I, 1)), axis=1)
            con = [A_1, b_1]
            all_iv_old = all_iv.copy()
            x, all_iv = self.compute_fp(simplex, con, theta_ind)
            A_2 = np.concatenate((A_2, np.zeros(self.I).reshape(self.I, 1)), axis=1)
            con = [A_1, b_1]'''
            """all_iv_old = all_iv.copy()
            # x, all_iv = self.compute_fp(simplex, con, theta_ind)
            # x = self.lt_algo(simplex, con, theta_ind)
            all_iv_old_keys = all_iv_old.keys()
            for (i, v_ind) in all_iv_old_keys:
                delta_p = self.get_outgoing_edges(self.V[v_ind])
                # delta_p_act = self.get_outgoing_active_edges(v)
                # delta_p_inact = [e for e in delta_p if e not in delta_p_act]
                if (i, v_ind) not in all_iv:
                    for e_ind in delta_p:
                        if self.fp[i][e_ind][-1][1] > self.eps:
                            self.fp[i][e_ind].append((theta, 0))
                            self.fp_ind[i][e_ind].append(theta_ind)
                else:
                    for e_ind in delta_p_inact:
                        if self.fp[i][e_ind][-1][1] > self.eps:
                            self.fp[i][e_ind].append(theta, 0)"""
            '''
            x_sum = np.zeros(self.m)
            if not isinstance(x, int):
                dim = len(x)  # - 1
                rows, col_pointer = simplex[1:]
                co = 0
                for d in range(dim):
                    try:
                        while col_pointer[co+1] <= d:
                            co += 1
                    except IndexError:
                        pass
                    v_ind = self.V.index(self.E[co][0])
                    x_total[rows[d], co] = x[d] * b[rows[d]][v_ind] * len(all_iv.keys())
                    # x_total[rows[d], co] = x[d] * b[rows[d]][v_ind] * com_no[v_ind]
                for e_ind in range(self.m):
                    x_sum[e_ind] = np.sum([x_total[i, e_ind] for i in range(self.I)])
            '''

            for ti in range(self.I):
                for v in top_ords[ti]:
                    active_paths = self.get_outgoing_active_edges(ti, v)
                    # betrachte aktive Kanten
                    for e_ind in active_paths:
                        if theta == 0:
                            # if len(self.u[ti][v_ind]) > 0 and self.u[ti][v_ind][0][0] == 0:
                            if x_total[ti, e_ind] > self.eps:
                                self.fp[ti][e_ind][0] = (0, x_total[ti, e_ind])
                                self.fp_ind[ti][e_ind].append(0)
                        # falls sich f^+ -Wert in dieser Phase ändert, aktualisiere 'self.fp'
                        elif abs(self.fp[ti][e_ind][-1][1] - x_total[ti, e_ind]) > self.bd_tol:
                            #if theta - self.fp[ti][e_ind][-1][0] < self.bd_tol:
                            #    self.fp[ti][e_ind].pop()
                            #    self.fp_ind[ti][e_ind].pop()
                            #else:
                            self.fp[ti][e_ind].append((theta, x_total[ti, e_ind]))
                            self.fp_ind[ti][e_ind].append(theta_ind)
                        if x_sum[e_ind] > self.eps:
                            flow_ratio = x_total[ti, e_ind] / x_sum[e_ind]
                        else:
                            flow_ratio = 0
                        if self.q_global[theta_ind][e_ind] > self.eps:
                            outflow_time = theta + self.q_global[theta_ind][e_ind]/self.nu[e_ind] + self.r[e_ind]
                            outflow = self.nu[e_ind] * flow_ratio
                        else:
                            outflow_time = theta + self.r[e_ind]
                            if x_sum[e_ind] < self.nu[e_ind] + self.eps:
                                outflow = x_total[ti, e_ind]
                            else:
                                outflow = self.nu[e_ind] * flow_ratio
                        fm_ind = self.last_fm_change(ti, e_ind, outflow_time)
                        # falls sich f_{ti}^- -Wert durch die Flussaufteilung in dieser Phase ändert, aktualisiere 'self.fm'
                        if abs(self.fm[ti][e_ind][fm_ind][1] - outflow) > self.bd_tol:
                            #if outflow_time - self.fm[ti][e_ind][fm_ind][0] < self.bd_tol:
                            #    # 'self.fm' ändert sich innerhalb eines Zeitintervalls der Länge 'self.eps' zweimal -> nur wegen Rundungsfehler -> wird stattdessen gar nicht
                            #    # geändert
                            #    self.fm[ti][e_ind].pop()
                            #else:
                            self.fm[ti][e_ind].append((outflow_time, outflow))

            firstcom = self.m * [True]
            for ti in range(self.I):
                for v in top_ords[ti]:
                    active_paths = self.get_outgoing_active_edges(ti, v)
                    for e_ind in active_paths:
                        last_fm_ind = self.last_fm_change(ti, e_ind, theta)
                        # bestimme, ob sich vor dem Ende der aktuellen Phase ein f^- -Wert ändert -> verkürze Phase
                        if len(self.fm[ti][e_ind]) > last_fm_ind + 1 and self.eps < self.fm[ti][e_ind][last_fm_ind + 1][0] - theta < next_phase - self.eps:
                            next_phase = self.fm[ti][e_ind][last_fm_ind + 1][0] - theta
                            # wird zurückgesetzt, da durch Verkürzung der Phase der f^- -Wert erst in einer späteren Phase
                            # neu gesetzt wird
                            fm_to_zero = []

                        change_of_q = self.change_of_cost(e_ind, theta, x_sum[e_ind]) * self.nu[e_ind]
                        # Prüfe, ob die nachstehenden Zeilen bereits für dieselbe Kante für ein anderes Gut ausgeführt wurden. Falls nein, werden sie ausgeführt, falls ja, werden
                        # sie übersprungen
                        if firstcom[e_ind]:
                            # Falls die Warteschlange von 'e' unter aktuellem Fluss abgebaut wird, bestimme Zeitpunkt, zu dem diese vollständig abgebaut ist (bei gleich bleibendem
                            # Fluss)
                            if change_of_q < -self.eps:
                                # 'phase_length': Dauer bis Warteschlangenlänge gleich 0
                                phase_length = - self.q_global[theta_ind][e_ind] / change_of_q
                                if phase_length < next_phase - self.eps:
                                    next_phase = phase_length
                                    # prüfe ob Zufluss 0 und somit 'fm' auf 0 gesetzt werden muss
                                    if change_of_q + self.nu[e_ind] < self.eps:
                                        fm_to_zero = [e_ind]
                                elif max([abs(phase_length - next_phase), change_of_q + self.nu[e_ind]]) < self.eps:
                                    fm_to_zero.append(e_ind)
                            firstcom[e_ind] = False

                # erneuter Durchlauf der vorherigen Schleife, diesmal werden die inaktiven Kanten betrachtet. Diese Aufteilung ist notwendig, damit die Änderungen der Knotenlabels
                # erst richtig gesetzt (siehe Schleife zuvor), und dann für weitere Rechnungen (siehe nachstehende Schleife) verwendet werden.
                for v in top_ords[ti]:
                    active_paths = self.get_outgoing_active_edges(ti, v)
                    delta_p = self.get_outgoing_edges(v)
                    inactive_paths = [e for e in delta_p if e not in active_paths]
                    for e_ind in inactive_paths:
                        # falls f^+ -Wert vorher > 0 war, so wird dieser hier auf 0 gesetzt, da Kante inaktiv
                        if abs(self.fp[ti][e_ind][-1][1]) > self.eps:
                            self.fp[ti][e_ind].append((theta, 0))
                            self.fp_ind[ti][e_ind].append(theta_ind)
                            if self.q_global[theta_ind][e_ind] > self.eps:
                                outflow_time = theta + self.q_global[theta_ind][e_ind]/self.nu[e_ind] + self.r[e_ind]
                            else:
                                outflow_time = theta + self.r[e_ind]
                            fm_ind = self.last_fm_change(ti, e_ind, outflow_time)
                            # setze 'self.fm'- Wert auf 0, falls dies noch nicht geschehen ist
                            if abs(self.fm[ti][e_ind][fm_ind][1]) > self.eps:
                                self.fm[ti][e_ind].append((outflow_time, 0))

                        # bestimme, ob sich vor dem Ende der aktuellen Phase ein f^- -Wert ändert -> verkürze Phase
                        last_fm_ind = self.last_fm_change(ti, e_ind, theta)
                        if len(self.fm[ti][e_ind]) > last_fm_ind + 1 and self.eps < self.fm[ti][e_ind][last_fm_ind + 1][0] - theta < next_phase - self.eps:
                            next_phase = self.fm[ti][e_ind][last_fm_ind + 1][0] - theta
                            # wird zurückgesetzt, da durch Verkürzung der Phase der f^- -Wert erst in einer späteren Phase
                            # neu gesetzt wird
                            fm_to_zero = []

                        change_of_q = self.change_of_cost(e_ind, theta, x_sum[e_ind]) * self.nu[e_ind]
                        # Falls Warteschlange existiert und abgebaut wird, bestimme Zeitpunkt, zu dem diese vollständig abgebaut ist
                        if change_of_q < -self.eps:
                            phase_length = -self.q_global[theta_ind][e_ind] / change_of_q
                            if phase_length < next_phase - self.eps:
                                next_phase = phase_length
                                # prüfe ob Zufluss 0 und somit 'fm' auf 0 gesetzt werden muss
                                if change_of_q + self.nu[e_ind] < self.eps:
                                    fm_to_zero = [e_ind]
                            elif max([abs(phase_length - next_phase), change_of_q + self.nu[e_ind]]) < self.eps:
                                fm_to_zero.append(e_ind)

                    # v_ind = self.V.index(v)
                    # if self.b[theta_ind][ti, v_ind] < self.eps:
                    #    # Ist kein Fluss von Gut 'ti' in 'v' vorhanden, so hat nachstehender Block keine Auswirkung auf aktuelle Flussaufteilung und wird daher übersprungen
                    #    continue
                    len_act = len(active_paths)
                    if len_act:
                        for i in range(len_act):
                            active_ind = active_paths[i]
                            # bestimme Kante, die während der gesamten Phase für Gut 'ti' aktiv bleibt
                            if x_total[ti, active_ind] > 0 or i == len_act - 1:
                                # Änderung der Kosten dieser Kante
                                active_change = self.change_of_cost(active_ind, theta, x_sum[active_ind])
                                break

                        for e_ind in inactive_paths:
                            change = self.change_of_cost(e_ind, theta, x_sum[e_ind])
                            # prüfe, wann inaktive Kanten unter momentanem Einfluss aktiv werden
                            tar_ind = self.V.index(self.E[e_ind][1])
                            act_ind = self.V.index(self.E[active_ind][1])
                            if self.labels[ti][theta_ind][tar_ind] + self.q_global[-1][e_ind]/self.nu[e_ind] + self.r[e_ind] + \
                                    (change + a[ti, tar_ind]) * next_phase < \
                                    self.labels[ti][theta_ind][act_ind] + self.q_global[-1][active_ind]/self.nu[active_ind] + \
                                    self.r[active_ind] + (active_change + a[ti, act_ind]) * next_phase and \
                                    abs(change + a[ti, tar_ind] - active_change - a[ti, act_ind]) \
                                    > self.eps:
                                time_ub = np.abs((self.labels[ti][theta_ind][tar_ind] + self.q_global[-1][e_ind]/self.nu[e_ind]
                                                  + self.r[e_ind] - self.labels[ti][theta_ind][act_ind]
                                                  - self.q_global[-1][active_ind]/self.nu[active_ind] - self.r[active_ind])
                                                 / (active_change + a[ti, act_ind] - change - a[ti, tar_ind]))
                                if time_ub < next_phase:
                                    next_phase = time_ub
                                    # wird zurückgesetzt, da durch Verkürzung der Phase f^- -Wert erst in einer späteren Phase neu gesetzt wird
                                    fm_to_zero = []

            if next_phase != T:
                # aktualisiere Warteschlangenlängen und Kosten
                new_q_global = []
                self.c.append(np.zeros(self.m))
                for e_ind in range(self.m):
                    next_q_len = self.q_global[theta_ind][e_ind] + self.change_of_cost(e_ind, theta, x_sum[e_ind]) * self.nu[e_ind] * next_phase
                    if next_q_len < self.eps:
                        next_q_len = 0
                    new_q_global.append(next_q_len)
                    self.c[-1][e_ind] = new_q_global[-1] / self.nu[e_ind] + self.r[e_ind]
                # speichere aktuelle Warteschlangenlängen
                self.q_global.append(new_q_global)

            """
            print("Kanten mit positivem Einfluss zum Zeitpunkt", theta, ":")
            for v_ind in range(self.n):
                out_neighbors = self.get_outgoing_edges(top_ord_act[v_ind])
                for e_ind in out_neighbors:
                    if x_total[e_ind] > 0:
                        print(self.E[e_ind], x_total[e_ind])
            """

            theta += next_phase
            if next_phase != T:
                # speichere Phase
                self.global_phase.append(theta)
                self.E_active = scipy.sparse.lil_matrix((self.I, self.m))
                for i in range(self.I):
                    # ti = "t{}".format(i+1)
                    # self.labels[i].append(self.dijkstra(self.graphReversed, ti, len(self.global_phase) - 1, visited=[], distances={}))
                    self.labels[i].append([])
                    for v_ind in range(self.n):
                        self.labels[i][-1].append(self.labels[i][-2][v_ind] + next_phase * a[i, v_ind])
                    for (v_ind, v) in enumerate(self.V):
                        outneighbors = self.G[v].keys()
                        for w in outneighbors:
                            w_ind = self.V.index(w)
                            edge = self.E.index((v, w))
                            if abs(self.labels[i][-1][v_ind] - self.labels[i][-1][w_ind] - self.c[-1][edge]) < self.eps + 10 * self.bd_tol:
                                self.E_active[i, edge] = 1

                fm_to_zero = list(set(fm_to_zero))
                for e_ind in fm_to_zero:
                    for i in range(self.I):
                        if self.fm[i][e_ind][-1][1] > self.eps:
                            self.fm[i][e_ind].append((theta + self.r[e_ind], 0))
                            stop_outflow.append(theta + self.r[e_ind])

        for i in range(self.I):
            # am Ende sind alle f^+ -Werte 0
            for e in range(self.m):
                if abs(self.fp[i][e][-1][1]) > self.eps:
                    self.fp[i][e].append((theta - next_phase, 0))
                    self.fp_ind[i][e].append(theta_ind)

        # erzeuge Ausgabe
        OutputTableMulti(self.G, self.V, self.E, self.I, self.nu, self.fp, self.fp_ind, self.fm, self.q_global, self.global_phase, self.c, self.labels, self.flow_vol)
        phasess = self.global_phase
        """
        # 2. A+0 in B2 von C'
        fp38 = self.fp[1][1638]
        fm38 = self.fm[1][1638]
        #q38 = self.q_global[:][1638]

        fp47 = self.fp[1][1647]
        fm47 = self.fm[1][1647]
        #q47 = self.q_global[:][1647]

        fp48 = self.fp[1][1648]
        fm48 = self.fm[1][1648]
        #q48 = self.q_global[:][1648]

        fp49 = self.fp[1][1649]
        fm49 = self.fm[1][1649]
        #q49 = self.q_global[:][1649]

        fp50 = self.fp[1][1650]
        fm50 = self.fm[1][1650]
        #qq50 = self.q_global[:][1650]

        fp51 = self.fp[1][1651]
        fm51 = self.fm[1][1651]
        #qq51 = self.q_global[:][1651]

        fp52 = self.fp[1][1652]
        fm52 = self.fm[1][1652]
        #qq52 = self.q_global[:][1652]

        fp53 = self.fp[1][1653]
        fm53 = self.fm[1][1653]
        #qq53 = self.q_global[:][1653]

        fp54 = self.fp[1][1654]
        fm54 = self.fm[1][1654]
        #qq54 = self.q_global[:][1654]

        fp55 = self.fp[1][1655]
        fm55 = self.fm[1][1655]
        #qq55 = self.q_global[:][1655]

        fp56 = self.fp[1][1656]
        fm56 = self.fm[1][1656]
        #qq56 = self.q_global[:][1656]

        fp57 = self.fp[1][1657]
        fm57 = self.fm[1][1657]
        #qq57 = self.q_global[:][1657]

        fp58 = self.fp[1][1658]
        fm58 = self.fm[1][1658]
        #qq58 = self.q_global[:][1658]
        """

        # 8. Gadget (A+4) aus B7+4 in C
        """fp96 = self.fp[0][1596]
        fm96 = self.fm[0][1596]

        fp97 = self.fp[0][1597]
        fm97 = self.fm[0][1597]

        fp01 = self.fp[0][1601]
        fm01 = self.fm[0][1601]

        fp02 = self.fp[0][1602]
        fm02 = self.fm[0][1602]

        fp03 = self.fp[0][1603]
        fm03 = self.fm[0][1603]

        fp98 = self.fp[0][1598]
        fm98 = self.fm[0][1598]

        fp05 = self.fp[0][1605]
        fm05 = self.fm[0][1605]

        fp06 = self.fp[0][1606]
        fm06 = self.fm[0][1606]

        fp87 = self.fp[0][1587]
        fm87 = self.fm[0][1587]

        fp99 = self.fp[0][1599]
        fm99 = self.fm[0][1599]

        fp00 = self.fp[0][1600]
        fm00 = self.fm[0][1600]

        fp04 = self.fp[0][1604]
        fm04 = self.fm[0][1604]

        fp07 = self.fp[0][1607]
        fm07 = self.fm[0][1607]"""

        end_time = time.time()
        timediff1 = end_time - self.time_vor_main
        timediff2 = end_time - self.start_time
        return 0
