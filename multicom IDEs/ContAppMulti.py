import numpy as np
from OutputTableMulti import OutputTableMulti
from scipy.sparse.linalg import spsolve
from scipy.optimize import linprog
import scipy
from operator import itemgetter
from random import random


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
        self.items = G.items()
        self.keys = G.keys()
        self.eps = 10**(-9)  # Für Rundungsfehler
        self.b = np.zeros((self.I, self.n))  # merke Flusswerte je Gut für jeden Knoten
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

        self.fp = [[] for _ in range(self.I)]  # f^+    [[],[]]
        self.fp_ind = [[] for _ in range(self.I)]  # Liste der Indizes der Phasen
        self.fm = [[] for _ in range(self.I)]  # f^-
        for i in range(self.I):
            for k in range(self.m):
                self.fp[i].append([(0, 0)])
                self.fp_ind[i].append([])
                self.fm[i].append([(0, 0)])

        self.c.append(np.zeros(self.m))
        for e in self.E:  # Initialisierung "self.c" (Kosten)
            self.c[0][self.E.index(e)] = self.r[self.E.index(e)]

        # Zeitpunkte, zu denen sich der Zufluss in mindestens einem Quellknoten ändert
        self.u_start = set()
        for com in self.u:
            for s_list in com:
                for t in s_list:
                    self.u_start.add(t[0])

        self.graphReversed = self.reverse_graph(G)

        self.E_active = self.I * [[np.ones(self.m)]]
        for i in range(self.I):
            self.labels.append([self.dijkstra(self.graphReversed, "t{}".format(i+1), 0,visited=[], distances={})])

        self.E_active = []
        for i in range(self.I):
            self.E_active.append([np.zeros(self.m)])
        for v in self.V:
            v_ind = self.V.index(v)
            outneighbors = self.G[v].keys()
            for w in outneighbors:
                w_ind = self.V.index(w)
                edge = self.E.index((v, w))
                for i in range(self.I):
                    if abs(self.labels[i][0][v_ind] - self.labels[i][0][w_ind] - self.c[0][edge]) < self.eps:
                        self.E_active[i][0][edge] = 1

        self.del_plus_label = self.I * [np.zeros(self.n)]
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
                self.topologicalSortUtil(self.V[k], visited, stack)
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

    def calc_b(self, i, v, phase):
        """
        Berechnet zum Zeitpunkt 'phase' im Knoten 'v' vorhandene Flussmenge b_{i,v}^- (phase) von Gut 'i'
        :param v: Knoten
        :param phase: Zeitpunkt
        :return: b_{i,v}^- (phase)
        """
        b = 0
        v_ind = self.V.index(v)
        if v_ind == self.V.index('t{}'.format(i+1)):
            return 0
        in_paths = self.get_ingoing_edges(v)
        for e_ind in in_paths:
            ind = self.last_fm_change(i, e_ind, phase)
            if len(self.fm[i][e_ind]) > ind + 1 and abs(self.fm[i][e_ind][ind + 1][0] - phase) < self.eps:
                ind += 1
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

    def get_outgoing_active_edges(self, i, v, theta_ind=-1):
        """
        bestimmt alle aus 'v' ausgehende, momentan für Gut 'i' aktive Kanten
        :param i: Index des Guts
        :param v: Knoten
        :param theta_ind: Index des Zeitpunkts. Standardmäßig '-1', also der letzte berechnete Zeitpunkt
        :return: Liste der Indizes der aktiven Kanten
        """
        delta_p = [self.E.index((v,u)) for u in self.G[v].keys()]
        return [e for e in delta_p if self.E_active[i][theta_ind][e]]

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

    def find_cover_simplex(self, theta_ind):
        simplex = [[], [], []]  # CCS - Format für Eckpunkte des Simplex
        # Dimension: Anzahl Einträge, die ungleich 0 sein können -> len(S[0])
        for e_ind in range(self.m):
            e_start_ind = self.V.index(self.E[e_ind][0])
            simplex[2].append(len(simplex[0]))
            for i in range(self.I):
                if self.b[i][e_start_ind] > self.eps and self.E_active[i][theta_ind][e_ind]:
                    simplex[0].append(self.b[i][e_start_ind])  # Wert
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
        return b

    '''def barycentric_trafo_sp(self, simplex, x):
        n = len(x)
        s1 =
        s2 = np.ones(n)
        b = spsolve(A,x)'''

    def find_ancestor_simplex(self, k, p):
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
        # so ok ???
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

        '''
        # if any(x < -self.eps):
        x[0] = np.ceil(q[0])
        lam[0] = round(x[0] - q[0], 9)
        for i in range(1, n-1):
            x[i] = np.ceil(q[i] - lam[i-1])
            lam[i] = round(x[i] - q[i] + lam[i-1], 9)
        x[n-1] = q[n-1] - lam[n-2]'''

        pi = range(1, n)
        sort = sorted(zip(lam, pi),  key=lambda y: (y[0], -y[1]), reverse=True)
        pi = [element for _, element in sort]
        return x, pi

    def in_k(self, x, con, all_iv, k):
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
        return False

    def g(self, e_ind, theta_ind, x):
        if self.q_global[theta_ind][e_ind] > self.eps:
            return x - self.nu[e_ind]
        else:
            return max([x - self.nu[e_ind], 0])

    def find_fv(self, node, theta_ind, all_iv, simplex):
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
                        xsum += self.b[rows[ro]][v_ind] * node[ro] / sum([node[d] for d in all_iv[(rows[ro], v_ind)]])
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
                if self.b[i][v_ind] < self.eps:
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
        return X[:, part], f_X[:, part]

    def compute_fp(self, simplex, con, theta_ind):

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
            print("k", k, "FPK", fpk)
            if len(fpk) > 2 and np.linalg.norm(fpk[-2] - fpk[-1]) < self.eps:
                print("FOLGE", fpk)
                return fpk[-1], all_iv
            # HIER WEITER: folgendes wurde immer für jedes k gemacht (ineffizient da gleich). jetzt lieber für späteres
            # k altes x0, pi (angepasst) verwenden um in der nähe des alten fixpunkts zu starten?
            if k == 3:
                start_node = cs_bz
            else:
                start_node = fpk[-1]
            #start_node = np.zeros(dim)
            #start_node[-1] = 1
            print("startnode", start_node)
            xc0 , pi_c = self.find_ancestor_simplex(k, start_node)
            Xc, f_Xc = self.fk_on_tk(xc0, pi_c, con, simplex, all_iv, theta_ind, cs_bz, k)

            # bestimme baryzentrische Koordinaten von cs_bz bzgl. Xc[:, 0], ..., Xc[:, dim]
            bz = self.bary_coords(Xc, start_node)
            '''
            while np.any(bz < self.eps):
                bz0 = np.where(bz < self.eps)[0][0]
                cs_bz[bz0] += 1/(k+3) 
                for i in range(dim):
                    if cs_bz[i] > 1/(k+3) + self.eps and i != bz0:
                        cs_bz[i] -= 1/(k+3)
                        break
                bz = self.bary_coords(Xc, cs_bz)'''
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
            '''
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
                        break'''
            # Perturbation von 'r'
            r += v_eps
            # HIER WEITER: perturbation von r so nicht korrekt!? da epsilon so klein wird r nicht ausreichend perturbiert!?
            #r = np.zeros(dim)
            #r[0] = -0.1
            #r[1] = -0.66
            #r[2] = 0.45
            #r[3] = 0.31
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
                # HIER WEITER: nochmal drüber schauen:
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
                # HIER WEITER: wieso geht z gegen 1??

                if zss < self.eps:
                    xk = np.zeros(dim)
                    for d in range(dim):
                        xk += yss[d] * X[:, d]
                    fpk.append(xk)
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
                    '''
                    # table 1 falsch!? eigener Versuch: Vertausche Spalten in L:
                    L_old = L.copy()
                    if t < dim - 2:
                        L = np.concatenate((L_old[:, :t], L_old[:, t+1].reshape((L_rows, 1)), L_old[:, t].reshape((L_rows, 1)), L_old[:, t+2:]), axis=1)
                    else:
                        # t = dim - 2
                        try:
                            L = np.concatenate((L_old[:, :t], L_old[:, t+1].reshape((L_rows, 1)), L_old[:, t].reshape((L_rows, 1))), axis=1)
                        except TypeError:
                            print("da")'''
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
            fv0 = self.fk_on_tk(v0, [], con, simplex, all_iv, theta_ind, cs_bz, k, part=[0])[1]
            lv0 = fv0 - v0 + np.ones(dim)
            base = np.eye(dim)
            base_v = np.zeros((dim, dim))
            base_y = np.ones(dim)
            wt_next = v0
            T = []
            gamma = []
            R = np.zeros(dim)
            L = [lv0]
            tauw = [v0]
            t = 0
            mu_no = dim
            while True:
                if mu_no == 0:
                    clss.append(tauw)
                    k += 1
                    break
                # Im Fall 'part=[0]' ist 'pi' irrelevant
                fvt_next = self.fk_on_tk(wt_next, [], con, simplex, all_iv, theta_ind, cs_bz, k, part=[0])[1]
                lvt_next = fvt_next - wt_next.reshape((dim, 1)) + np.ones(dim).reshape((dim, 1))
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
                    wt_next = base_v[:, out]
                    base[:, out] = lvt_next.reshape(dim)
                    base_v[:, out] = wt_next
                else:
                    mu_no -= 1
                    qi = np.zeros(dim)
                    i = np.where(base[:, out] > self.eps)[0][0]
                    qi[i] = -1
                    qi[i+1] = 1
                    tauw.append(tauw[-1] + qi)
                    base[:, out] = lvt_next.reshape(dim)
                    base_v[:, out] = wt_next
                    wt_next = tauw[-1]
                    T.append(i)
                    gamma.append(i)
                    # HIER WEITER: w^1 enthält Eintrag -1? R wird nie angepasst und somit wird einziger eintrag aus gamma gelöscht -> indexError
                    t += 1
                    continue
                # if lvt_next in L:
                    # s = L.index(lvt_next)
                while True:
                    if out == 0:
                        qgamma1 = np.zeros(dim)
                        qgamma1[gamma[0]] = -1
                        qgamma1[gamma[0]+1] = 1
                        tauw[0] += qgamma1
                        g0 = gamma[0]
                        gamma = np.delete(gamma, 0)
                        gamma = np.append(gamma, g0)
                        R[g0] += 1
                        qj = np.zeros(dim)
                        qj[gamma[0]] = -1
                        qj[gamma[0]+1] = 1
                        tauw.append(tauw[-1] + qj)
                        wt_next = tauw[-1]
                        del tauw[0]
                        del L[0]
                        L.append(lvt_next)
                    elif out == t:
                        qgammat = np.zeros(dim)
                        gt = gamma[-1]
                        gamma = np.delete(gamma, -1)
                        gamma = np.insert(gamma, 0, gt)
                        qgammat[gt] = -1
                        qgammat[gt + 1] = 1
                        tauw[0] -= qgammat
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
                                if base_y[j] - lam[j] * delt < -self.eps:
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
                                tauw.append(tauw[-1] + qi)
                                wt_next = tauw[-1]
                                T.append(i)
                                gamma.append(i)
                                t += 1
                                break
                        else:
                            R[gt] -= 1
                            del tauw[-1]
                            tauw.append(tauw[0])
                            del L[-1]
                            L.append(lvt_next)
                            break
                    else:
                        gs = gamma[out]
                        gamma[out] = gamma[out-1]
                        gamma[out-1] = gs
                        qj = np.zeros(dim)
                        qj[gamma[out]] = -1
                        qj[gamma[out]+1] = 1
                        tauw.append(tauw[out] + qj)
                        del tauw[out]
                        del L[out]
                        L.append(lvt_next)
                        break
                # else:
                #    tauw.append(wt_next)
                #    L.append(lvt_next)
                #    T.append(lvt_next)
                #    gamma.append(lvt_next)

    def main(self):
        """
        Hauptteil: Konstuiert schrittweise einen kontinuierlichen IDE-Fluss zu den Eingabedaten aus 'Math.cont_data'. Erzeugt Tabelle der
        Klasse 'Graphics.OutputTable' mit den entsprechenden Daten des Flusses.
        :return: 0
        """
        theta = 0
        # Obergrenze für theta
        T = 10000
        stop_outflow = []
        simplex = [[], [], []]
        all_iv = {}
        while theta < T:
            # in der Zukunft liegende Zeitpunkte aus der Liste 'self.u_start'
            start_points = [t for t in self.u_start if t > theta]
            # in der Zukunft liegende Zeitpunkte, zu denen der f^- -Wert von mindestens einer Kante auf 0 springt (wird
            # während der Laufzeit aktualisiert)
            stop_outflow = [t for t in stop_outflow if t > theta]
            theta_ind = self.global_phase.index(theta)
            # Liste aller Kanten, deren Warteschlange in der aktuellen Phase 0 wird, und deren f^- -Werte auf 0
            # gesetzt werden müssen
            fm_to_zero = []
            # Aufteilung des Flusses
            x_total = np.zeros((self.I, self.m))
            # 'next_phase' bestimmt am Ende der Schleife, wie lange aktuelle Phase andauert
            if len(start_points) > 0 or len(stop_outflow) > 0:
                next_phase = np.min(start_points + stop_outflow) - theta
            else:
                next_phase = T
            self.flow_vol = np.zeros(self.n)
            # Anz. Güter je Knoten, die aktuell eine Entscheidung über Flussaufteilung treffen müssen
            com_no = np.zeros(self.n)
            for v_ind in range(self.n):
                v = self.V[v_ind]
                for i in range(self.I):
                    # speichere b - Wert
                    self.b[i][v_ind] = self.calc_b(i, v, theta)
                    if self.b[i][v_ind] > self.eps:
                        com_no[v_ind] += 1
                        self.flow_vol[v_ind] += self.b[i][v_ind]
            simplex = self.find_cover_simplex(theta_ind)
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
            b_1 = np.ones(A_1.shape[0])

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
                    # b_1[i * (self.n-1) + vi_ind] = self.b[i][v_ind]
                    # nicht b_{i,v}, sondern seine bary. Koordinaten, da wir in Standardsimplex arbeiten
                    if self.b[i][v_ind] > self.eps:
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
                    if self.b[i][v_ind] < self.eps:
                        delta_p = self.get_outgoing_edges(self.V[v_ind])
                        for e_ind in delta_p:
                            zero_cols.append(i * self.m + e_ind)
            col_changes_cpy = col_changes.copy()
            col_changes = [c for c in col_changes_cpy if c not in zero_cols]
            A_1 = A_1[:, col_changes]
            A_1 = np.concatenate((A_1, np.zeros(self.I * (self.n - 1)).reshape(self.I * (self.n - 1), 1)), axis=1)
            A_2 = A_2[:, col_changes]
            A_2 = np.concatenate((A_2, np.zeros(self.I).reshape(self.I, 1)), axis=1)'''
            con = [A_1, b_1]
            all_iv_old = all_iv.copy()
            # x, all_iv = self.compute_fp(simplex, con, theta_ind)
            x = self.lt_algo(simplex, con, theta_ind)
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
                '''
                else:
                    for e_ind in delta_p_inact:
                        if self.fp[i][e_ind][-1][1] > self.eps:
                            self.fp[i][e_ind].append(theta, 0)'''

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
                    x_total[rows[d], co] = x[d] * self.b[rows[d]][v_ind] * len(all_iv.keys())
                    # x_total[rows[d], co] = x[d] * self.b[rows[d]][v_ind] * com_no[v_ind]
                for e_ind in range(self.m):
                    x_sum[e_ind] = np.sum([x_total[i, e_ind] for i in range(self.I)])

            top_ord = []
            all_edge_outflows = [[] for _ in range(self.m)]
            for ti in range(self.I):
                top_ord.append(self.topologicalSort(ti))
                for v in top_ord[ti]:
                    v_ind = self.V.index(v)
                    active_paths = self.get_outgoing_active_edges(ti, v)
                    # if self.flow_vol[v_ind] > self.eps:
                    #     flow_ratio = self.b[ti][v_ind]/self.flow_vol[v_ind]
                    # else:
                    #     flow_ratio = 0
                    # betrachte aktive Kanten
                    for e_ind in active_paths:
                        if theta == 0:
                            if len(self.u[ti][v_ind]) > 0 and self.u[ti][v_ind][0][0] == 0:
                                self.fp[ti][e_ind][0] = (0, x_total[ti, e_ind])
                                self.fp_ind[ti][e_ind].append(0)
                                if x_total[ti, e_ind] > self.eps:
                                    all_edge_outflows[e_ind].append(ti)
                        # falls sich f^+ -Wert in dieser Phase ändert, aktualisiere 'self.fp'
                        elif abs(self.fp[ti][e_ind][-1][1] - x_total[ti, e_ind]) > self.eps:
                            self.fp[ti][e_ind].append((theta, x_total[ti, e_ind]))
                            self.fp_ind[ti][e_ind].append(theta_ind)
                            if x_total[ti, e_ind] > self.eps:
                                all_edge_outflows[e_ind].append(ti)
                        if x_sum[e_ind] > self.eps:
                            flow_ratio = x_total[ti, e_ind] / x_sum[e_ind]
                        else:
                            flow_ratio = 0
                        if self.q_global[theta_ind][e_ind] > self.eps:
                            outflow_time = theta + float(self.q_global[theta_ind][e_ind])/self.nu[e_ind] + self.r[e_ind]
                            outflow = self.nu[e_ind] * flow_ratio
                        else:
                            outflow_time = theta + self.r[e_ind]
                            if x_sum[e_ind] < self.nu[e_ind] + self.eps:
                                outflow = x_total[ti, e_ind]
                            else:
                                outflow = self.nu[e_ind] * flow_ratio
                        fm_ind = self.last_fm_change(ti, e_ind, outflow_time)
                        # falls sich f_{ti}^- -Wert durch die Flussaufteilung in dieser Phase ändert, aktualisiere 'self.fm'
                        if abs(self.fm[ti][e_ind][fm_ind][1] - outflow) > self.eps:
                            self.fm[ti][e_ind].append((outflow_time, outflow))

            # Alles falsch?
            '''
            for e_ind in range(self.m):
                if len(all_edge_outflows[e_ind]) > 1:
                    total_outflow = 0
                    for i in all_edge_outflows[e_ind]:
                        total_outflow += self.fm[i][e_ind][-1][1]
                    for i in all_edge_outflows[e_ind]:
                        outflow_time = self.fm[i][e_ind][-1][0]
                        outflow = self.fm[i][e_ind][-1][1] / total_outflow
                        self.fm[i][e_ind][-1] = (outflow_time, outflow)'''

            firstcom = self.m * [True]
            for ti in range(self.I):
                for v in top_ord[ti]:
                    v_ind = self.V.index(v)
                    active_paths = self.get_outgoing_active_edges(ti, v)
                    firstedge = True
                    for e_ind in active_paths:
                        # bestimme, ob sich vor dem Ende der aktuellen Phase ein f^- -Wert ändert -> verkürze Phase
                        last_fm_ind = self.last_fm_change(ti, e_ind, theta)
                        if len(self.fm[ti][e_ind]) > last_fm_ind + 1 and \
                                self.eps < self.fm[ti][e_ind][last_fm_ind + 1][0] - theta < next_phase - self.eps:
                            next_phase = self.fm[ti][e_ind][last_fm_ind + 1][0] - theta
                            # wird zurückgesetzt, da durch Verkürzung der Phase der f^- -Wert erst in einer späteren Phase
                            # neu gesetzt wird
                            fm_to_zero = []

                        change_of_c = self.change_of_cost(e_ind, theta, x_sum[e_ind])
                        change_of_q = change_of_c * self.nu[e_ind]
                        new_del_plus = change_of_c + self.del_plus_label[ti][self.V.index(self.E[e_ind][1])]
                        # prüfe, ob Änderung des Labels von Knoten 'v' angepasst werden muss
                        if firstedge or new_del_plus < self.del_plus_label[ti][v_ind]:
                            self.del_plus_label[ti][v_ind] = new_del_plus
                            firstedge = False
                        # Prüfe, ob die nachstehenden Zeilen bereits für die selbe Kante für ein anderes Gut ausgeführt wurden. Falls
                        # nein, werden sie ausgeführt, falls ja, werden sie übersprungen
                        if firstcom[e_ind]:
                            # Falls die Warteschlange von 'e' unter aktuellem Fluss abgebaut wird, bestimme Zeitpunkt, zu dem diese
                            # vollständig abgebaut ist (bei gleich bleibendem Fluss)
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

                # erneuter Durchlauf der vorherigen Schleife, diesmal werden die inaktiven Kanten betrachtet. Diese Aufteilung ist
                # notwendig, damit die Änderungen der Knotenlabels erst richtig gesetzt (siehe Schleife zuvor), und dann für weitere
                # Rechnungen (siehe nachstehende Schleife) verwendet werden.
                for v in top_ord[ti]:
                    active_paths = self.get_outgoing_active_edges(ti, v)
                    delta_p = self.get_outgoing_edges(v)
                    inactive_paths = [e for e in delta_p if e not in active_paths]
                    for e_ind in inactive_paths:
                        # falls f^+ -Wert vorher > 0 war, so wird dieser hier auf 0 gesetzt, da Kante inaktiv
                        if abs(self.fp[ti][e_ind][-1][1]) > self.eps:
                            self.fp[ti][e_ind].append((theta, 0))
                            self.fp_ind[ti][e_ind].append(theta_ind)
                            if self.q_global[theta_ind][e_ind] > self.eps:
                                outflow_time = theta + self.q_global[theta_ind][e_ind] + self.r[e_ind]
                            else:
                                outflow_time = theta + self.r[e_ind]
                            fm_ind = self.last_fm_change(ti, e_ind, outflow_time)
                            # setze 'self.fm'- Wert auf 0, falls dies noch nicht geschehen ist
                            if abs(self.fm[ti][e_ind][fm_ind][1]) > self.eps:
                                self.fm[ti][e_ind].append((outflow_time, 0))

                        # bestimme, ob sich vor dem Ende der aktuellen Phase ein f^- -Wert ändert -> verkürze Phase
                        last_fm_ind = self.last_fm_change(ti, e_ind, theta)
                        if len(self.fm[ti][e_ind]) > last_fm_ind + 1 and \
                                self.eps < self.fm[ti][e_ind][last_fm_ind + 1][0] - theta < next_phase - self.eps:
                            next_phase = self.fm[ti][e_ind][last_fm_ind + 1][0] - theta
                            # wird zurückgesetzt, da durch Verkürzung der Phase der f^- -Wert erst in einer späteren Phase
                            # neu gesetzt wird
                            fm_to_zero = []

                        change_of_q = self.change_of_cost(e_ind, theta, x_sum[e_ind]) * self.nu[e_ind]
                        # Falls Warteschlange existiert und abgebaut wird, bestimme Zeitpunkt, zu dem diese vollständig abgebaut ist
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
                                    (change + self.del_plus_label[ti][tar_ind]) * next_phase < \
                                    self.labels[ti][theta_ind][act_ind] + self.q_global[-1][active_ind]/self.nu[active_ind] + \
                                    self.r[active_ind] + (active_change + self.del_plus_label[ti][act_ind]) * next_phase and \
                                    abs(change + self.del_plus_label[ti][tar_ind] - active_change - self.del_plus_label[ti][act_ind]) \
                                    > self.eps:
                                time_ub = np.abs((self.labels[ti][theta_ind][tar_ind] + self.q_global[-1][e_ind]/self.nu[e_ind]
                                                  + self.r[e_ind] - self.labels[ti][theta_ind][act_ind]
                                                  - self.q_global[-1][active_ind]/self.nu[active_ind] - self.r[active_ind])
                                                 / (active_change + self.del_plus_label[ti][act_ind] - change -
                                                    self.del_plus_label[ti][tar_ind]))
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
                                 self.change_of_cost(e_ind, theta, x_sum[e_ind]) * self.nu[e_ind] * next_phase
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
                for i in range(self.I):
                    ti = "t{}".format(i+1)
                    self.labels[i].append(self.dijkstra(self.graphReversed, ti, len(self.global_phase) - 1, visited=[], distances={}))
                    self.E_active[i].append(np.zeros(self.m))
                    for v in self.V:
                        v_ind = self.V.index(v)
                        outneighbors = self.G[v].keys()
                        for w in outneighbors:
                            w_ind = self.V.index(w)
                            edge = self.E.index((v, w))
                            if abs(self.labels[i][-1][v_ind] - self.labels[i][-1][w_ind] - self.c[-1][edge]) < self.eps:
                                self.E_active[i][-1][edge] = 1

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
        OutputTableMulti(self.V, self.E, self.I, self.nu, self.fp, self.fp_ind, self.fm, self.q_global, self.global_phase, self.c,
                         self.labels, self.flow_vol)
        return 0
