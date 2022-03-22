import numpy as np
from OutputTableMulti import OutputTableMulti
from scipy.sparse.linalg import spsolve
from scipy.optimize import linprog


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
        self.labels = self.I*[[]]  # Knotenlabels
        self.V = list(G.keys())  # Liste der Knoten
        self.n = len(self.V)  # Anzahl Knoten
        self.items = G.items()
        self.keys = G.keys()
        self.eps = 10**(-8)  # Für Rundungsfehler
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

        self.fp = self.I * [[]]  # f^+    [[],[]]
        self.fp_ind = self.I * [[]]  # Liste der Indizes der Phasen
        self.fm = self.I * [[]]  # f^-
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
            self.labels[i].append(self.dijkstra(self.graphReversed, "t{}".format(i+1), 0,visited=[], distances={}))

        self.E_active = self.I * [[np.zeros(self.m)]]
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
        ones = np.ones(n)
        ones = ones.reshape((1, n))
        v = np.concatenate((v, ones))
        x = np.append(x, 1)
        q, r = np.linalg.qr(v)
        z = q.T.dot(x)
        b = np.linalg.solve(r, z)
        # b = np.linalg.solve(v, x)  # anschauen: multigrid verfahren
        return b

    '''def barycentric_trafo_sp(self, simplex, x):
        n = len(x)
        s1 =
        s2 = np.ones(n)
        b = spsolve(A,x)'''

    def find_ancestor_simplex(self, k, p):
        # kann Fall xi < 0 auftreten? - Wenn ja -> noch behandeln
        n = len(p) - 1
        q = k * p
        x = np.zeros(n+1)
        lam = np.zeros(n)
        x[0] = np.floor(q[0] + 1)
        lam[0] = x[0] - q[0]
        for i in range(1, n):
            x[i] = np.floor(q[i] - lam[i-1] +1)
            lam[i] = x[i] - q[i] + lam[i-1]
        x[n] = q[n] - lam[n-1]
        pi = range(1, n+1)
        zipped = sorted(zip(lam, pi), reverse=True)
        pi = [element for _, element in zipped]
        return x, pi

    def in_k(self, x, con):
        A, b = con
        if np.linalg.norm(A.dot(x) - b) < self.eps:  # ?
            return True
        return False

    def g(self, e_ind, theta_ind, x):
        if self.q_global[theta_ind][e_ind] > self.eps:
            return x - self.nu[e_ind]
        else:
            return max([x - self.nu[e_ind], 0])

    def fk_on_tk(self,x0, pi, con, simplex, theta_ind, c, k):

        [vals, rows, col_pointer] = simplex
        dim = len(vals)
        U1 = np.eye(dim+1, k=-1) - np.eye(dim+1)
        U2 = np.zeros(dim+1)
        U2[-1] = 1
        U = np.concatenate((U1, U2.reshape((1,dim+1))))
        X = np.zeros((dim+1, dim+1))
        X[:][0] = 1/k * x0
        f_X = -np.ones((dim+1, dim+1))
        for j in range(1, dim+1):
            X[:][j] = X[:][j-1] + U[:][pi[j-1]]
            g_e = np.zeros(self.m)
            for e_ind in range(self.m):
                if col_pointer[e_ind] - col_pointer[e_ind + 1]:
                    g_e[e_ind] = self.g(e_ind, theta_ind, sum(X[col_pointer[e_ind]:col_pointer[e_ind+1], j]))
                else:
                    g_e[e_ind] = self.g(e_ind, theta_ind, 0)
            if self.in_k(X[:][j], con):
                # A, b = con
                # gamma_A = A.copy()
                # gamma_b = b.copy()
                # gamma auswerten
                gamma_zero = np.zeros(dim+1)
                for i in range(self.I):
                    act_edges = [e_ind for e_ind in range(self.m) if self.E_active[i][theta_ind][e_ind]]
                    top_ord = self.topologicalSort(i)[1:]
                    a = np.zeros((self.I, self.n))
                    for v in top_ord:
                        if i == 2:
                            print("flasc")
                        delta_p_act = self.get_outgoing_active_edges(i, v)
                        a[i][v] = np.Inf
                        for e_ind in delta_p_act:
                            w = self.V.index(self.E[e_ind][1])
                            if a[i][v] > g_e[e_ind] + a[i][w]:
                                a[i][v] = g_e[e_ind] + a[i][w]

                    for e_ind in act_edges:
                        v, w = self.E[e_ind][:]
                        v_ind = self.V.index(v)
                        w_ind = self.V.index(w)
                        if a[i][v_ind] < g_e[e_ind]/self.nu[e_ind] + a[i][w_ind]:
                            col_start = col_pointer[e_ind]
                            col_end = col_pointer[e_ind+1]
                            col_skip = len([r for r in rows[col_start:col_end] if r < i])
                            gamma_zero[col_start + col_skip] = 1
                            # e_i = np.zeros(dim+1)
                            # e_i[col_start + col_skip] = 1
                            # np.concatenate((gamma_A, e_i))
                            # np.concatenate((gamma_b, 0))
                            f_X[col_start + col_skip][j] = 0
                e_ind = 0
                for d in range(dim + 1):
                    if f_X[d][j] > -self.eps:
                        continue
                    try:
                        while col_pointer[e_ind+1] <= d:
                            e_ind += 1
                    except IndexError:
                        pass
                    i = rows[d]
                    v_ind = self.V.index(self.E[e_ind][0])
                    f_X[d][j] = self.b[i][v_ind]
                    if i == 2:
                        print("flasc")
                    delta_p = self.get_outgoing_active_edges(i, self.V[v_ind]) - [e_ind]
                    for e in delta_p:
                        elt_ind = rows[col_pointer[e]:].index(e)
                        f_X[elt_ind][j] = 0
            else:
                f_X[:][j] = c
        return X, f_X

    def compute_fp(self, simplex, con, theta_ind):

        [vals, rows, col_pointer] = simplex
        dim = len(vals)
        c = np.zeros(dim+1)
        for d in range(dim):
            v_ind = self.V.index(self.E[col_pointer[d]][0])
            if rows[d] == 2:
                print("flasc")
            len_dp_act = len(self.get_outgoing_active_edges(rows[d], self.V[v_ind], theta_ind=theta_ind))
            c[d] = self.b[rows[d]][v_ind] / len_dp_act

        # transformiere cs auf Standardsimplex durch Normieren
        c /= np.linalg.norm(c)
        k = 2
        fpk = []
        while True:
            if len(fpk) > 1 and np.linalg.norm(fpk[-2] - fpk[-1]) < self.eps:
                return fpk[-1]
            xc0 , pi_c = self.find_ancestor_simplex(k, c)
            Xc, f_Xc = self.fk_on_tk(xc0, pi_c, con, simplex, theta_ind, c, k)

            # bestimme baryzentrische Koordinaten von c bzgl. X[:][0], ..., X[:][dim]
            bz = self.bary_coords(Xc, c)
            # berechne f_k(c) mit den bary. Koordinaten von c bzgl. X und deren Funktionswerten f_k(X)
            fk_c = np.zeros(dim+1)
            for i in range(dim+1):
                # HIER WEITER: Welche der Dimensionen ist falsch. Habe heute dim von c auf dim gesetzt (statt dim+1)
                fk_c += bz[i] * f_Xc[:, i]

            r = fk_c - c  # cs = c + kappa * r, mit kappa max., s.d. cs >= 0.
            epsilon = 0.5
            for i in range(dim+1):
                if abs(r[i]) < epsilon:
                    epsilon = abs(r[i])
            epsilon *= 0.5
            v_eps = np.zeros(dim+1)
            for i in range(0, dim+1, 2):
                v_eps[i] = epsilon/(i+1)
            for i in range(1, dim+1, 2):
                v_eps[i] = -epsilon/i
            # falls 'dim' gerade, so wird zweite 'for'-Schleife einmal weniger durchlaufen als die erste
            if not dim % 2:
                v_eps[dim-1] -= epsilon/(dim+1)
            # Perturbation von 'r'
            r += v_eps
            kappa = np.Inf
            for i in range(dim+1):
                if c[i] + kappa * r[i] < -self.eps:
                    kappa = -c[i]/r[i]
            cs = c + kappa * r
            x0, pi = self.find_ancestor_simplex(k, cs)
            X, f_X = self.fk_on_tk(x0, pi, con, simplex, theta_ind, c, k)

            L = np.zeros((dim+2, dim+1))
            for d in range(dim+1):
                L[:dim+1, d] = f_X[:, d] - X[:, d]
                L[dim+1, d] = 1

            tf = np.concatenate((np.zeros(dim+2), np.ones(dim+2)))
            r = np.append(r, 1).reshape((dim+2, 1))
            A_total = np.concatenate((L, r, np.eye(dim+2)), axis=1)
            b_total = np.zeros(dim+2)
            b_total[-1] = 1
            basic_sol = linprog(tf, A_eq=A_total, b_eq=b_total)  # [:2*dim+5]
            y = basic_sol[:dim+1]
            z = basic_sol[dim+1]
            fun = basic_sol[-1]

            if fun > self.eps:
                raise TypeError('Basic feasible solution can\'t be found')
            if z < self.eps:
                xk = np.zeros(dim+1)
                for d in range(dim+1):
                    xk += y[d] * X[:][d]
                fpk.append(xk)
                k += 1
                continue
                # return x0, pi, y


            s = y.index(0)
            while True:
                if s == dim:
                    L_s = L[:][:s]
                else:
                    L_s = np.concatenate((L[:][:s], L[:][s+1:]))
                L_s = np.concatenate((L_s, r_total))
                lam = np.linalg.solve(L_s, L[:][s])  # Darstellung von 'L[:][s]' als Linearkomb. der Spalten von 'L_s'

                y_range = range(dim+1) - [s]
                delt = np.Inf
                t = np.Inf
                for col in y_range:
                    if y[col] - lam[col] * delt < self.eps:
                        delt = y[col]/lam[col]
                        t = col
                if z - lam[dim+1] * delt < self.eps:
                    delt = z/lam[dim+1]

                delta_y = np.concatenate((-lam[:s] * delt, delt, -lam[s:-1] * delt))
                delta_z = -lam[-1] * delt
                yss = y + delta_y
                zss = z + delta_z

                if zss < self.eps:
                    xk = np.zeros(dim+1)
                    for d in range(dim+1):
                        xk += y[d] * X[:][d]
                    fpk.append(xk)
                    break
                    # return x0, pi, yss

                if t == 0:
                    x0[pi[0]-1] -= 1.0/k
                    x0[pi[0]] += 1.0/k
                    j1 = pi.pop(0)
                    pi.append(j1)
                    yss.pop(0)
                    yss.append(0)
                    L = np.concatenate((L[:][1:], L[:][-1]))
                    s = dim
                elif t == dim:
                    x0[pi[-1]-1] += 1.0/k
                    x0[pi[-1]] -= 1.0/k
                    jn = pi.pop(-1)
                    pi.insert(0, jn)
                    yss.pop(-1)
                    yss.insert(0, 0)
                    L = np.concatenate((L[:][0], L[:][:-1]))
                    s = 0
                else:
                    pi_old = pi
                    pi[t-1] = pi[t]
                    pi[t] = pi_old[t-1]
                    s = t
                y = yss
                z = zss

            k += 1

        # wollen eine erste bfs bestimmen. phase 1 simplex? oder kann nicht-degeneriertheit verwendet werden? z.b. lin.
        # abh. Spalte finden -> bringt nicht so viel, weil dann rest des systems immer noch gelöst werden muss? also
        # einfach doch sofort lösen? -> phase 1 simplex. Aber wie funktioniert dann pivot schritt? pivotvariable leicht
        # zu bestimmen, da einzige 0-variable. dann linearkombination von dieser spalte mithilfe der anderen spalten
        # finden. Diese liefert infos über verhalten der anderen variablen bei erhöhung der pivotvariable. Dann einfach die
        # finden die zuerst 0 wird -> neue pivotvariable. so lange machen bis z=0 dann fertig! was passier dann in table
        # 1? keine änderung außer die eine vertauschung in pi? was soll "deleting X^t" bedeuten?



    def main(self):
        """
        Hauptteil: Konstuiert schrittweise einen kontinuierlichen IDE-Fluss zu den Eingabedaten aus 'Math.cont_data'. Erzeugt Tabelle der
        Klasse 'Graphics.OutputTable' mit den entsprechenden Daten des Flusses.
        :return: 0
        """
        theta = 0
        # Obergrenze für theta
        T = 10000
        # Aufteilung des Flusses
        x_total = self.I * [np.zeros(self.m)]
        stop_outflow = []
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
            # 'next_phase' bestimmt am Ende der Schleife, wie lange aktuelle Phase andauert
            if len(start_points) > 0 or len(stop_outflow) > 0:
                next_phase = np.min(start_points + stop_outflow) - theta
            else:
                next_phase = T
            self.flow_vol = np.zeros(self.n)
            for v_ind in range(self.n):
                v = self.V[v_ind]
                for i in range(self.I):
                    # speichere b - Wert
                    self.b[i][v_ind] = self.calc_b(i, v, theta)
                    self.flow_vol[v_ind] += self.b[i][v_ind]
            for ti in range(self.I):
                top_ord = self.topologicalSort(ti)
                for v in top_ord:
                    if ti == 2:
                        print("flasc")
                    simplex = self.find_cover_simplex(theta_ind)
                    A_1 = np.zeros((self.I * (self.n-1), self.I * self.m))
                    A_2 = np.zeros((self.I, self.I * self.m))

                    b_1 = np.zeros(self.I * (self.n-1))
                    b_2 = np.zeros(self.I)

                    for i in range(self.I):
                        V_i = list(set(self.V) - set(['t{}'.format(i+1)]))
                        for v in V_i:
                            vi_ind = V_i.index(v)
                            b_1[i * (self.n-1) + vi_ind] = self.b[i][self.V.index(v)]
                            delta_v = self.get_outgoing_edges(v)  # Liste der Kantenindizes
                            for e_ind in delta_v:
                                A_1[i * (self.n-1) + vi_ind][i * self.m + e_ind] = 1

                        for v in self.V:
                            if i == 2:
                                print("flasc")
                            delta_v_inact = list(set(self.get_outgoing_edges(v)) - set(self.get_outgoing_active_edges(i, v)))
                            for e_ind in delta_v_inact:
                                A_2[i][i * self.m + e_ind] = 1
                    col_changes = [0]
                    len_cp = len(simplex[2])
                    for cp in range(1, len_cp):
                        if simplex[2][cp] != simplex[2][cp-1]:
                            col_changes.append(cp)
                    # col_changes.pop(-1)
                    A_1 = A_1[:, col_changes]
                    A_2 = A_2[:, col_changes]
                    A = np.concatenate((A_1, A_2))
                    b = np.concatenate((b_1, b_2))
                    con = [A, b]
                    x = self.compute_fp(simplex, con, theta_ind)
                    for e_ind in active_paths:
                        x_total[ti][e_ind] += x[active_paths.index(e_ind)]
            x_sum = np.zeros(self.m)
            for e_ind in range(self.m):
                x_sum[e_ind] += np.sum([x_total[i][e_ind] for i in range(self.I)])

            firstcom = self.m * [True]
            for ti in range(self.I):
                top_ord = self.topologicalSort(ti)
                for v in top_ord:
                    v_ind = self.V.index(v)
                    if ti == 2:
                        print("flasc")
                    active_paths = self.get_outgoing_active_edges(ti, v)
                    bi = self.calc_b(ti, v, theta)
                    if self.flow_vol[v_ind] > self.eps:
                        flow_ratio = float(bi)/self.flow_vol[v_ind]
                    else:
                        flow_ratio = 0
                    firstedge = True
                    # betrachte aktive Kanten
                    for e_ind in active_paths:
                        e = self.E[e_ind]

                        if theta == 0:
                            start_ind = self.V.index(e[0])
                            if len(self.u[ti][start_ind]) > 0 and self.u[ti][start_ind][0][0] == 0:
                                self.fp[ti][e_ind][0] = (0, x_total[ti][e_ind])
                                self.fp_ind[ti][e_ind].append(0)
                        # falls sich f^+ -Wert in dieser Phase ändert, aktualisiere 'self.fp'
                        elif abs(self.fp[ti][e_ind][-1][1] - x_total[ti][e_ind]) > self.eps:
                            self.fp[ti][e_ind].append((theta, x_total[ti][e_ind]))
                            self.fp_ind[ti][e_ind].append(theta_ind)
                        if self.q_global[theta_ind][e_ind] > self.eps:
                            outflow_time = theta + self.r[e_ind]
                            outflow = self.nu[e_ind] * flow_ratio
                        else:
                            outflow_time = theta + float(self.q_global[theta_ind][e_ind])/self.nu[e_ind] + self.r[e_ind]
                            outflow = np.min([x_total[ti][e_ind], self.nu[e_ind]]) * flow_ratio
                        fm_ind = self.last_fm_change(ti, e_ind, outflow_time)
                        # falls sich f_{ti}^- -Wert durch die Flussaufteilung in dieser Phase ändert, aktualisiere 'self.fm'
                        if abs(self.fm[ti][e_ind][fm_ind][1] - outflow) > self.eps:
                            self.fm[ti][e_ind].append((outflow_time, outflow))

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
                        new_del_plus = change_of_c + self.del_plus_label[ti][self.V.index(e[1])]
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
                for v in top_ord:
                    if ti == 2:
                        print("flasc")
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
                            if x_total[ti][active_ind] > 0 or i == len_act - 1:
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

                for e_ind in fm_to_zero:
                    for i in range(self.I):
                        self.fm[i][e_ind].append((theta + self.r[e_ind], 0))
                        stop_outflow[i].append(theta + self.r[e_ind])

            print(theta)

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
