import tkinter as tk
import numpy as np
import functools


class OutputTable(object):

    def __init__(self, V, E, nu, fp, fp_ind, fm, q, phases, c, labels, vol):

        self.V = V  # Knotenmenge
        self.E = E  # Kantenmenge
        self.n = len(V)  # Anzahln Knoten
        self.m = len(E)  # Anzahl Kanten
        self.q = q  # Warteschlangenlängen aller Kanten für jede Phase
        self.phases = phases  # Zeitpunkte der Phasenänderungen
        self.nu = nu  # Kantenkapazitäten
        self.fp = fp  # f^+ -Werte
        self.fp_ind = fp_ind  # Liste der Indizes der Phasen
        self.fm = fm  # f^- -Werte
        self.c = c  # Kantenkosten
        self.labels = labels  # Knotenlabels
        self.flow_vol = vol  # Flussvolumen zu jedem Zeitpunkt in den einzelnen Knoten
        self.table = tk.Tk()  # Ausgabetabelle
        self.menu = tk.Menu(self.table)  # Menüleiste
        self.table.config(menu=self.menu)
        self.table.title("Zeit: 0")
        self.filemenuR = tk.Menu(self.menu)  # Menüpunkt Zeilen
        self.filemenuC = tk.Menu(self.menu)  # Menüpunkt Spalten
        self.menuStartingNodes = tk.Menu(self.filemenuR)  # Untermenü von 'self.filemenuR
        self.menu.add_cascade(label="Zeilen", menu=self.filemenuR)
        self.menu.add_cascade(label="Spalten", menu=self.filemenuC)

        # Variablen der Checkboxen zu den Zeilengruppen
        self.CheckVarAll = tk.BooleanVar()  # wenn 'True', zeige alle Zeilen
        self.CheckVarPosFlow = tk.BooleanVar()  # wenn 'True', zeige alle Zeilen mit momentan positivem Flussvolumen im
                                                # Startknoten der entsprechenden Kante
        self.CheckVarPosQ = tk.BooleanVar()  # wenn 'True', zeige alle Zeilen mit momentan positiver Warteschlangen-
                                             # länge
        self.CheckVarAll.set(True)
        self.CheckVarPosFlow.set(True)
        self.CheckVarPosQ.set(True)

        self.filemenuR.add_checkbutton(
            label="alle", variable=self.CheckVarAll, command=self.check_all_rows)
        self.filemenuR.add_checkbutton(label="Kanten mit positivem Flussvolumen im Startknoten",
                                       variable=self.CheckVarPosFlow, command=self.check_pos_flow)
        self.filemenuR.add_checkbutton(label="Kanten mit positiver Warteschlangenlänge",variable=self.CheckVarPosQ,
                                       command=self.check_pos_q)

        self.filemenuR.add_cascade(label="Kanten mit Startknoten ...", menu=self.menuStartingNodes)
        self.CheckVarsV = []

        def get_bool(f):
            return f.get()

        for i in range(len(self.V)):
            newCheckVar = tk.BooleanVar()
            self.CheckVarsV.append(newCheckVar)
            self.CheckVarsV[-1].set(True)
            self.menuStartingNodes.add_checkbutton(label="{}".format(self.V[i]), variable=self.CheckVarsV[-1],
                                                   command=functools.partial(self.check_v, i,
                                                                             functools.partial(get_bool, self.CheckVarPosFlow),
                                                                             functools.partial(get_bool, self.CheckVarPosQ)))

        # Variablen der Checkboxen zu den einzelnen Spalten
        self.CheckVarFp = tk.BooleanVar()  # wenn 'True', zeige Spalte 'f^+'
        self.CheckVarFm = tk.BooleanVar()  # wenn 'True', zeige Spalte 'f^-'
        self.CheckVarQ = tk.BooleanVar()  # wenn 'True', zeige Spalte für Warteschlangenlänge
        self.CheckVarDelta = tk.BooleanVar()  # wenn 'True', zeige Verhältnis 'q' / 'nu'
        self.CheckVarC = tk.BooleanVar()  # wenn 'True', zeige Kantenkosten zu Beginn der aktuellen Phase
        self.CheckVarLabel = tk.BooleanVar()  # wenn 'True', zeige Label des Endknotens zu Beginn der aktuellen Phase
        self.CheckVarFp.set(True)
        self.CheckVarFm.set(True)
        self.CheckVarQ.set(True)
        self.CheckVarDelta.set(True)
        self.CheckVarC.set(True)
        self.CheckVarLabel.set(True)
        self.CheckVars = [self.CheckVarFp, self.CheckVarFm, self.CheckVarQ, self.CheckVarDelta, self.CheckVarC,
                          self.CheckVarLabel]

        # Checkboxen
        self.filemenuC.add_checkbutton(
            label="f^+", variable=self.CheckVarFp, command=functools.partial(self.check_column, 1))
        self.filemenuC.add_checkbutton(
            label="f^-", variable=self.CheckVarFm, command=functools.partial(self.check_column, 2))
        self.filemenuC.add_checkbutton(
            label="q zu Beginn", variable=self.CheckVarQ, command=functools.partial(self.check_column, 3))
        self.filemenuC.add_checkbutton(
            label="Änderung q / nu", variable=self.CheckVarDelta, command=functools.partial(self.check_column, 4))
        self.filemenuC.add_checkbutton(
            label="c zu Beginn", variable=self.CheckVarC, command=functools.partial(self.check_column, 5))
        self.filemenuC.add_checkbutton( label="label des Endknotens zu Beginn", variable=self.CheckVarLabel,
            command=functools.partial(self.check_column, 6))
        self.eps = 10**(-8)
        self.col_heads = ["f^+", "f^-", "q zu Beginn", "Änderung q / nu", "c zu Beginn",
                          "label des Endknotens \n zu Beginn"]
        self.cols = len(self.col_heads)
        self.phase_ind = 0

        # Gitter der Größe <#Kanten + 2> x <#Spalten + 1>. Die ersten zwei Zeilen sind für Überschrift und Trennlinie,
        # die erste Spalte für die Bezeichnungen der Kanten. Der Rest des Gitters der Größe <#Kanten> x <#Spalten> ent-
        # hält alle entsprechenden Werte.
        self.grid_entries = [[None for col in range(self.cols + 1)] for ro in range(self.m + 2)]

        # Beschriftung der ersten zwei Zeilen der ersten Spalte
        self.grid_entries[0][0] = tk.Text(self.table, width=22, height=2)
        self.grid_entries[0][0].insert('end', "Kante")
        self.grid_entries[0][0].config(state='disabled')
        self.grid_entries[0][0].grid(row=0, column=0)
        self.grid_entries[1][0] = tk.Text(self.table, width=22, height=1)
        self.grid_entries[1][0].insert('end', "--------------------------")
        self.grid_entries[1][0].config(state='disabled')
        self.grid_entries[1][0].grid(row=1, column= 0)

        # Beschriftung der ersten zwei Zeilen der Spalten 2 - 'self.cols'+1
        for co in range(1, self.cols + 1):
            self.grid_entries[0][co] = tk.Text(self.table, width=22, height=2)
            self.grid_entries[0][co].insert('end', self.col_heads[co - 1])
            self.grid_entries[0][co].config(state='disabled')
            self.grid_entries[0][co].grid(row=0, column=co)
            self.grid_entries[1][co] = tk.Text(self.table, width=22, height=1)
            self.grid_entries[1][co].insert('end', "--------------------------")
            self.grid_entries[1][co].config(state='disabled')
            self.grid_entries[1][co].grid(row=1, column=co)

        # Beschriftung der Zeilen 3 - 'self.m'+2 für alle Spalten
        for ro in range(self.m):
            grid_ro = ro + 2
            self.grid_entries[grid_ro][0] = tk.Text(self.table, width=22, height=1)
            self.grid_entries[grid_ro][0].insert('end', "({}, {})".format(self.E[ro][0], self.E[ro][1]))
            self.grid_entries[grid_ro][0].config(state='disabled')
            self.grid_entries[grid_ro][0].grid(row=grid_ro, column=0)
            for co in range(1, self.cols + 1):
                self.grid_entries[grid_ro][co] = tk.Text(self.table, width=22, height=1)

            self.grid_entries[grid_ro][1].insert('end', self.fp[ro][0][1])
            self.grid_entries[grid_ro][2].insert('end', 0)
            self.grid_entries[grid_ro][3].insert('end', 0)
            self.grid_entries[grid_ro][4].insert('end', self.change_of_cost(ro, 0, self.fp[ro][0][1]))
            self.grid_entries[grid_ro][5].insert('end', self.c[0][ro])
            self.grid_entries[grid_ro][6].insert('end', self.labels[0][self.V.index(self.E[ro][1])])

            for co in range(1, self.cols + 1):
                self.grid_entries[grid_ro][co].config(state='disabled')
                self.grid_entries[grid_ro][co].grid(row=grid_ro, column=co)

        self.next = tk.Button(self.table, text="Weiter", padx=68, command=self.next)
        self.prev = tk.Button(self.table, text="Zurück", padx=68, command=self.previous)
        self.next.grid(row=self.m+2, column=self.cols)
        self.prev.grid(row=self.m+2, column=0)
        # Hilfsvariablen zur Platzierung des 'Weiter' - Buttons
        self.next_btn_col = self.cols
        self.hidden_cols = []
        # Hilfsvariable zur Verwaltung der Zeilen
        self.hidden_rows = []

        self.table.mainloop()

    def change_of_cost(self, e_ind, phase_ind, z):
        """
        Gibt die momentane Kostenänderung der Kante mit Index 'e_ind' bei Zufluss 'z' an.
        :param e_ind: Index der betrachteten Kante
        :param phase_ind: Index des betrachteten Zeitpunkts
        :param z: Einflussrate in Kante
        :return: Änderungsrate der Kosten
        """
        if self.q[phase_ind][e_ind] > self.eps:
            return (z - self.nu[e_ind])/self.nu[e_ind]
        return np.max([(z - self.nu[e_ind])/self.nu[e_ind], 0])

    def check_all_rows(self):
        """
        Überprüft Status der Checkbox zum Zeilenmenüpunkt 'alle'. Ist die zugehörige Variable 'self.CheckVarAll' auf
        True, so werden alle Zeilen angezeigt. Andernfalls werden durch aufrufen der Funktionen 'self.check_pos_q' und
        'self.check_pos_vol' die Status der anderen Zeilenfilter überprüft.
        :return: kein Rückgabewert
        """
        non_hidden_cols = list(set(range(self.cols + 1)) - set(self.hidden_cols))
        if self.CheckVarAll.get():
            self.CheckVarPosFlow.set(True)
            self.CheckVarPosQ.set(True)
            for i in range(len(self.CheckVarsV)):
                self.CheckVarsV[i].set(True)
            for ro in self.hidden_rows:
                for co in non_hidden_cols:
                    self.grid_entries[ro + 2][co].grid()
            self.hidden_rows = []
        else:
            self.hidden_rows = list(range(self.m))
            for ro in range(2, self.m + 2):
                for co in non_hidden_cols:
                    self.grid_entries[ro][co].grid_remove()
            # Im Fall das genau eine der Variablen 'self.CheckVarPosFlow' und 'self.CheckVarPosQ' 'False' ist, muss die entsprechende
            # Funktion 'self.check_pos_flow', bzw. 'self.check_pos_q' zuerst aufgerufen werden, um Korrektheit zu garantieren (sonst
            # können Zeilen die eigentlich angezeigt werden sollen versteckt werden).
            bool_q = self.CheckVarPosQ.get()
            if self.CheckVarPosFlow.get():
                self.check_pos_q()
                self.check_pos_flow()
                if bool_q:
                    for i in range(self.n):
                        self.check_v(i, check_vol=True, check_q=True)
                else:
                    for i in range(self.n):
                        self.check_v(i, check_vol=True, check_q=False)
            else:
                self.check_pos_flow()
                self.check_pos_q()
                if bool_q:
                    for i in range(self.n):
                        self.check_v(i, check_vol=False, check_q=True)
                else:
                    for i in range(self.n):
                        self.check_v(i, check_vol=False, check_q=False)
        return

    def check_v(self, v_ind, check_vol, check_q):
        """
        Aufruf bei Änderung des Status der Checkbox zum Knoten mit Index 'v_ind'. Ist die zur Checkbox gehörende Variable
        'self.CheckVarsV[v_ind]' "True", so werden im Menü alle Kanten mit diesem Startknoten in der Tabelle angezeigt. Umgekehrt werden
        diese Kanten aus der Tabelle entfernt, falls die Checkbox auf "False" gesetzt ist.
        :param v_ind: Index des betrachteten Knotens
        :param check_vol: Momentaner Wahrheitswert der Checkbox 'self.CheckVarPosFlow' oder Funktion, die diesen Wahrheitswert bestimmt
        :param check_q: Momentaner Wahrheitswert der Checkbox 'self.CheckVarPosQ' oder Funktion, die diesen Wahrheitswert bestimmt
        :return: Kein Rückgabewert
        """
        if self.CheckVarAll.get():
            return
        # 'check_vol', 'check_q' sind hier entweder der Wahrheitswert der entsprechenden Checkbox, oder eine Funktion der Form
        # 'functools.partial', durch deren Aufruf diese Wahrheitswerte bestimmt werden.
        if callable(check_vol):
            check_vol = check_vol()
            check_q = check_q()
        v = self.V[v_ind]
        delta_plus = [self.E.index(e) for e in self.E if e[0] == v]
        non_hidden_cols = list(set(range(self.cols + 1)) - set(self.hidden_cols))
        if self.CheckVarsV[v_ind].get():
            # leer, falls 'self.CheckVarsV[v_ind]' vorher bereits "True"
            newly_visible = list(set(self.hidden_rows).intersection(set(delta_plus)))
            for ro in newly_visible:
                self.hidden_rows.remove(ro)
                for co in non_hidden_cols:
                    self.grid_entries[ro + 2][co].grid()
        elif not (check_vol and self.flow_vol[self.phase_ind][v_ind] > self.eps):
            if check_q:
                # entferne alle Kanten in 'delta_plus' mit Warteschlangenlänge 0
                for e_ind in delta_plus:
                    if self.q[self.phase_ind][e_ind] < self.eps:
                        self.hidden_rows.append(e_ind)
                        for co in non_hidden_cols:
                            self.grid_entries[e_ind + 2][co].grid_remove()
            else:
                # entferne alle Kanten in 'delta_plus'
                for ro in delta_plus:
                    if ro not in self.hidden_rows:
                        self.hidden_rows.append(ro)
                        for co in non_hidden_cols:
                            self.grid_entries[ro + 2][co].grid_remove()
        return

    def check_pos_flow(self):
        """
        Prüft den Status der Variable 'self.CheckVarPosFlow'. Ist dieser 'True', so werden immer alle Kanten angezeigt,
        in deren Startknoten sich eine positive Flussmenge befindet. Ist der Wert 'False', so werden diese Kanten
        versteckt (else - Fall), außer die Variable 'self.CheckVarPosQ' ist 'True' und die entsprechende Kante hat eine
        momentan positive Warteschlangenlänge (elif - Fall).
        :return: kein Rückgabewert
        """
        if self.CheckVarAll.get():
            return
        non_hidden_cols = list(set(range(self.cols + 1)) - set(self.hidden_cols))
        if self.CheckVarPosFlow.get():
            init_rows = self.hidden_rows.copy()
            for ro in init_rows:
                v = self.E[ro][0]
                v_ind = self.V.index(v)
                if self.flow_vol[self.phase_ind][v_ind] > self.eps:
                    self.hidden_rows.remove(ro)
                    for co in non_hidden_cols:
                        self.grid_entries[ro + 2][co].grid()
        elif self.CheckVarPosQ.get():
            # Liste aller Startknoten von Kanten, welche momentan angezeigt werden (ohne Wiederholungen)
            non_hidden_rows = list(set(range(self.m)) - set(self.hidden_rows))
            non_hidden_start_nodes = list(set([self.E[e_ind][0] for e_ind in non_hidden_rows]))
            for v in non_hidden_start_nodes:
                v_ind = self.V.index(v)
                if self.flow_vol[self.phase_ind][v_ind] > self.eps:
                    self.check_v(v_ind, check_vol=False, check_q=True)
        else:
            # Liste aller Startknoten von Kanten, welche momentan angezeigt werden (ohne Wiederholungen)
            non_hidden_rows = list(set(range(self.m)) - set(self.hidden_rows))
            non_hidden_start_nodes = list(set([self.E[e_ind][0] for e_ind in non_hidden_rows]))
            for v in non_hidden_start_nodes:
                v_ind = self.V.index(v)
                if self.flow_vol[self.phase_ind][v_ind] > self.eps:
                    self.check_v(v_ind, check_vol=False, check_q=False)
        return

    def check_pos_q(self):
        """
        Prüft den Status der Variable 'self.CheckVarPosQ'. Ist dieser 'True', so werden immer alle Kanten angezeigt, die
        eine positive Warteschlangenlänge besitzen. Ist der Wert 'False', so werden diese Kanten versteckt
        (else - Fall), außer die Variable 'self.CheckVarPosFlow' ist 'True' und die entsprechende Kante hat eine
        momentan positive Flussmenge in ihrem Startknoten (elif - Fall).
        :return: kein Rückgabewert
        """
        if self.CheckVarAll.get():
            return
        non_hidden_cols = list(set(range(self.cols + 1)) - set(self.hidden_cols))
        if self.CheckVarPosQ.get():
            init_rows = self.hidden_rows.copy()
            for ro in init_rows:
                if self.q[self.phase_ind][ro] > self.eps:
                    self.hidden_rows.remove(ro)
                    for co in non_hidden_cols:
                        self.grid_entries[ro + 2][co].grid()
        else:
            non_hidden_rows = list(set(range(self.m)) - set(self.hidden_rows))
            # Liste aller Startknoten von Kanten, welche momentan angezeigt werden und positive Warteschlangenlänge haben (ohne
            # Wiederholungen)
            start_nodes_pos_q = list(set([self.E[e_ind][0] for e_ind in non_hidden_rows if self.q[self.phase_ind][e_ind] > self.eps]))
            if self.CheckVarPosFlow.get():
                for v in start_nodes_pos_q:
                    self.check_v(self.V.index(v), check_vol=True, check_q=False)
            else:
                for v in start_nodes_pos_q:
                    self.check_v(self.V.index(v), check_vol=False, check_q=False)
        return

    def check_column(self, col_ind):
        """
        Wird aufgerufen bei Änderung des Status einer der Checkboxen im Menü 'Spalten'. Ist die Checkbox zur Spalte
        'col_ind' aktiv, so wird die entsprechende Spalte angezeigt, andernfalls wird die Spalte versteckt.
        :param col_ind: Index der Spalte im grid - layout der Tabelle
        :return: kein Rückgabewert
        """
        if self.CheckVars[col_ind - 1].get():
            self.hidden_cols.remove(col_ind)
            self.grid_entries[0][col_ind].grid()
            self.grid_entries[1][col_ind].grid()
            if self.CheckVarAll.get():
                for ro in range(self.m + 2):
                    self.grid_entries[ro][col_ind].grid()
            else:
                self.check_all_rows()
            if col_ind > self.next_btn_col:
                self.next_btn_col = col_ind
                self.next.grid(row=self.m+2, column=self.next_btn_col)
        else:
            for ro in range(self.m + 2):
                self.grid_entries[ro][col_ind].grid_remove()
            if col_ind == self.next_btn_col:
                leftover = set(range(1, col_ind)) - set(self.hidden_cols)
                if len(leftover) > 0:
                    self.next_btn_col = np.max(list(leftover))
                else:
                    self.next_btn_col = 1
                self.next.grid(row=self.m+2, column=self.next_btn_col)
            self.hidden_cols.append(col_ind)
        return

    def next(self):
        """
        wird aufgerufen beim Klicken des Benutzers auf den Button 'Weiter'. Aktualisiert die Tabelle und zeigt die Werte
        des nächsten Zeitpunkts an
        :return: kein Rückgabewert
        """
        if self.phase_ind == len(self.phases) - 1:
            return
        self.phase_ind += 1
        theta = self.phases[self.phase_ind]
        self.table.title("Zeit: {}".format(theta))
        for ro in range(self.m):
            grid_ro = ro + 2
            for co in [3, 5, 6]:
                entry = self.grid_entries[grid_ro][co]
                entry.config(state='normal')
                entry.delete(1.0, tk.END)

            if self.phase_ind in self.fp_ind[ro]:
                ind = self.fp_ind[ro].index(self.phase_ind)
                if self.fp[ro][0] == (0, 0):
                    ind += 1
                change = self.change_of_cost(ro, self.phase_ind, self.fp[ro][ind][1])
                fp_entry = self.grid_entries[grid_ro][1]
                fp_entry.config(state='normal')
                fp_entry.delete(1.0, tk.END)
                fp_entry.insert('end', self.fp[ro][ind][1])
                q_entry = self.grid_entries[grid_ro][4]
                q_entry.config(state='normal')
                q_entry.delete(1.0, tk.END)
                q_entry.insert('end', change)
            elif abs(self.q[self.phase_ind][ro]) < self.eps:
                q_entry = self.grid_entries[grid_ro][4]
                q_entry.config(state='normal')
                q_entry.delete(1.0, tk.END)
                q_entry.insert('end', 0)

            fm_times = [t for (t, v) in self.fm[ro]]
            if theta in fm_times:
                fm_entry = self.grid_entries[grid_ro][2]
                fm_entry.config(state='normal')
                fm_entry.delete(1.0, tk.END)
                fm_entry.insert('end', self.fm[ro][fm_times.index(theta)][1])

            self.grid_entries[grid_ro][3].insert('end', self.q[self.phase_ind][ro])
            self.grid_entries[grid_ro][5].insert('end', self.c[self.phase_ind][ro])
            self.grid_entries[grid_ro][6].insert('end', self.labels[self.phase_ind][self.V.index(self.E[ro][1])])

            for co in range(1, self.cols + 1):
                self.grid_entries[grid_ro][co].config(state='disable')
            self.check_all_rows()
        return

    def previous(self):
        """
        wird aufgerufen beim Klicken des Benutzers auf den Button 'Zurück'. Aktualisiert die Tabelle und zeigt die Werte
        des vorherigen Zeitpunkts an
        :return: kein Rückgabewert
        """
        if self.phase_ind == 0:
            return
        old_theta = self.phases[self.phase_ind]
        self.phase_ind -= 1
        theta = self.phases[self.phase_ind]
        self.table.title("Zeit: {}".format(theta))
        for ro in range(self.m):
            grid_ro = ro + 2
            for co in [3, 5, 6]:
                entry = self.grid_entries[grid_ro][co]
                entry.config(state='normal')
                entry.delete(1.0, tk.END)

            if self.phase_ind + 1 in self.fp_ind[ro]:
                ind = self.fp_ind[ro].index(self.phase_ind + 1) - 1
                if self.fp[ro][0] == (0, 0):
                    ind += 1
                fp_entry = self.grid_entries[grid_ro][1]
                fp_entry.config(state='normal')
                fp_entry.delete(1.0, tk.END)
                fp_entry.insert('end', self.fp[ro][ind][1])
                q_entry = self.grid_entries[grid_ro][4]
                q_entry.config(state='normal')
                q_entry.delete(1.0, tk.END)
                q_entry.insert('end', self.change_of_cost(ro, self.phase_ind, self.fp[ro][ind][1]))
            elif abs(self.q[self.phase_ind][ro]) > self.eps > abs(self.q[self.phase_ind + 1][ro]):
                q_entry = self.grid_entries[grid_ro][4]
                q_entry.config(state='normal')
                q_entry.delete(1.0, tk.END)
                # berechne Rate, mit der Warteschlange abgebaut wird (abhängig von Einfluss < 'self.nu[ro]')
                dq_dnu = - self.q[self.phase_ind][ro] / ((old_theta - theta) * self.nu[ro])
                if -1 - self.eps < dq_dnu < -1 + self.eps:
                    dq_dnu = -1.0
                q_entry.insert('end', dq_dnu)

            fm_times = [t for (t, v) in self.fm[ro]]
            if theta in fm_times:
                fm_entry = self.grid_entries[grid_ro][2]
                fm_entry.config(state='normal')
                fm_entry.delete(1.0, tk.END)
                fm_entry.insert('end', self.fm[ro][fm_times.index(theta)][1])
            elif old_theta in fm_times:
                fm_entry = self.grid_entries[grid_ro][2]
                fm_entry.config(state='normal')
                fm_entry.delete(1.0, tk.END)
                fm_entry.insert('end', self.fm[ro][fm_times.index(old_theta) - 1][1])

            self.grid_entries[grid_ro][3].insert('end', self.q[self.phase_ind][ro])
            self.grid_entries[grid_ro][5].insert('end', self.c[self.phase_ind][ro])
            self.grid_entries[grid_ro][6].insert('end', self.labels[self.phase_ind][self.V.index(self.E[ro][1])])

            for co in range(1, self.cols + 1):
                self.grid_entries[grid_ro][co].config(state='disable')
            self.check_all_rows()
        return
