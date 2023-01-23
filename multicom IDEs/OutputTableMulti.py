import tkinter as tk
import numpy as np
import functools


class OutputTableMulti(object):

    def __init__(self, G, V, E, I, nu, fp, fp_ind, fm, q, phases, c, labels, vol):

        self.G = G
        self.V = V  # Knotenmenge
        self.E = E  # Kantenmenge
        self.I = I  # Anzahl Güter
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
        # self.menu.option_add('*tearOff', False)
        self.table.config(menu=self.menu)
        self.table.title("Zeit: 0")
        self.filemenuR = tk.Menu(self.menu)  # Menüpunkt Zeilen
        self.filemenuC = tk.Menu(self.menu)  # Menüpunkt Spalten
        self.menu.add_cascade(label="Zeilen", menu=self.filemenuR)
        self.menu.add_cascade(label="Spalten", menu=self.filemenuC)

        # Variablen der Checkboxen zu den Zeilengruppen
        # self.CheckVarAll = tk.BooleanVar()  # wenn 'True', zeige alle Zeilen
        # self.CheckVarPosFlow = tk.BooleanVar()  # wenn 'True', zeige alle Zeilen mit momentan positivem Flussvolumen im Startknoten der entsprechenden Kante
        # self.CheckVarPosQ = tk.BooleanVar()  # wenn 'True', zeige alle Zeilen mit momentan positiver Warteschlangenlänge
        # self.CheckVarAll.set(True)
        # self.CheckVarPosFlow.set(True)
        # self.CheckVarPosQ.set(True)

        # self.filemenuR.add_checkbutton(label="alle", variable=self.CheckVarAll, command=self.check_all_rows)
        # self.filemenuR.add_checkbutton(label="Kanten mit positivem Flussvolumen im Startknoten", variable=self.CheckVarPosFlow, command=self.check_pos_flow)
        # self.filemenuR.add_checkbutton(label="Kanten mit positiver Warteschlangenlänge", variable=self.CheckVarPosQ, command=self.check_pos_q)

        self.singles = tk.Menu(self.filemenuR)
        self.filemenuR.add_cascade(label="Kanten mit Startknoten...", menu=self.singles)
        self.CheckVarNodes = []
        for (v_ind, v) in enumerate(self.V):
            self.CheckVarNodes.append(tk.BooleanVar())
            self.CheckVarNodes[v_ind].set(True)
            self.singles.add_checkbutton(label="{}".format(v), variable=self.CheckVarNodes[v_ind], command=functools.partial(self.check_node, v_ind))

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
        self.CheckVars = [self.CheckVarFp, self.CheckVarFm, self.CheckVarQ, self.CheckVarDelta, self.CheckVarC, self.CheckVarLabel]

        # Checkboxen
        self.filemenuC.add_checkbutton(label="f^+", variable=self.CheckVarFp, command=functools.partial(self.check_column, 2))
        self.filemenuC.add_checkbutton(label="f^-", variable=self.CheckVarFm, command=functools.partial(self.check_column, 3))
        self.filemenuC.add_checkbutton(label="q zu Beginn", variable=self.CheckVarQ, command=functools.partial(self.check_column, 4))
        self.filemenuC.add_checkbutton(label="Änderung q / nu", variable=self.CheckVarDelta, command=functools.partial(self.check_column, 5))
        self.filemenuC.add_checkbutton(label="c zu Beginn", variable=self.CheckVarC, command=functools.partial(self.check_column, 6))
        self.filemenuC.add_checkbutton( label="label des Endknotens zu Beginn", variable=self.CheckVarLabel, command=functools.partial(self.check_column, 7))
        self.eps = 10**(-12)
        self.bd_tol = 10**(-6)
        self.col_heads = ["Gut", "f^+", "f^-", "q zu Beginn", "Änderung q / nu", "c zu Beginn", "label des Endknotens \n zu Beginn"]
        self.cols = len(self.col_heads)
        self.single_cols = [4, 5, 6]
        self.multi_cols = [1, 2, 3, 7]
        self.phase_ind = 0

        # Gitter der Größe <#Kanten + 2> x <#Spalten + 1>. Die ersten zwei Zeilen sind für Überschrift und Trennlinie,
        # die erste Spalte für die Bezeichnungen der Kanten. Der Rest des Gitters der Größe <#Kanten> x <#Spalten> ent-
        # hält alle entsprechenden Werte.
        self.grid_entries = np.full((self.m + 2, self.cols + 1), None)

        # Beschriftung der ersten zwei Zeilen der ersten Spalte
        self.grid_entries[0, 0] = tk.Text(self.table, width=16, height=2)
        self.grid_entries[0, 0].insert('end', "Kante")
        self.grid_entries[0, 0].config(state='disabled')
        self.grid_entries[0, 0].grid(row=0, column=0)
        self.grid_entries[1, 0] = tk.Text(self.table, width=16, height=1)
        self.grid_entries[1, 0].insert('end', "--------------------------")
        self.grid_entries[1, 0].config(state='disabled')
        self.grid_entries[1, 0].grid(row=1, column= 0)

        # Beschriftung der ersten zwei Zeilen der Spalten 2 - 'self.cols'+1
        for co in range(1, self.cols + 1):
            self.grid_entries[0, co] = tk.Text(self.table, width=22, height=2)
            self.grid_entries[0, co].insert('end', self.col_heads[co - 1])
            self.grid_entries[0, co].config(state='disabled')
            self.grid_entries[0, co].grid(row=0, column=co)
            self.grid_entries[1, co] = tk.Text(self.table, width=22, height=1)
            self.grid_entries[1, co].insert('end', "--------------------------")
            self.grid_entries[1, co].config(state='disabled')
            self.grid_entries[1, co].grid(row=1, column=co)

        self.frame_entries = np.full((self.m, len(self.multi_cols), self.I), None)
        # Beschriftung der Zeilen 3 - 'self.m'+2 für alle Spalten
        for ro in range(self.m):
            grid_ro = ro + 2
            if grid_ro % 2 == 0:
                bg_color = "lightgrey"
            else:
                bg_color = "white"
            self.grid_entries[grid_ro, 0] = tk.Text(self.table, width=16, height=self.I, bg=bg_color)
            self.grid_entries[grid_ro, 0].insert('end', "({}, {})".format(self.E[ro][0], self.E[ro][1]))
            self.grid_entries[grid_ro, 0].config(state='disabled')
            self.grid_entries[grid_ro, 0].grid(row=grid_ro, column=0, sticky="ns")
            for co in self.single_cols:
                self.grid_entries[grid_ro, co] = tk.Text(self.table, width=22, height=self.I, bg=bg_color)
            len_multi_cols = len(self.multi_cols)
            for ind in range(len_multi_cols):
                co = self.multi_cols[ind]
                self.grid_entries[grid_ro, co] = tk.Frame(self.table, bg=bg_color)
                for i in range(self.I):
                    self.frame_entries[ro, ind, i] = tk.Text(self.grid_entries[grid_ro, co], width=22, height=1, bg=bg_color)
                    self.frame_entries[ro, ind, i].grid(row=i, column=0)

            self.grid_entries[grid_ro, 4].insert('end', 0)
            total_inflow = 0
            for i in range(self.I):
                total_inflow += self.fp[i][ro][0][1]
            self.grid_entries[grid_ro, 5].insert('end', self.change_of_cost(ro, 0, total_inflow))
            self.grid_entries[grid_ro, 6].insert('end', self.c[0][ro])

            for i in range(self.I):
                self.frame_entries[ro, 0, i].insert('end', i)
                self.frame_entries[ro, 0, i].config(state='disabled')
                self.frame_entries[ro, 1, i].insert('end', self.fp[i][ro][0][1])
                self.frame_entries[ro, 1, i].config(state='disabled')
                self.frame_entries[ro, 2, i].insert('end', 0)
                self.frame_entries[ro, 2, i].config(state='disabled')
                self.frame_entries[ro, 3, i].insert('end', self.labels[i][0][self.V.index(self.E[ro][1])])
                self.frame_entries[ro, 3, i].config(state='disabled')

            for co in self.single_cols:
                self.grid_entries[grid_ro, co].config(state='disabled')

            for co in range(1, self.cols + 1):
                self.grid_entries[grid_ro, co].grid(row=grid_ro, column=co, sticky="ns")

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
        Gibt die momentane Kostenänderung der Kante mit Index 'e_ind' bei (Gesamt-)Zufluss 'z' an.
        :param e_ind: Index der betrachteten Kante
        :param phase_ind: Index des betrachteten Zeitpunkts
        :param z: Einflussrate in Kante
        :return: Änderungsrate der Kosten
        """
        dif = z - self.nu[e_ind]
        if abs(dif) < self.bd_tol:
            return 0
        if self.q[phase_ind][e_ind] > 0:
            return dif / self.nu[e_ind]
        return np.max([dif / self.nu[e_ind], 0])

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
            for ro in self.hidden_rows:
                for co in non_hidden_cols:
                    self.grid_entries[ro + 2, co].grid()
            self.hidden_rows = []
        else:
            self.hidden_rows = list(range(self.m))
            for ro in range(self.m):
                for co in non_hidden_cols:
                    self.grid_entries[ro + 2, co].grid_remove()
            # Im Fall dass genau eine der Variablen 'self.CheckVarPosFlow' und 'self.CheckVarPosQ' 'False' ist, muss
            # die entsprechende Funktion 'self.check_pos_flow', bzw. 'self.check_pos_q' zuerst aufgerufen werden,
            # um Korrektheit zu garantieren (sonst können Zeilen die eigentlich angezeigt werden sollen versteckt
            # werden).
            if self.CheckVarPosFlow.get():
                self.check_pos_q()
                self.check_pos_flow()
            else:
                self.check_pos_flow()
                self.check_pos_q()
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
                if self.flow_vol[self.phase_ind][v_ind] > 0:
                    self.hidden_rows.remove(ro)
                    for co in non_hidden_cols:
                        self.grid_entries[ro + 2, co].grid()
        elif self.CheckVarPosQ.get():
            non_hidden_rows = list(set(range(self.m)) - set(self.hidden_rows))
            for ro in non_hidden_rows:
                v = self.E[ro][0]
                v_ind = self.V.index(v)
                if self.flow_vol[self.phase_ind][v_ind] > 0 and self.q[self.phase_ind][ro] < self.eps:
                    self.hidden_rows.append(ro)
                    for co in non_hidden_cols:
                        self.grid_entries[ro + 2, co].grid_remove()
        else:
            non_hidden_rows = list(set(range(self.m)) - set(self.hidden_rows))
            for ro in non_hidden_rows:
                v = self.E[ro][0]
                v_ind = self.V.index(v)
                if self.flow_vol[self.phase_ind][v_ind] > 0:
                    self.hidden_rows.append(ro)
                    for co in non_hidden_cols:
                        self.grid_entries[ro + 2, co].grid_remove()
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
                if self.q[self.phase_ind][ro] > 0:
                    self.hidden_rows.remove(ro)
                    for co in non_hidden_cols:
                        self.grid_entries[ro + 2, co].grid()
        elif self.CheckVarPosFlow.get():
            non_hidden_rows = list(set(range(self.m)) - set(self.hidden_rows))
            for ro in non_hidden_rows:
                v_ind = self.V.index(self.E[ro][0])
                if self.q[self.phase_ind][ro] > 0 and self.flow_vol[self.phase_ind][v_ind] < self.eps:
                    self.hidden_rows.append(ro)
                    for co in non_hidden_cols:
                        self.grid_entries[ro + 2, co].grid_remove()
        else:
            non_hidden_rows = list(set(range(self.m)) - set(self.hidden_rows))
            for ro in non_hidden_rows:
                if self.q[self.phase_ind][ro] > 0:
                    self.hidden_rows.append(ro)
                    for co in non_hidden_cols:
                        self.grid_entries[ro + 2, co].grid_remove()
        return

    def check_node(self, v_ind):
        """

        :param v_ind:
        :return:
        """
        non_hidden_cols = list(set(range(self.cols + 1)) - set(self.hidden_cols))
        if self.CheckVarNodes[v_ind].get():
            v = self.V[v_ind]
            delta_p = [self.E.index((v,u)) for u in self.G[v].keys()]
            for ro in delta_p:
                for co in non_hidden_cols:
                    self.grid_entries[ro + 2, co].grid()
        else:
            v = self.V[v_ind]
            delta_p = [self.E.index((v,u)) for u in self.G[v].keys()]
            for ro in delta_p:
                for co in non_hidden_cols:
                    self.grid_entries[ro + 2, co].grid_remove()

    def check_column(self, col_ind):
        """
        Wird aufgerufen bei Änderung des Status einer der Checkboxen im Menü 'Spalten'. Ist die Checkbox zur Spalte
        'col_ind' aktiv, so wird die entsprechende Spalte angezeigt, andernfalls wird die Spalte versteckt.
        :param col_ind: Index der Spalte im grid - layout der Tabelle
        :return: kein Rückgabewert
        """
        if self.CheckVars[col_ind - 2].get():
            self.hidden_cols.remove(col_ind)
            self.grid_entries[0, col_ind].grid()
            self.grid_entries[1, col_ind].grid()
            if self.CheckVarAll.get():
                if col_ind in self.multi_cols:
                    self.grid_entries[0, col_ind].grid()
                    self.grid_entries[1, col_ind].grid()
                for ro in range(self.m + 2):
                    self.grid_entries[ro, col_ind].grid()
            else:
                self.check_all_rows()
            if col_ind > self.next_btn_col:
                self.next_btn_col = col_ind
                self.next.grid(row=self.m+2, column=self.next_btn_col)
        else:
            if col_ind in self.multi_cols:
                self.grid_entries[0, col_ind].grid_remove()
                self.grid_entries[1, col_ind].grid_remove()
            for ro in range(self.m + 2):
                self.grid_entries[ro, col_ind].grid_remove()
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
            for co in [4, 6]:
                entry = self.grid_entries[grid_ro, co]
                entry.config(state='normal')
                entry.delete(1.0, tk.END)

            multi_ind = self.multi_cols.index(7)
            for i in range(self.I):
                self.frame_entries[ro, multi_ind, i].config(state='normal')
                self.frame_entries[ro, multi_ind, i].delete(1.0, tk.END)

            # fp_frame = self.frame_entries[ro, self.multi_cols.index(2), :]
            fp_frame = self.frame_entries[ro, 1, :]
            q_delta = self.grid_entries[grid_ro, 5]
            fp_sum = 0
            fp_changed = False
            for i in range(self.I):
                if self.phase_ind in self.fp_ind[i][ro]:
                    fp_changed = True
                    ind = self.fp_ind[i][ro].index(self.phase_ind)
                    if self.fp[i][ro][0] == (0, 0):
                        ind += 1
                    fp_sum += self.fp[i][ro][ind][1]
                    fp_frame[i].config(state='normal')
                    fp_frame[i].delete(1.0, tk.END)
                    fp_frame[i].insert('end', self.fp[i][ro][ind][1])

                fm_times = [t for (t, v) in self.fm[i][ro]]
                fm_where = np.where([0 <= theta - t < self.bd_tol * 10 and not (0 <= self.phases[self.phase_ind-1] - t < self.bd_tol * 10) for t in fm_times])[0]
                for fm_ind in fm_where:
                    fm_entry = self.frame_entries[ro, self.multi_cols.index(3), i]
                    # fm_entry = self.grid_entries[grid_ro, 3][i]
                    fm_entry.config(state='normal')
                    fm_entry.delete(1.0, tk.END)
                    fm_entry.insert('end', self.fm[i][ro][fm_ind][1])

            if fp_changed:
                change = self.change_of_cost(ro, self.phase_ind, fp_sum)
                q_delta.config(state='normal')
                q_delta.delete(1.0, tk.END)
                q_delta.insert('end', change)
            elif abs(self.q[self.phase_ind][ro]) < self.eps:
                q_delta.config(state='normal')
                q_delta.delete(1.0, tk.END)
                q_delta.insert('end', 0)

            self.grid_entries[grid_ro, 4].insert('end', self.q[self.phase_ind][ro])
            self.grid_entries[grid_ro, 6].insert('end', self.c[self.phase_ind][ro])
            for i in range(self.I):
                self.frame_entries[ro, 3, i].insert('end', self.labels[i][self.phase_ind][self.V.index(self.E[ro][1])])

            for co in self.single_cols:
                self.grid_entries[grid_ro, co].config(state='disable')
            for co in self.multi_cols:
                multi_ind = self.multi_cols.index(co)
                for i in range(self.I):
                    self.frame_entries[ro, multi_ind, i].config(state='disable')

            # self.check_all_rows()
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
            for co in [4, 6]:
                entry = self.grid_entries[grid_ro, co]
                entry.config(state='normal')
                entry.delete(1.0, tk.END)

            multi_ind = self.multi_cols.index(7)
            for i in range(self.I):
                self.frame_entries[ro, multi_ind, i].config(state='normal')
                self.frame_entries[ro, multi_ind, i].delete(1.0, tk.END)

            # fp_frame = self.frame_entries[ro, self.multi_cols.index(3), :]
            fp_frame = self.frame_entries[ro, 1, :]
            q_delta = self.grid_entries[grid_ro, 5]
            fp_sum = 0
            fp_changed = False
            for i in range(self.I):
                if self.phase_ind + 1 in self.fp_ind[i][ro]:
                    fp_changed = True
                    ind = self.fp_ind[i][ro].index(self.phase_ind + 1) - 1
                    if self.fp[i][ro][0] == (0, 0):
                        ind += 1
                    fp_sum += self.fp[i][ro][ind][1]
                    fp_frame[i].config(state='normal')
                    fp_frame[i].delete(1.0, tk.END)
                    fp_frame[i].insert('end', self.fp[i][ro][ind][1])

                fm_times = [t for (t, v) in self.fm[i][ro]]
                fm_where = np.where([0 <= theta - t < self.bd_tol *10 and (self.phase_ind == 0 or not (0 <= self.phases[self.phase_ind-1] < self.bd_tol *10)) for t in fm_times])[0]
                fm_where_old = np.where([0 < old_theta - t < self.bd_tol * 10 and not (0 <= theta - t < self.bd_tol * 10) for t in fm_times])[0]
                for fm_ind in fm_where:
                    fm_entry = self.frame_entries[ro, self.multi_cols.index(3), i]
                    # fm_entry = self.grid_entries[grid_ro, 3][i]
                    fm_entry.config(state='normal')
                    fm_entry.delete(1.0, tk.END)
                    fm_entry.insert('end', self.fm[i][ro][fm_ind][1])
                for fm_ind in fm_where_old:
                    fm_entry = self.frame_entries[ro, self.multi_cols.index(3), i]
                    # fm_entry = self.grid_entries[grid_ro, 3][i]
                    fm_entry.config(state='normal')
                    fm_entry.delete(1.0, tk.END)
                    fm_entry.insert('end', self.fm[i][ro][fm_ind - 1][1])

            if fp_changed:
                q_delta.config(state='normal')
                q_delta.delete(1.0, tk.END)
                q_delta.insert('end', self.change_of_cost(ro, self.phase_ind, fp_sum))
            elif abs(self.q[self.phase_ind][ro]) > self.eps > abs(self.q[self.phase_ind + 1][ro]):
                # berechne Rate, mit der Warteschlange abgebaut wird (abhängig von Einfluss < 'self.nu[ro]')
                dq_dnu = - self.q[self.phase_ind][ro] / ((old_theta - theta) * self.nu[ro])
                if -1 - self.eps < dq_dnu < -1 + self.eps:
                    dq_dnu = -1.0
                q_delta.config(state='normal')
                q_delta.delete(1.0, tk.END)
                q_delta.insert('end', dq_dnu)

            self.grid_entries[grid_ro, 4].insert('end', self.q[self.phase_ind][ro])
            self.grid_entries[grid_ro, 6].insert('end', self.c[self.phase_ind][ro])
            for i in range(self.I):
                self.frame_entries[ro, 3, i].insert('end', self.labels[i][self.phase_ind][self.V.index(self.E[ro][1])])

            for co in self.single_cols:
                self.grid_entries[grid_ro, co].config(state='disable')
            len_multi_cols = len(self.multi_cols)
            for co in range(len_multi_cols):
                for i in range(self.I):
                    self.frame_entries[ro, co, i].config(state='disable')
            # self.check_all_rows()
        return
