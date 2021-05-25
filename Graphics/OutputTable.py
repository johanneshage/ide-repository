import tkinter as tk
import numpy as np

class OutputTable(object):

    def __init__(self, V, E, nu, fp, fp_ind, fm, q, phases, end, c, labels):
        self.V = V
        self.E = E
        self.m = len(E)
        self.q = q
        self.phases = phases
        self.nu = nu
        self.fp = fp
        self.fp_ind = fp_ind
        self.fm = fm
        self.end = end
        self.c = c
        self.labels = labels
        self.table = tk.Tk()
        self.table.title("Zeit: 0")
        self.rows = self.m
        self.col_heads = ["f^+", "f^-", "q zu Beginn", "Änderung q / nu", "c zu Beginn",
                          "label des Endknotens \n zu Beginn"]
        self.cols = len(self.col_heads)
        self.phase_ind = 0

        first_entry = tk.Text(self.table, width=22, height=2)
        first_entry.insert('end', "Kante")
        first_entry.config(state='disabled')
        first_entry.grid(row=0, column=0)
        first_newline = tk.Text(self.table, width=22, height=1)
        first_newline.insert('end', "--------------------------")
        first_newline.config(state='disabled')
        first_newline.grid(row=1, column= 0)

        for co in range(self.cols):
            entry = tk.Text(self.table, width=22, height=2)
            entry.insert('end', self.col_heads[co])
            entry.config(state='disabled')
            entry.grid(row=0, column=co+1)
            newline = tk.Text(self.table, width=22, height=1)
            newline.insert('end', "--------------------------")
            newline.config(state='disabled')
            newline.grid(row=1, column=co+1)

        self.grid_values = [[None for col in range(self.cols)] for ro in range(self.m)]
        for ro in range(self.m):
            first = tk.Text(self.table, width=22, height=1)
            first.insert('end', "({}, {})".format(self.E[ro][0], self.E[ro][1]))
            first.config(state='disabled')
            first.grid(row=ro+2, column=0)
            for co in range(self.cols):
                self.grid_values[ro][co] = tk.Text(self.table, width=22, height=1)

            self.grid_values[ro][0].insert('end', self.fp[ro][0][1])
            self.grid_values[ro][1].insert('end', 0)
            self.grid_values[ro][2].insert('end', 0)
            self.grid_values[ro][3].insert('end', self.change_of_cost(self.E[ro], 0, self.fp[ro][0][1]))
            self.grid_values[ro][4].insert('end', self.c[0][ro])
            self.grid_values[ro][5].insert('end', self.labels[0][self.V.index(self.E[ro][1])])

            for co in range(self.cols):
                self.grid_values[ro][co].config(state='disabled')
                self.grid_values[ro][co].grid(row=ro+2, column=co+1)

        self.next = tk.Button(self.table, text="Weiter", padx=68, command=self.next)
        self.prev = tk.Button(self.table, text="Zurück", padx=68)  # , command=self.previous)
        self.next.grid(row=self.m+2, column=self.cols)
        self.prev.grid(row=self.m+2, column=0)

        self.table.mainloop()

    def change_of_cost(self, e, phase_ind, z):
        """
        Gibt die momentane Kostenänderung der Kante 'e' bei Zufluss 'z' an.
        :param e: Betrachtete Kante
        :param phase_ind: Index des betrachteten Zeitpunkts
        :param z: Einflussrate in Kante 'e'
        :return: Änderungsrate der Kosten
        """
        e_ind = self.E.index(e)
        if self.q[phase_ind][e_ind] > 0:
            return (z - self.nu[e_ind])/self.nu[e_ind]
        return np.max([(z - self.nu[e_ind])/self.nu[e_ind], 0])

    def next(self):
        if self.phase_ind == len(self.phases) - 1:
            return
        self.phase_ind += 1
        theta = self.phases[self.phase_ind]
        self.table.title("Zeit: {}".format(theta))
        for ro in range(self.m):
            for co in [2, 4, 5]:
                entry = self.grid_values[ro][co]
                entry.config(state='normal')
                entry.delete(1.0, tk.END)

            if self.phase_ind in self.fp_ind[ro]:
                ind = self.fp_ind[ro].index(self.phase_ind)
                if self.fp[ro][0] == (0, 0):
                    ind += 1
                fp_entry = self.grid_values[ro][0]
                fp_entry.config(state='normal')
                fp_entry.delete(1.0, tk.END)
                fp_entry.insert('end', self.fp[ro][ind][1])
                q_entry = self.grid_values[ro][3]
                q_entry.config(state='normal')
                q_entry.delete(1.0, tk.END)
                q_entry.insert('end', self.change_of_cost(self.E[ro], self.phase_ind, self.fp[ro][ind][1]))

            fm_times = [t for (t, v) in self.fm[ro]]
            print("FM TIMES", fm_times)
            if theta in fm_times:
                fm_entry = self.grid_values[ro][1]
                fm_entry.config(state='normal')
                fm_entry.delete(1.0, tk.END)
                fm_entry.insert('end', self.fm[ro][fm_times.index(theta)][1])

            self.grid_values[ro][2].insert('end', self.q[self.phase_ind][ro])
            self.grid_values[ro][4].insert('end', self.c[self.phase_ind][ro])
            self.grid_values[ro][5].insert('end', self.labels[self.phase_ind][self.V.index(self.E[ro][1])])

            for co in range(self.cols):
                self.grid_values[ro][co].config(state='disable')
        return
