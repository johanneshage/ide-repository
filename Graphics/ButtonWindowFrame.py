import tkinter as tk
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import networkx as nx
from collections import deque
from matplotlib.widgets import Button
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import sys


class ButtonWindowFrame(tk.Frame):
    """
    Diese Klasse verwaltet die Button-Leiste am unteren Bidlschirmrand, mit den Buttons:
    "Zurück", "Pause", bzw. "Fortsetzen", und "Weiter"
    """

    def __init__(self, E, V, I, ST, r, posit):
        """
        Erzeugt Buttons und Layout der Button-Leiste
        :param app: Objekt der Klasse "Application", wird übergeben, damit durch Klicken der Buttons auf die Haupt-
         anwendung zugegriffen werden kann
        """
        self.E = E.copy()
        self.V = V.copy()
        self.I = I.copy()
        self.ST = ST.copy()
        self.r = r.copy()
        self.master = tk.Tk()
        self.fig = plt.figure(figsize=(3,3), dpi=100)
        self.ax = self.fig.add_subplot(111)
        plt.ylim(-1.5, 1.5)  # setzte Höhe des Plots, damit Rechtecke in jedem Plot gleich groß
        plt.xlim(-1.25, 1.25)
        #self.app = app
        self.screenwidth = self.master.winfo_screenwidth()
        self.screenheight = self.master.winfo_screenheight()
        #self.master.wm_attributes('-topmost', True)
        #self.master.overrideredirect(True)
        self.master.geometry("{}x{}".format(int(self.screenwidth), int(self.screenheight)))
        #self.master.config(background='black')
        #self.master.update_idletasks()
        #self.master.config(width = self.master.winfo_screenwidth(), height = self.master.winfo_screenheight()/6) #,
        # background="#0000ee"
        super().__init__(self.master)
        self.master.protocol('WM_DELETE_WINDOW', sys.exit)  # beendet Programm bei Klicken des 'X'-Buttons
        #self.pack()
        #self.pack(fill="both", expand=True)
        canvas = FigureCanvasTkAgg(self.fig, master=self.master)
        #canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)


        self.colorMap = []
        # Farben werden für unterschiedliche Zielknoten verwendet (gibt es mehr Zielknoten als "len(self.colors)", so
        # werden die Farben mehrfach verwendet)
        self.colors = ["limegreen", "orange", "magenta", "turquoise", "gold", "darkviolet"]
        self.graph = nx.DiGraph()  # erzeuge networkx-Graph
        self.posit = posit
        for v in self.V:
            self.graph.add_node(v)
            name = str(v)
            if (name.startswith("s") or name.startswith("t")) and len(name) > 1:
                try:
                    # bestimme Nummer des Start-/Zielknotens; ValueError, falls keine gültige Nummer vorhanden
                    farbe = int(name[1:])
                    if name.startswith("s"):
                        # weise Quellknoten die jeweilige Farbe zu
                        self.colorMap.append("deepskyblue")
                        continue
                    else:
                        # weise Zielknoten die jeweilige Farbe zu
                        self.colorMap.append(self.colors[(farbe-1) % len(self.colors)])
                        continue
                except ValueError:
                    # Knoten ist kein Start- oder Zielknoten, sondern ein anderer Knoten dessen Name mit "s" oder "t"
                    # beginnt
                    self.colorMap.append("paleturquoise")  # Standardfarbe für alle anderen Knoten
                    continue
            self.colorMap.append("paleturquoise")  # Standardfarbe für alle anderen Knoten
        for e in self.E:
            self.graph.add_edge(e)

        self.spieler = []  # Rechtecke für Spieler
        self.numbers = []  # Nummerierung der Spieler
        self.spielerV = []  # Rechtecke für virtuelle Spieler
        self.rec_height = 0.05
        self.rec_width = 0.04
        for i in self.I:
            # Rechtecke repräsentieren Spieler
            self.spieler.append(Rectangle(self.posit['s{}'.format(self.ST[i][0])], self.rec_width,
                                          self.rec_height, linewidth=1, edgecolor='black',
                                          facecolor=self.colors[self.ST[i][1]-1 % len(self.colors)]))
            # hinzufügen Spieler zu subplot
            self.ax.add_patch(self.spieler[-1])
            # hinzufügen Nummern (für Spieler) zu subplot
            self.numbers.append(self.ax.text(self.posit['s{}'.format(self.ST[i][0])][0],
                                             self.posit['s{}'.format(self.ST[i][0])][1], str(i+1), fontsize=6))
        for e in self.E:
            self.spielerV.append(deque())

        # erste Liste in "first" enthält alle vorkommenden Positionen, zweite Liste den Spieler mit kleinstem Index in
        # dieser Position
        first = [[], []]
        for i in self.I:  # Plot vor Start des Programms
            currentPos = self.ST[i][0]  # Spieler befinden sich zum Zeitpunkt '0' in ihren Quellen
            if currentPos not in first[0]:  # erster Spieler in Position "currentPos"
                first[0].append(currentPos)
                first[1].append(i)
                self.spieler[i].set_xy(self.posit[currentPos]) # neue Position Spieler
                self.numbers[i].set_x(self.posit[currentPos][0]+self.rec_width/3)  # neue Position Nummer
                self.numbers[i].set_y(self.posit[currentPos][1]+self.rec_height/4)
            else:  # alle weiteren Spieler in Position "currentPos"
                count = self.ST[:i][0].count(currentPos)
                x, y = self.spieler[first[1][first[0].index(currentPos)]].get_xy()
                # setze Rechtecke der Spieler in gleicher Position übereinander
                self.spieler[i].set_xy((x, y+count*self.rec_height))
                self.numbers[i].set_x(x + self.rec_width/3)
                self.numbers[i].set_y(y+count*self.rec_height + self.rec_height/4)

        nx.draw(self.graph, self.posit, node_color=self.colorMap, with_labels=True, alpha=0.8)  # zeichnen des Graphen
        plt.title('Theta = 0')  # setze Titel des Plots
        self.fig.canvas.draw()

        #toolbar = NavigationToolbar2Tk(canvas, self.master)
        #toolbar.update()
        #canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        #canvas.get_tk_widget().grid(row=1, column=1)

        self.unpaused = False
        self.zeitpunkt = 0

        # Rahmen für die Buttons
        button_frame = tk.Frame(self)
        #button_frame.grid(row=0, column=0, sticky="n")
        #button_frame.pack(fill=None, expand=True, side=tk.TOP, padx=int(self.screenwidth/8))
        #button_frame.pack(fill=None, padx=int(self.screenwidth/8))
        # button_frame.config(background="#0000ee")
        # self.master.update_idletasks()
        #self.framewidth = button_frame.winfo_width()
        #self.frameheight = button_frame.winfo_height()
        # button_frame.update()
        # self.master.update_idletasks()

        # Button: Zurück
        self.prev = tk.Button(canvas.get_tk_widget(), text="Zurück", state="disabled")
        #self.prev["bg"] = "#14ADA0"  # Farbe
        #self.prev["fg"] = "#FFFFFF"  # Schriftfarbe: weiß
        self.prev["relief"] = "flat"
        self.prev["height"] = 2
        self.prev["width"] = 19
        self.prev["font"] = "Arial 12 bold"
        #self.prev.pack(padx = int(1/3 * self.screenwidth), pady = int(11/12 * self.screenheight))
        self.prev.pack(side=tk.BOTTOM, padx=0, pady=0)
        #self.prev.grid(row=2,column=1,sticky="we")
        #button1_window = canvas.create_window(10, 10, anchor=NW, window=self.prev)

        # Button: Pause
        self.pause = tk.Button(canvas.get_tk_widget(), text="Pause", state="disabled")
        #self.pause["bg"] = "#1894CE"
        #self.pause["fg"] = "#FFFFFF"
        self.pause["relief"] = "flat"
        self.pause["height"] = 2
        self.pause["width"] = 19
        self.pause["font"] = "Arial 12 bold"
        self.pause.pack(side=tk.LEFT, padx=0, pady=0)

        # Button: Weiter
        self.nex = tk.Button(canvas.get_tk_widget(), text="Weiter", state="disabled")
        #self.nex["bg"] = "#14ADA0"  # Farbe
        #self.nex["fg"] = "#FFFFFF"  # Schriftfarbe: weiß
        self.nex["relief"] = "flat"  # entfernt Rand
        self.nex["height"] = 2
        self.nex["width"] = 19
        self.nex["font"] = "Arial 12 bold"
        self.nex.pack(side=tk.RIGHT, padx=0, pady=0)

    def redraw(self, theta):
        """
        :param theta: Zeitpunkt der Zeichnung
        :return: kein Rückgabewert
        """
        plt.title('Theta = {}'.format(theta))  # setze Titel des Plots
        self.fig.canvas.draw()
        return

    def remove_artist(self, ind):
        """
        entfernt virtuellen Spieler aus Plot
        :param ind: Index des zu entfernenden virtuellen Spielers
        :return: kein Rückgabewert
        """
        artist = self.spielerV[ind].popleft()
        artist.remove()
        return

    def restore_players(self, positionenV, red_list):
        """
        Wird aufgerufen von "Application.back()". Spieler, die bei einem Rückwärtsschritt zurück in eine Warteschlange
        kommen, müssen auch in Plot berücksichtigt werden: Viruelle Spieler werden an passende Position gesetzt, reale
        Spieler werden rot gefärbt.
        :param positionenV: Liste welche die Reihenfolge der virtuellen Spieler enthält, in der diese gezeichnet werden
        :param red_list: Liste der realen Spieler, die rot gefärbt werden
        :return: kein Rückgabewert
        """
        for e in self.E:
            for pos in positionenV:  # zeichne virtuelle Spieler
                self.spielerV[self.E.index(e)].append(Rectangle((0.91 * self.posit[e[0]][0] +
                                                                 0.09 * self.posit[e[1]][0],
                                                                 0.91 * self.posit[e[0]][1] +
                                                                 0.09 * self.posit[e[1]][1] +
                                                                 self.rec_height * pos),
                                                                 self.rec_width,
                                                                 self.rec_height,
                                                                 linewidth=1,
                                                                 edgecolor='black',
                                                                 facecolor='white', alpha=0.5))
                self.ax.add_patch(self.spielerV[self.E.index(e)][-1])
        for out in red_list:  # Färbe Spieler in Warteschlange rot
            self.spieler[self.spieler.index(out)].set_facecolor('red')
        return

    def add_virtual_player(self, kante, deq_pos):
        """
        speichert virtuellen Spieler für Kante 'kante' und Warteschlangenposition 'deq_pos' in 'self.SpielerV' und fügt
        diesen zu Plot hinzu
        :param kante: Kante
        :param deq_pos: Warteschlangenposition
        :return: kein Rückgabewert
        """
        # speichere virtuelle Spieler in dict., um diese aus Plot entfernen zu können
        self.spielerV[self.E.index(kante)].appendleft(
            Rectangle((0.91*self.posit[kante[0]][0] + 0.09*self.posit[kante[1]][0],
                       0.91 * self.posit[kante[0]][1] + 0.09*self.posit[kante[1]][1] + self.rec_height * deq_pos),
                        self.rec_width, self.rec_height, linewidth=1, edgecolor='black', facecolor= 'white', alpha=0.5))
        # füge virtuellen Spieler zu Plot hinzu
        self.ax.add_patch(self.spielerV[self.E.index(kante)][0])
        return

    def change_color(self, player_color, t_index):
        """
        Ändert die Farbe der übergebenen Spieler im Plot
        :param player_color: Liste der Spieler
        :param t_index: Liste mit gleicher Länge wie 'player_color'. Enthält für jeden Spieler den Index des
        entsprechenden Zielknotens (und bestimmt somit die Farbe)
        :return: kein Rückgabewert
        """
        for p in range(len(player_color)):
            self.spieler[self.spieler.index(player_color[p])].set_facecolor(self.colors[t_index[p] % len(self.colors)])
        return

    def updates(self, visibility, nex_config, prev_config):
        """
        Setzt Spieler in Liste 'visibility' auf sichtbar. Aktiviert (deaktiviert) gegebenenfalls den Button 'self.nex'
        ('self.prev')
        :param visibility: Liste der Spielerindizes die sichtbar werden sollen
        :param nex_config: Boolean. 'True', falls Button 'self.nex' aktiviert werden soll
        :param prev_config: Boolean. 'True', falls Button 'self.prev' deaktiviert werden soll
        :return: kein Rückgabewert
        """
        for i in visibility:
            self.spieler[i].set_visible(True)  # aktualisiere Sichtbarkeit
            self.numbers[i].set_visible(True)

        if nex_config:
            self.nex.config(state="normal")

        if prev_config:
            self.prev.config(state="disabled")
        return

    def set_nex_status(self, status):
        """
        set-Funktion zur Aktivierung/Deaktivierung des Buttons 'self.nex'
        :param status: neuer Status, 2 mögliche Werte: "normal", "disabled"
        :return: kein Rückgabewert
        """
        self.nex.config(state="{}".format(status))
        return

    def set_prev_status(self, status):
        """
        set-Funktion zur Aktivierung/Deaktivierung des Buttons 'self.prev'
        :param status: neuer Status, 2 mögliche Werte: "normal", "disabled"
        :return: kein Rückgabewert
        """
        self.prev.config(state="{}".format(status))
        return

    def set_pause_text(self, tex):
        """
        set-Funktion zur Bearbeitung des Texts des Buttons 'self.pause'
        :param tex: neuer Text, 2 mögliche Werte: "Pause", "Fortsetzen"
        :return:
        """
        self.pause.config(text="{}".format(tex))
        return

    def pop_right(self, edge):
        """
        entfernt nächsten virtuellen Spieler aus der Warteschlange der übergebenen Kante
        :param edge: Index der betrachteten Kante
        :return: kein Rückgabewert
        """
        if len(self.spielerV[edge]) > 0:
            v = self.spielerV[edge].pop()  # entferne virtuellen Spieler
            v.remove()  # entferne virtuellen Spieler aus plot
        return

#nicht mehr benötigt?
    """def set_player_position(self, i, position):
        
        aufgerufen von "Application.py", setzt Position des übergebenen Spielers für Plot
        :param i: Index des Spielers
        :param position: aktuelle Position des Spielers (bspw. in Knoten, Punkt im Durchlauf einer Kante, etc.)
        :return: kein Rückgabewert
        
        # neue Position Spieler
        self.spieler[i].set_xy(self.posit[position])
        # neue Position Nummer
        self.numbers[i].set_x(self.posit[position][0]+self.rec_width/3)
        self.numbers[i].set_y(self.posit[position][1]+self.rec_height/4)
        return"""

    def draw_new_positions(self, positions):
        # erste Liste in "first" enthält alle vorkommenden Positionen, zweite Liste den Spieler mit kleinstem Index in
        # dieser Position
        first = [[], []]
        for i in self.I:
            if positions[i] == "t{}".format(self.ST[i][1]):  # Spieler werden nicht mehr angezeigt, wenn Senke erreicht
                self.spieler[i].set_visible(False)
                self.numbers[i].set_visible(False)
            else:
                if positions[i] in self.V:
                    if positions[i] not in first[0]:
                        first[0].append(positions[i])
                        first[1].append(i)
                        # neue Position Spieler
                        self.spieler[i].set_xy(self.posit[positions[i]])
                        # neue Position Nummer
                        self.numbers[i].set_x(self.posit[positions[i]][0]+self.rec_width/3)
                        self.numbers[i].set_y(self.posit[positions[i]][1]+self.rec_height/4)
                    else:
                        # zähle Spieler in gleicher Position mit kleinerem Index
                        count = positions[:i].count(positions[i])
                        # ziehe Spieler die an ihrem Zielknoten angekommen (und nicht sichtbar) sind, ab
                        for sp in self.I[:i]:
                            if positions[sp] == positions[i] and positions[sp] == "t{}".format(self.ST[sp][1]):
                                count -= 1
                        x, y = self.spieler[first[1][first[0].index(positions[i])]].get_xy()
                        # setze Rechtecke der Spieler in gleicher Position übereinander
                        self.spieler[i].set_xy((x, y+count*self.rec_height))
                        self.numbers[i].set_x(x + self.rec_width/3)
                        self.numbers[i].set_y(y+count*self.rec_height + self.rec_height/4)
                elif positions[i][1] != 0:
                    if positions[i] not in first[0]:
                        first[0].append(positions[i])
                        first[1].append(i)
                        self.spieler[i].set_xy((1-positions[i][1]/self.r[self.E.index(positions[i][0])]) *
                                               self.posit[positions[i][0][0]] +
                                               (positions[i][1]/self.r[self.E.index(positions[i][0])]) *
                                               self.posit[positions[i][0][1]])  # neue Position Spieler auf Kante
                        self.numbers[i].set_x(((1-positions[i][1]/self.r[self.E.index(positions[i][0])]) *
                                               self.posit[positions[i][0][0]] +
                                               (positions[i][1]/self.r[self.E.index(positions[i][0])]) *
                                               self.posit[positions[i][0][1]])[0]+self.rec_width/3)
                        self.numbers[i].set_y(((1-positions[i][1]/self.r[self.E.index(positions[i][0])]) *
                                               self.posit[positions[i][0][0]] +
                                               (positions[i][1]/self.r[self.E.index(positions[i][0])]) *
                                               self.posit[positions[i][0][1]])[1]+self.rec_height/4)
                    else:
                        count = positions[:i].count(positions[i])
                        # Koordinaten des ersten Spielers
                        x, y = self.spieler[first[1][first[0].index(positions[i])]].get_xy()
                        # setze Rechtecke der Spieler in gleicher Position übereinander
                        self.spieler[i].set_xy((x,y+count*self.rec_height))
                        self.numbers[i].set_x(x+self.rec_width/3)
                        self.numbers[i].set_y(y+count*self.rec_height+self.rec_height/4)
        return

    def draw_new_queue_positions(self, deques):
        for e in self.E:  # setzen der Positionen von Spielern in Warteschlange
            j = len(deques[self.E.index(e)]) -1  # erster Index in der deque
            vcount = 0
            while j >= 0:
                try:
                    # TypeError, falls Spieler "self.deques[self.E.index(e)][j]" KEIN virtueller Spieler
                    int(deques[self.E.index(e)][j])
                    # setze Position virtueller Spieler
                    self.spielerV[self.E.index(e)][vcount].set_xy((0.91*self.posit[e[0]][0] + 0.09*self.posit[e[1]][0],
                                                                   0.91*self.posit[e[0]][1] + 0.09*self.posit[e[1]][1] +
                                                                   self.rec_height *
                                                                   (len(deques[self.E.index(e)]) -1 -j)))
                    vcount += 1
                except TypeError:
                    deques[self.E.index(e)][j].set_facecolor('red')
                    # setze Position von Spieler "self.deques[self.E.index(e)][j]"
                    self.spieler[self.spieler.index(deques[self.E.index(e)][j])].set_xy((0.91*self.posit[e[0]][0] +
                                                                                              0.09*self.posit[e[1]][0],
                                                                                              0.91*self.posit[e[0]][1] +
                                                                                              0.09*self.posit[e[1]][1] +
                                                                self.rec_height*(len(deques[self.E.index(e)]) -1 -j)))
                    self.numbers[self.spieler.index(deques[self.E.index(e)][j])].set_x(0.91*self.posit[e[0]][0] +
                                                                                            0.09*self.posit[e[1]][0] +
                                                                                            self.rec_width/3)
                    self.numbers[self.spieler.index(deques[self.E.index(e)][j])].set_y(0.91*self.posit[e[0]][1] +
                                                                                            0.09*self.posit[e[1]][1] +
                                                                                            self.rec_height/4 +
                                                                                            self.rec_height *
                                                                                (len(deques[self.E.index(e)]) -1 -j))
                j -= 1
        return
