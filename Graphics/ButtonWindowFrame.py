import tkinter as tk
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import networkx as nx
from collections import deque
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import sys
#import time


class ButtonWindowFrame(tk.Frame):
    """
    Diese Klasse verwaltet die Button-Leiste am unteren Bildschirmrand, mit den Buttons:
    "Zurück", "Pause", bzw. "Fortsetzen", und "Weiter". Objekte werden erzeugt in der Initialisierung von
    "Application.py".
    """

    def __init__(self, E, V, ST, r, posit):
        """
        Erzeugt Buttons und ersten Plot mit Spielern in deren Startknoten.
        :param E: Liste der Kanten, benötigt u.a. für Positionen der Spieler
        :param V: Liste der Knoten des Graphs
        :param ST: Liste von Tupeln, wobei k-tes Tupel (i,j) beschreibt, dass Spieler k +1 Quelle s_i und Senke t_j
         besitzt
        :param r: Liste der Reisezeiten aller Kanten
        :param posit: Knotenpositionen. Werden bestimmt/berechnet in "Main.py" und von dort an "Application.py"
         übergeben
        """
        self.E = E.copy()
        self.V = V.copy()
        self.ST = ST.copy()
        self.I = range(len(ST))
        self.r = r.copy()
        self.master = tk.Tk(className="graph")
        self.fig = plt.figure(figsize=(3,3), dpi=100)
        self.ax = self.fig.add_subplot(111)
        plt.ylim(-1.5, 1.5)  # setzte Höhe des Plots, damit Rechtecke in jedem Plot gleich groß
        plt.xlim(-1.25, 1.25)
        self.screenwidth = self.master.winfo_screenwidth()
        self.screenheight = self.master.winfo_screenheight()
        self.master.geometry("{}x{}".format(int(self.screenwidth), int(self.screenheight) -100))
        super().__init__(self.master)
        self.master.protocol('WM_DELETE_WINDOW', sys.exit)  # beendet Programm bei Klicken des 'X'-Buttons
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.master)

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
                    farbe = int(name[1:])  # ValueError, falls "v" kein Start- oder Zielknoten ist
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
            self.graph.add_edge(e[0], e[1])

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
        # erste Liste in "last_pos" enthält für jede vorkommende Position, falls vorhanden, die Position des letzten
        # Spielers im Plot (alle Spieler darüber passen nicht mehr hinein und werden nicht angezeigt). Zweite Liste
        # enthält den Spielerindex des letzten Spielers in dieser Position (!= letzter im Plot), dieser soll im Plot
        # angezeigt werden, und andeuten, dass alle Spieler bis zu diesem letzten Spieler in der Warteschlange sind,
        # aber aus Platzgründen nicht angezeigt werden.
        last_pos = [[], []]
        for i in self.I:  # Plot vor Start des Programms
            currentPos = 's{}'.format(self.ST[i][0])  # Spieler befinden sich zum Zeitpunkt '0' in ihren Quellen
            if currentPos not in first[0]:  # erster Spieler in Position "currentPos"
                first[0].append(currentPos)
                first[1].append(i)
                last_pos[0].append(None)
                last_pos[1].append(None)
                self.spieler[i].set_xy(self.posit[currentPos])  # neue Position Spieler
                self.numbers[i].set_x(self.posit[currentPos][0]+self.rec_width/4)  # neue Position Nummer
                self.numbers[i].set_y(self.posit[currentPos][1]+self.rec_height/4)
            else:  # alle weiteren Spieler in Position "currentPos"
                # zähle Spieler die sich in gleicher Position befinden wie Spieler "i" und kleineren Index haben
                count = [s[0] for s in self.ST[:i]].count(int(currentPos[1:]))
                # Koordinaten des ersten Spielers in gleicher Position wie Spieler "i"
                x, y = self.spieler[first[1][first[0].index(currentPos)]].get_xy()
                if y + count*self.rec_height < 1.45:  # prüfe, ob Spieler in den Plot passt
                    # setze Rechtecke der Spieler in gleicher Position übereinander
                    self.spieler[i].set_xy((x, y+count*self.rec_height))
                    self.numbers[i].set_x(x + self.rec_width/4)
                    self.numbers[i].set_y(y+count*self.rec_height + self.rec_height/4)
                else:  # falls nicht, so wird er auch nicht gezeichnet
                    if last_pos[0][first[0].index(currentPos)] is None:
                        last_pos[0][first[0].index(currentPos)] = (x, y + count*self.rec_height)
                        last_pos[1][first[0].index(currentPos)] = i
                    else:
                        last_pos[1][first[0].index(currentPos)] = i
                    self.spieler[i].set_visible(False)
                    self.numbers[i].set_visible(False)

        for last in last_pos[1]:
            if last is not None:
                x, y = last_pos[0][last_pos[1].index(last)]
                self.spieler[last].set_xy((x, y))
                self.numbers[last].set_x(x + self.rec_width/4)
                self.numbers[last].set_y(y + self.rec_height/4)
                self.spieler[last].set_visible(True)
                self.numbers[last].set_visible(True)

        self.zeitpunkt = 0

        # Rahmen für die Buttons
        button_frame = tk.Frame(self.master)

        # Button: Zurück
        self.prev = tk.Button(button_frame, text="Zurück", state="disabled")
        #self.prev["bg"] = "#14ADA0"  # Farbe
        #self.prev["fg"] = "#FFFFFF"  # Schriftfarbe: weiß
        self.prev["relief"] = "flat"
        self.prev["height"] = 2
        self.prev["width"] = 19
        self.prev["font"] = "Arial 12 bold"
        self.prev.pack(side=tk.LEFT, padx=0, pady=0)

        # Button: Pause
        self.pause = tk.Button(button_frame, text="Pause", state="disabled")
        #self.pause["bg"] = "#1894CE"
        #self.pause["fg"] = "#FFFFFF"
        self.pause["relief"] = "flat"
        #self.pause["height"] = 2
        #self.pause["width"] = 19
        self.pause["font"] = "Arial 12 bold"
        self.pause.pack(side=tk.LEFT, padx=0, pady=0)

        # Button: Weiter
        self.nex = tk.Button(button_frame, text="Weiter", state="disabled")
        #self.nex["bg"] = "#14ADA0"  # Farbe
        #self.nex["fg"] = "#FFFFFF"  # Schriftfarbe: weiß
        self.nex["relief"] = "flat"  # entfernt Rand
        self.nex["height"] = 2
        self.nex["width"] = 19
        self.nex["font"] = "Arial 12 bold"
        self.nex.pack(side=tk.LEFT, padx=0, pady=0)

        button_frame.pack(side=tk.BOTTOM)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        nx.draw(self.graph, self.posit, node_color=self.colorMap, with_labels=True, alpha=0.8)  # zeichnen des Graphen
        plt.title('Theta = 0')  # setze Titel des Plots

    def redraw(self, theta):
        """
        Wird aufgerufen, um Plot zu aktualisieren
        :param theta: Zeitpunkt des Plots
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
        kommen, müssen auch in Plot berücksichtigt werden: Virtuelle Spieler werden an passende Position gesetzt, reale
        Spieler werden rot gefärbt.
        :param positionenV: Liste welche die Reihenfolge der virtuellen Spieler enthält, in der diese gezeichnet werden
        :param red_list: Liste der realen Spieler, die rot gefärbt werden
        :return: kein Rückgabewert
        """
        for e in self.E:
            for pos in positionenV[self.E.index(e)]:  # zeichne virtuelle Spieler
                self.spielerV[self.E.index(e)].append(Rectangle((0.91 * self.posit[e[0]][0] +
                                                                 0.09 * self.posit[e[1]][0],
                                                                 0.91 * self.posit[e[0]][1] +
                                                                 0.09 * self.posit[e[1]][1] + self.rec_height * pos),
                                                                 self.rec_width, self.rec_height, linewidth=1,
                                                                 edgecolor='black', facecolor='white', alpha=0.5))
                self.ax.add_patch(self.spielerV[self.E.index(e)][-1])
        for out in red_list:  # Färbe Spieler in Warteschlange rot
            self.spieler[out[0]].set_facecolor('red')
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
            player_index = player_color[p][0]
            self.spieler[player_index].set_facecolor(self.colors[t_index[p] % len(self.colors)])
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

    def draw_new_positions(self, positions):
        """
        Dient dazu, Positionen aller realer Spieler, welche sich zum aktuellen Zeitpunkt nicht in einer Wartschlange
        befinden, korrekt für den Plot zu setzen. (Für Spieler in Warteschlange siehe self.draw_new_queue_positions )
        :param positions: Liste der Positionen aller Spieler, diese Funktion berücksichtigt nur 2 Fälle von möglichen
        Positionen: 1. Spieler befindet sich in einem Knoten -> Position 'v'
                    2. Spieler befindet sich auf einer Kante, aber nicht in deren Wartschlange -> Position '(e,k)',
                       wobei 'e' Kante und 'k' > 0
        :return: kein Rückgabewert
        """
        # erste Liste in "first" enthält alle vorkommenden Positionen, zweite Liste den Spieler mit kleinstem Index in
        # dieser Position, dritte Liste die momentane Anzahl Spieler in dieser Position
        first = [[], [], []]
        # erste Liste in "last_pos" enthält für jede vorkommende Position, falls vorhanden, die Position des letzten
        # Spielers im Plot (alle Spieler darüber passen nicht mehr hinein und werden nicht angezeigt). Zweite Liste
        # enthält den Spielerindex des letzten Spielers in dieser Position (!= letzter im Plot), dieser soll im Plot
        # angezeigt werden, und andeuten, dass alle Spieler bis zu diesem letzten Spieler in der Warteschlange sind,
        # aber aus Platzgründen nicht angezeigt werden.
        last_pos = [[], []]
        for i in self.I:
            if positions[i] == "t{}".format(self.ST[i][1]):  # Spieler werden nicht mehr angezeigt, wenn Senke erreicht
                self.spieler[i].set_visible(False)
                self.numbers[i].set_visible(False)
            else:
                if positions[i] in self.V:
                    if positions[i] not in first[0]:
                        first[0].append(positions[i])
                        first[1].append(i)
                        first[2].append(0)
                        last_pos[0].append(None)
                        last_pos[1].append(None)
                        # neue Position Spieler
                        self.spieler[i].set_xy(self.posit[positions[i]])
                        # neue Position Nummer
                        self.numbers[i].set_x(self.posit[positions[i]][0]+self.rec_width/4)
                        self.numbers[i].set_y(self.posit[positions[i]][1]+self.rec_height/4)
                        self.spieler[i].set_visible(True)
                        self.numbers[i].set_visible(True)
                    else:
                        pos = first[0].index(positions[i])
                        first[2][pos] += 1  # zählt Spieler welche sich momentan in Position "positions[i]" befinden
                        count = first[2][pos]
                        x, y = self.spieler[first[1][first[0].index(positions[i])]].get_xy()
                        if y+count*self.rec_height < 1.45:  # prüfe, ob Spieler noch in den Plot passt
                            # setze Rechtecke der Spieler in gleicher Position übereinander
                            self.spieler[i].set_xy((x, y+count*self.rec_height))
                            self.numbers[i].set_x(x + self.rec_width/4)
                            self.numbers[i].set_y(y+count*self.rec_height + self.rec_height/4)
                            self.spieler[i].set_visible(True)
                            self.numbers[i].set_visible(True)
                        else:  # Spieler passt nicht in den Plot
                            # merke letzten Spieler in dieser Position, um anzudeuten, wie viele Spieler sich noch in
                            # dieser Position befinden
                            if last_pos[0][first[0].index(positions[i])] is None:
                                last_pos[0][first[0].index(positions[i])] = (x, y + count*self.rec_height)
                                last_pos[1][first[0].index(positions[i])] = i
                            else:
                                last_pos[1][first[0].index(positions[i])] = i
                            self.spieler[i].set_visible(False)
                            self.numbers[i].set_visible(False)
                elif positions[i][1] != 0:
                    if positions[i] not in first[0]:
                        first[0].append(positions[i])
                        first[1].append(i)
                        first[2].append(0)
                        last_pos[0].append(None)
                        last_pos[1].append(None)
                        self.spieler[i].set_xy((1-positions[i][1]/self.r[self.E.index(positions[i][0])]) *
                                               self.posit[positions[i][0][0]] +
                                               (positions[i][1]/self.r[self.E.index(positions[i][0])]) *
                                               self.posit[positions[i][0][1]])  # neue Position Spieler auf Kante
                        self.numbers[i].set_x(((1-positions[i][1]/self.r[self.E.index(positions[i][0])]) *
                                               self.posit[positions[i][0][0]] +
                                               (positions[i][1]/self.r[self.E.index(positions[i][0])]) *
                                               self.posit[positions[i][0][1]])[0]+self.rec_width/4)
                        self.numbers[i].set_y(((1-positions[i][1]/self.r[self.E.index(positions[i][0])]) *
                                               self.posit[positions[i][0][0]] +
                                               (positions[i][1]/self.r[self.E.index(positions[i][0])]) *
                                               self.posit[positions[i][0][1]])[1]+self.rec_height/4)
                        self.spieler[i].set_visible(True)
                        self.numbers[i].set_visible(True)
                    else:
                        pos = first[0].index(positions[i])
                        first[2][pos] += 1  # zählt Spieler welche sich momentan in Position "positions[i]" befinden
                        count = first[2][pos]
                        # Koordinaten des ersten Spielers
                        x, y = self.spieler[first[1][first[0].index(positions[i])]].get_xy()
                        if y+count*self.rec_height < 1.45:  # prüfe, ob Spieler noch in den Plot passt
                            # setze Rechtecke der Spieler in gleicher Position übereinander
                            self.spieler[i].set_xy((x,y+count*self.rec_height))
                            self.numbers[i].set_x(x+self.rec_width/4)
                            self.numbers[i].set_y(y+count*self.rec_height+self.rec_height/4)
                            self.spieler[i].set_visible(True)
                            self.numbers[i].set_visible(True)
                        else:  # Spieler passt nicht in den Plot
                            # merke letzten Spieler in dieser Position, um anzudeuten, wie viele Spieler sich noch in
                            # dieser Position befinden
                            if last_pos[0][first[0].index(positions[i])] is None:
                                last_pos[0][first[0].index(positions[i])] = (x, y + count*self.rec_height)
                                last_pos[1][first[0].index(positions[i])] = i
                            else:
                                last_pos[1][first[0].index(positions[i])] = i
                            self.spieler[i].set_visible(False)
                            self.numbers[i].set_visible(False)

        # zeichne alle Spieler in "last_pos[1]"
        for last in last_pos[1]:
            if last is not None:
                x, y = last_pos[0][last_pos[1].index(last)]
                self.spieler[last].set_xy((x, y))
                self.numbers[last].set_x(x + self.rec_width/4)
                self.numbers[last].set_y(y + self.rec_height/4)
                self.spieler[last].set_visible(True)
                self.numbers[last].set_visible(True)

        return

    def draw_new_queue_positions(self, deques):
        """
        Dient dazu, Positionen aller Spieler (real + virtuell), die sich momentan in einer Wartschlange befinden, für
        den Plot zu setzen. Für alle anderen Spieler, siehe "self.draw_new_positions".
        :param deques: Liste aller aktuellen Warteschlangen
        :return: kein Rückgabewert
        """
        for e in self.E:  # setzen der Positionen von Spielern in Warteschlange
            # in "last_pos" wird für die aktuelle Warteschlange die letzte Position gespeichert, in welcher die Spieler
            # aus dem Plot noch angezeigt werden (danach wird abgeschnitten, da Spieler nicht mehr in den Plot passen).
            # Außerdem wird der letzte Spielerindex in der aktuellen Wartschlange gespeichert, damit der entsprechende
            # Spieler als letztes im Plot dargestellt werden kann.
            last_pos = None
            j = len(deques[self.E.index(e)]) -1  # erster Index in der deque
            vcount = 0
            while j >= 0:
                try:
                    pos = 0.91*self.posit[e[0]][1] + 0.09*self.posit[e[1]][1]
                    # TypeError, falls Spieler "self.deques[self.E.index(e)][j]" KEIN virtueller Spieler
                    int(deques[self.E.index(e)][j])
                    # prüfe, ob Spieler noch in den Plot passt
                    if pos + self.rec_height * (len(deques[self.E.index(e)]) -1 -j) < 1.45:
                        # setze Position virtueller Spieler
                        self.spielerV[self.E.index(e)][vcount].set_xy((0.91*self.posit[e[0]][0] +
                                                                       0.09*self.posit[e[1]][0],
                                                                       0.91*self.posit[e[0]][1] +
                                                                       0.09*self.posit[e[1]][1] +
                                                                       self.rec_height *
                                                                       (len(deques[self.E.index(e)]) -1 -j)))
                        self.spielerV[self.E.index(e)][vcount].set_visible(True)
                    else:  # virtueller Spieler passt nicht mehr in den Plot
                        self.spielerV[self.E.index(e)][vcount].set_visible(False)
                        #self.spielerV[self.E.index(e)][vcount].set_facecolor('blue')
                    vcount += 1
                except TypeError:
                    # prüfe, ob Spieler noch in den Plot passt
                    if pos + self.rec_height * (len(deques[self.E.index(e)]) -1 -j) < 1.45:
                        self.spieler[deques[self.E.index(e)][j][0]].set_facecolor('red')
                        # setze Position von Spieler "self.deques[self.E.index(e)][j]"
                        self.spieler[deques[self.E.index(e)][j][0]].set_xy((0.91*self.posit[e[0]][0] +
                                                                                  0.09*self.posit[e[1]][0],
                                                                                  0.91*self.posit[e[0]][1] +
                                                                                  0.09*self.posit[e[1]][1] +
                                                                self.rec_height*(len(deques[self.E.index(e)]) -1 -j)))
                        self.numbers[deques[self.E.index(e)][j][0]].set_x(0.91*self.posit[e[0]][0] +
                                                                          0.09*self.posit[e[1]][0] +
                                                                          self.rec_width/4)
                        self.numbers[deques[self.E.index(e)][j][0]].set_y(0.91*self.posit[e[0]][1] +
                                                                                0.09*self.posit[e[1]][1] +
                                                                                self.rec_height/4 +
                                                                                self.rec_height *
                                                                                (len(deques[self.E.index(e)]) -1 -j))
                        self.spieler[deques[self.E.index(e)][j][0]].set_visible(True)
                        self.numbers[deques[self.E.index(e)][j][0]].set_visible(True)
                    else:  # Spieler passt nicht mehr in den Plot
                        # merke letzten Spieler in dieser Position, um anzudeuten, wie viele Spieler sich noch in
                        # dieser Warteschlange befinden
                        if last_pos is None:
                            last_pos = [deques[self.E.index(e)][j][0], j]
                        else:
                            last_pos[0] = deques[self.E.index(e)][j][0]
                        self.spieler[deques[self.E.index(e)][j][0]].set_visible(False)
                        self.numbers[deques[self.E.index(e)][j][0]].set_visible(False)
                j -= 1

            # zeichne Spieler in "last_pos"
            if last_pos is not None:
                self.spieler[last_pos[0]].set_facecolor('red')
                self.spieler[last_pos[0]].set_xy((0.91*self.posit[e[0]][0] + 0.09*self.posit[e[1]][0],
                                                  0.91*self.posit[e[0]][1] +0.09*self.posit[e[1]][1] +
                                                  self.rec_height*(len(deques[self.E.index(e)]) -1 -last_pos[1])))
                self.numbers[last_pos[0]].set_x(0.91*self.posit[e[0]][0] + 0.09*self.posit[e[1]][0] + self.rec_width/4)
                self.numbers[last_pos[0]].set_y(0.91*self.posit[e[0]][1] + 0.09*self.posit[e[1]][1] +
                                                                  self.rec_height/4 +
                                                                  self.rec_height *
                                                                  (len(deques[self.E.index(e)]) -1 -last_pos[1]))
                self.spieler[last_pos[0]].set_visible(True)
                self.numbers[last_pos[0]].set_visible(True)

        return
