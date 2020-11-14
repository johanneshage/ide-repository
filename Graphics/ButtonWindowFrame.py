import tkinter as tk
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.widgets import Button


class ButtonWindowFrame(tk.Frame):
    """
    Diese Klasse verwaltet die Button-Leiste am unteren Bidlschirmrand, mit den Buttons:
    "Zurück", "Pause", bzw. "Fortsetzen", und "Weiter"
    """

    def __init__(self, app):
        """
        Erzeugt Buttons und Layout der Button-Leiste
        :param app: Objekt der Klasse "Application", wird übergeben, damit durch Klicken der Buttons auf die Haupt-
         anwendung zugegriffen werden kann
        """
        self.master = tk.Tk()
        self.app = app
        self.screenwidth = self.master.winfo_screenwidth()
        self.screenheight = self.master.winfo_screenheight()
        self.master.wm_attributes('-topmost', True)
        self.master.overrideredirect(True)
        self.master.geometry("{}x{}+0+{}".format(int(self.screenwidth), int(self.screenheight/6),
                                                 int(self.screenheight * 7/8)))
        self.master.config(background='black')
        #self.master.update_idletasks()
        #self.master.config(width = self.master.winfo_screenwidth(), height = self.master.winfo_screenheight()/6) #,
        # background="#0000ee"
        super().__init__(self.master)
        #self.pack()
        self.pack(fill="both", expand=True)

        self.unpaused = False
        self.zeitpunkt = 0

        # Rahmen für die Buttons
        button_frame = tk.Frame(self)
        button_frame.pack(fill=None, expand=True, side=tk.TOP, padx=int(self.screenwidth/8))
        #button_frame.pack(fill=None, padx=int(self.screenwidth/8))
        # button_frame.config(background="#0000ee")
        # self.master.update_idletasks()
        self.framewidth = button_frame.winfo_width()
        self.frameheight = button_frame.winfo_height()
        # button_frame.update()
        # self.master.update_idletasks()

        self.prev = tk.Button(button_frame, text="Zurück", command=self.back, state="disabled")  # Button: Zurück
        self.prev["bg"] = "#14ADA0"  # Farbe
        self.prev["fg"] = "#FFFFFF"  # Schriftfarbe: weiß
        self.prev["relief"] = "flat"
        self.prev["height"] = 2
        self.prev["width"] = 19
        self.prev["font"] = "Arial 12 bold"
        # self.prev.pack(padx = int(1/3 * self.screenwidth), pady = int(11/12 * self.screenheight))
        self.prev.pack(side=tk.LEFT, padx=0, pady=0)

        self.pause = tk.Button(button_frame, text="Pause", command=self.pause, state="disabled")  # Button: Pause
        self.pause["bg"] = "#1894CE"
        self.pause["fg"] = "#FFFFFF"
        self.pause["relief"] = "flat"
        self.pause["height"] = 2
        self.pause["width"] = 19
        self.pause["font"] = "Arial 12 bold"
        self.pause.pack(side=tk.LEFT, padx=0, pady=0)

        self.nex = tk.Button(button_frame, text ="Weiter", command = self.weiter, state="disabled")  # Button: Weiter
        self.nex["bg"] = "#14ADA0"  # Farbe
        self.nex["fg"] = "#FFFFFF"  # Schriftfarbe: weiß
        self.nex["relief"] = "flat"  # entfernt Rand
        self.nex["height"] = 2
        self.nex["width"] = 19
        self.nex["font"] = "Arial 12 bold"
        self.nex.pack(side=tk.LEFT, padx=0, pady=0)

    def pause(self):
        """
        Bei Aufruf wird Pause-Status geändert (Fortsetzen falls pausiert und andersrum)
        :return: Kein Rückgabewert
        """
        if self.get_unpaused():
            self.set_unpaused(False)
            self.pause.config(text="Fortsetzen")
        else:
            self.set_unpaused(True)
            self.pause.config(text="Pause")
        return

    def back(self):
        """
        Bei Aufruf wird in "app" ein Zeitschritt zurück gegangen. Dazu werden Warteschlangen entsprechend angepasst,
         wobei dazu virtuelle und reale Spieler gesondert betrachtet werden. Weiter werden Spieler, die bereits an ihrem
         Ziel angekommen sind und dieses nun wieder verlassen, wieder sichtbar gemacht. Anschließend wird die Funktion
         "app.runback()" aufgerufen, die alle Positionen neu berechnet und entsprechend anpasst.
        :return: Kein Rückgabewert
        """
        for e in self.app.E:
            # füge Spieler wieder zu Queue hinzu
            for out in reversed(self.app.leaveQueue[self.zeitpunkt -1][self.app.E.index(e)]):
                try:
                    int(out)  # TypeError, falls "out" KEIN virtueller Spieler
                    self.app.deques[self.app.E.index(e)].append(out)  # virtuelle Spieler wieder zur Queue hinzufügen
                    self.app.spielerV[self.app.E.index(e)].append(patches.Rectangle((0.91*self.app.posit[e[0]][0] +
                                                                                     0.09*self.app.posit[e[1]][0],
                                                                                     0.91*self.app.posit[e[0]][1] +
                                                                                     0.09*self.app.posit[e[1]][1] +
                                                                                     self.app.rec_height *
                                                                        (len(self.app.deques[self.app.E.index(e)]) -1)),
                                                                                    self.app.rec_width,
                                                                                    self.app.rec_height, linewidth=1,
                                                                                    edgecolor='black',
                                                                                    facecolor='white', alpha=0.5))
                    self.app.ax.add_patch(self.app.spielerV[self.app.E.index(e)][-1])
                    continue
                except TypeError:  # nicht-virtueller Spieler
                    posit = self.app.pos(self.app.spieler.index(out), self.zeitpunkt -1)
                    if posit not in self.app.V:  # prüfe, ob Spieler zuvor in Warteschlange war
                        # falls ja, wird er wieder hinzugefügt
                        self.app.deques[self.app.E.index(posit[0])].append(out)
                        # Farbe wieder auf 'rot', da Spieler zurück in Queue
                        self.app.spieler[self.app.spieler.index(out)].set_facecolor('red')

        for e in self.app.E:
            while True:
                try:
                    # IndexError falls deque leer. Entferne Spieler aus Queue, falls sie diese zum Zeitpunkt
                    # "self.zeitpunkt" betreten haben
                    latest = self.app.deques[self.app.E.index(e)].popleft()
                    # TypeError, falls "latest" KEIN virtueller Spieler ist, ansonsten: prüft, ob virtueller Spieler aus
                    # Plot und "self.app.spielerV" entfernt werden muss
                    if int(latest) == self.zeitpunkt or int(latest) == self.zeitpunkt -1:
                        # entferne Spieler auch aus "self.app.spielerV" und aus Plot
                        self.app.remove_artist(self.app.E.index(e))
                        continue
                    # latest ist virtueller Spieler und schon vor "self.zeitpunkt" in deque gewesen
                    else:
                        # "latest" wird wieder zur deque hinzugefügt und es wird abgebrochen
                        self.app.deques[self.app.E.index(e)].appendleft(latest)
                        break
                except IndexError:  # deque leer
                    break
                except TypeError:  # "latest" ist normaler Spieler
                    # prüft, ob Spieler "latest" sich im Startknoten von "e" befindet, also die deque genau zu
                    # "self.zeitpunkt" betreten hat und entfernt werden muss
                    if self.app.pos(self.app.spieler.index(latest), self.zeitpunkt) == e[0]:
                        continue
                    # prüft, ob Spieler "latest" Kante "e" zum Zeitpunkt "self.zeitpunkt-1" betreten hat
                    elif self.app.fp[self.zeitpunkt-1][self.app.spieler.index(latest)] == e:
                        # wenn ja, umfärben von rot auf ursprüngliche Farbe und weiter mit nächstem Spieler
                        self.app.spieler[self.app.spieler.index(latest)].set_facecolor(
                            self.app.colors[self.app.ST[self.app.spieler.index(latest)][1]-1 % len(self.app.colors)])
                        continue
                    else:
                        self.app.deques[self.app.E.index(e)].appendleft(latest)  # wenn nein, Abbruch
                        break

        for i in self.app.I:
            if self.zeitpunkt == self.app.z[i]:
                self.app.spieler[i].set_visible(True)  # aktualisiere Sichtbarkeit und Abbruchkriterium
                self.app.numbers[i].set_visible(True)
                self.app.ankunft -= 1

        self.zeitpunkt -= 1
        if self.nex['state'] == 'disabled':
            self.nex.config(state="normal")
            self.after(1000, self.app.runner)  # runner neustarten, da Programm bereits terminiert hat
        self.app.runback(self.zeitpunkt)
        if self.zeitpunkt == 0:
            self.prev.config(state="disabled")
        return

    def weiter(self):
        """
        Bei Aufruf wird in "app" ein Zeitschritt weiter gegangen, d.h. "app.run()" aufgerufen.
        :return: Kein Rückgabewert
        """
        if self.app.ankunft >= self.app.num:
            self.nex.config(state="disabled")
            return 0
        if self.prev['state'] == 'disabled':
            self.prev.config(state="normal")
        # Vergrößere Liste "self.app.leaveQueue", falls nötig (also falls "self.app.run(self.zeitpunkt + 1)" vorher noch
        # nicht aufgerufen wurde)
        if len(self.app.leaveQueue) < self.zeitpunkt + 2:
            initp = []
            initm = []
            for i in self.app.I:
                initp.append(None)
                initm.append(None)
            self.app.fm.append(initm)
            self.app.fp.append(initp)
            self.app.leaveQueue.append([])
        self.zeitpunkt += 1
        self.app.run(self.zeitpunkt)
        return

    def get_unpaused(self):
        """
        :return: True, falls aktuell nicht pausiert, False sonst
        """
        return self.unpaused

    def set_unpaused(self, status):
        """
        :param status: neuer Wahrheitswert für "unpaused"
        :return: kein Rückgabewert
        """
        self.unpaused = status
        return

    def get_zeit(self):
        """
        :return: aktueller Zeitpunkt der Application
        """
        return self.zeitpunkt

    def set_zeit(self, t):
        """
        :param t: neuer Zeitpunkt der Application
        :return: kein Rückgabewert
        """
        self.zeitpunkt = t
        return
