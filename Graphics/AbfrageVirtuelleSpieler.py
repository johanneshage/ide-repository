import tkinter as tk

"""
erzeugt Fenster zur Eingabe zusätzlicher virtueller Spieler zu Beginn des Programmaufrufs. Mit Klicken des Buttons
"Start" wird das Programm gestartet.
"""


class AbfrageVirtuelleSpieler(object):

    def __init__(self, app):
        """
        Erzeugt Fenster
        :param app: Objekt der Klasse "Application". Wird übergeben, um die Kantenkosten nach Einfügen neuer virtueller
         Spieler zu aktualisieren. Weiter wird "app" durch Klicken des "Start" Buttons gestartet, d.h. es werden
         "app.run(0)" und "app.button_win.set_unpaused(True)" aufgerufen.
        """
        self.app = app
        self.message = None

        def disable_event():  # Fenster soll nicht durch Klicken des 'X'-Buttons geschlossen werden können
            pass

        self.abfrage = tk.Toplevel(app.button_win)
        self.abfrage.attributes('-topmost', True)
        # Fenster kann nicht durch Klicken des 'X'-Buttons geschlossen werden
        self.abfrage.protocol("WM_DELETE_WINDOW", disable_event)
        self.abfrage.title("Einfügen virtueller Warteschlangen")
        self.abfrage.geometry("500x500")

        # Beschriftung Eingabemöglichkeiten
        tk.Label(self.abfrage, text="Kante: ").grid(row=0)
        tk.Label(self.abfrage, text="Startzeitpunkt: ").grid(row=1)
        tk.Label(self.abfrage, text="Endzeitpunkt: ").grid(row=2)
        tk.Label(self.abfrage, text="y0-Wert: ").grid(row=3)

        # Eingabemöglichkeiten
        self.eingabe_kante = tk.Entry(self.abfrage)
        self.eingabe_start = tk.Entry(self.abfrage)
        self.eingabe_ende = tk.Entry(self.abfrage)
        self.eingabe_y0 = tk.Entry(self.abfrage)

        # Anordnung
        self.eingabe_kante.grid(row=0, column=1)
        self.eingabe_start.grid(row=1, column=1)
        self.eingabe_ende.grid(row=2, column=1)
        self.eingabe_y0.grid(row=3, column=1)

        # Button zum Aufruf von self.add_parameter()
        tk.Button(self.abfrage, text='Werte hinzufügen', command=self.add_parameter).grid(row=4, column=0, sticky=tk.W,
                                                                                         pady=4)
        # Button zum Aufruf von start()
        tk.Button(self.abfrage, text='Start', command=self.start).grid(row=4, column=2, sticky=tk.W, pady=4)

    def add_parameter(self):
        """
        Aufruf bei Klicken des Buttons "Werte hinzufügen": Testet auf Gültigkeit der Eingaben und fügt virtuelle
         Warteschlangen hinzu. Für Gültigkeit wird beachtet:
            Kante, Startzeitpunkt und y0-Wert müssen angegeben werden (Bei fehlendem Endzeitpunkt, wird dieser gleich
            dem Startzeitpunkt gesetzt)
            Kante muss von der Form "(v,w)" sein (Leerzeichen werden nicht berücksichtigt
            Kante muss im Graph existieren (also in "app.E" enthalten sein)
            Es muss gelten: 0 <= Startzeitpunkt <= Endzeitpunkt, sowie y0 > 0 (insbesondere Werte Integer)
        :return: kein Rückgabewert
        """
        try:
            kante = self.eingabe_kante.get()
            start = int(self.eingabe_start.get())
            y0 = int(self.eingabe_y0.get())
        except:
            if self.message is not None:
                self.message.configure(text="Ungültige Eingabe: Fehlender Parameter")
            else:
                self.message = tk.Message(self.abfrage, text="Ungültige Eingabe: Fehlender Parameter", width=150)
                self.message.grid(row=5, column=1, columnspan=7)
            raise TypeError('Ungültige Eingabe')

        if len(self.eingabe_ende.get()) == 0:
            ende = start
        else:
            ende = int(self.eingabe_ende.get())

        e = kante.replace(" ", "")  # entferne Leerzeichen aus Eingabe
        if e[0] != "(" or e[-1] != ")" or "," not in e:  # prüfe, ob eingegebene Kante von der Form: "(v,w)" ist
            if self.message is not None:
                self.message.configure(text="Ungültige Kante!")
            else:
                self.message = tk.Message(self.abfrage, text="Ungültige Kante!", width=150)
                self.message.grid(row=5, column=1, columnspan=7)
            raise TypeError('Ungültige Kante')
        e = e[1:-1].split(",")
        edge = ('{}'.format(e[0]), '{}'.format(e[1]))
        if len(e) != 2 or edge not in self.app.E:  # prüfe, ob eingegebene Kante existiert
            if self.message is not None:
                self.message.configure(text="Ungültige Kante: Kante existiert nicht")
            else:
                self.message = tk.Message(self.abfrage, text="Ungültige Kante: Kante existiert nicht", width=150)
                self.message.grid(row=5, column=1, columnspan=7)
            raise TypeError('Ungültige Kante: Kante existiert nicht')
        try:
            if start > ende or ende < 0 or y0 <= 0 or start < 0:  # prüfe, ob Start-, End-, y0-Werte sinnvoll
                if self.message is not None:
                    self.message.configure(text="Ungültige Eingabewerte!")
                else:
                    self.message = tk.Message(self.abfrage, text="Ungültige Eingabewerte!", width=150)
                    self.message.grid(row=5, column=1, columnspan=7)
                raise TypeError('Ungültiger Zeitraum')
        except ValueError:
            if self.message is not None:
                self.message.configure(text="Ungültige Eingabewerte!")
            else:
                self.message = tk.Message(self.abfrage, text="Ungültige Eingabewerte!", width=150)
                self.message.grid(row=5, column=1, columnspan=7)
            raise TypeError('Ungültige Eingabewerte')

        self.app.kanten_queue.append(edge)  # füge neue Werte zu den entsprechenden Listen hinzu
        self.app.start_queue.append(start)
        self.app.ende_queue.append(ende)
        self.app.y0_queue.append(y0)

        if start == ende:
            if self.message is not None:
                self.message.configure(text="Virtuelle Warteschlange der Länge {} hinzugefügt für Kante {} zum "
                                              "Zeitpunkt {}".format(y0, edge, start))
            else:
                self.message = tk.Message(self.abfrage, text="Virtuelle Warteschlange der Länge {} hinzugefügt für "
                                                               "Kante {} zum Zeitpunkt {}".format(y0, edge, start),
                                          width=150)
                self.message.grid(row=5, column=1, columnspan=7)
        else:
            if self.message is not None:
                self.message.configure(text="Virtuelle Warteschlange der Länge {} hinzugefügt für Kante {} zu jedem "
                                              "Zeitpunkt von {} bis {}".format(y0, edge, start, ende))
            else:
                self.message = tk.Message(self.abfrage, text = "Virtuelle Warteschlange der Länge {} hinzugefügt für "
                                                               "Kante {} zu jedem Zeitpunkt von {} bis {}".format(y0,
                                                                                                                  edge,
                                                                                                                  start,
                                                                                                                  ende),
                                          width=150)
                self.message.grid(row=5, column=1, columnspan=7)

        # beginnt neue virtuelle Queue bei "theta" = 0, so muss diese bei der Initialisierung von "self.c" beachtet
        # werden
        if start == 0:
            for i in self.app.I:
                self.app.c[0][i][self.app.E.index(edge)] += \
                    (self.app.alpha[i]/(1 - self.app.alpha[i])) * y0/float(self.app.nu[self.app.E.index(edge)])

        self.eingabe_kante.delete(0,20)
        self.eingabe_start.delete(0,20)
        self.eingabe_ende.delete(0,20)
        self.eingabe_y0.delete(0,20)
        return

    def start(self):
        """
        Aufruf beim Klicken des Buttons "Start". Startet die Application "app".
        :return: kein Rückgabewert
        """
        self.app.button_win.prev.config(state="normal")
        self.app.button_win.pause.config(state="normal")
        self.app.button_win.nex.config(state="normal")
        self.app.run(0)  # Fluss zum Zeitpunkt 0 nach eventuellem Einfügen virtueller Spieler
        self.abfrage.destroy()
        self.app.button_win.set_unpaused(True)
        return
