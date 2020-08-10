def concatenate(T,Z):
    """
    fügt einzelne Graphen zu einem zusammenhängenden Graphen zusammen
    erhalt als Eingabe zwei Listen T und Z
    kein Teilgraph darf einen Knoten enthalten, in dessen Name ein "_" vorkommt
    :param T: enthält alle Teilgraphen in Form von dictionarys, hat also die Form: T = [A,A,B,...,C], mit A,B,C
     gerichtete Graphen
    :param Z: enthält die Informationen, wie die Teilgraphen verbunden werden sollen und zwar auf folgende Weise:
        Z = [D1,...,Dk]. Z enthält für jeden Teilgraph in T ein dictionary. Das erste dictionary in Z korrespondiert zum
        ersten Teilgraphen in T, usw.. Für jeden Teilgraphen, der vom i-ten Teilgraphen über eine Kante erreichbar sein
        soll, existiert im i-ten dictionary von Z ein key, dessen Wert wieder ein dicitionary ist, nach folgendem
        Beispiel:
            D1 = {'2': {'(a,a)': (2,2)}, '4': {'(a,b)': (1,3), '(b,a)': (1,1)}}.
        Die Schlüssel (hier '2' und '4') bedeuten dabei, dass von dem aktuellen Teilgraph Kanten zum Teilgraph mit der
        Nummer des Schlüssels existieren (hier existieren also Kanten zum 2. und zum 4. Teilgraph in T). Zu jedem dieser
        Schlüssel gibt es als Wert ein weiteres Dictionary, welches alle Kanten vom aktuellen Teilgraph zum Zielgraph
        enthält (die Namen der Start- und Endknoten entsprechen den Namen innerhalb des jeweiligen Teilgraphen),
        zusammen mit den zugehörigen Kapazitäten und Reisezeiten.
        Ist also "D1" wie oben das erste Dicitonary in der Liste Z, so bedeutet dies, dass vom Knoten "a" im ersten
        Teilgraphen ("D1") eine Kante zum Knoten "a" im zweiten Teilgraphen ("D2") mit Kapazität und Reisedauer gleich 2
        exsitiert. Weiter existieren in diesem Beispiel zwei Kanten zwischen den Teilgraphen "D1" und "D4", nämlich die
        Kante von Knoten "a" (aus "D1") zu Knoten "b" (aus "D4") mit Kapazität 1 und Reisedauer 3, und die Kante von
        Knoten "b" (aus "D1) zu Knoten "a" (aus "D4") mit Kapazität und Reisedauer gleich 1.
    :return: graph: der zusammengefügte Graph als Dictionary
    """
    if len(T) != len(Z):
        raise TypeError('T und Z muessen die gleiche Laenge haben!')
    graph = {}
    for i in range(len(T)):
        knoten = list(T[i].keys())
        for v in knoten:
            insert={}
            for w in list(T[i][v].keys()):
                insert['{}_{}'.format(w,i+1)] = T[i][v][w]
            # füge alle Knoten und Kanten der einzelnen Teilgraphen zu "graph" hinzu (Resultat noch nicht
            # zusammenhängend)
            graph['{}_{}'.format(v,i+1)] = insert

    for dic in range(len(Z)):
        for i in list(Z[dic].keys()):
            ziel = Z[dic][i]
            kanten = list(ziel.keys())
            for e in kanten:  # füge Kanten zwischen einzelnen Zusammenhangskomponenten ein
                dazu = []
                start = ""
                for stelle in range(1,e.index(",")):
                    start += e[stelle]
                dazu.append('{}_{}'.format(start,dic+1))
                ende = ""
                for stelle in range(e.index(",") + 1,len(e)-1):
                    ende += e[stelle]
                dazu.append('{}_{}'.format(ende,i))
                dazu.append(ziel[e])
                graph[dazu[0]].update({dazu[1]: dazu[2]})

    quellen_nummer = 0
    senken_nummer = 0
    # Quellen und Senken müssen nun auf die passende Form für "main"-Funktion umbenannt werden, also ohne "_"
    for v in list(graph.keys()):
        name = str(v)
        if "_" not in name:
            continue
        is_quelle = False
        is_senke = False
        if name.startswith("s"):  # prüfe, ob aktueller Knoten Quellknoten
            if name.index("_") != 1:
                is_quelle = True
            for stelle in range(1, name.index("_")):
                try:
                    int(name[stelle])
                except:
                    is_quelle = False
                    break
        elif name.startswith("t"):  # prüfe, ob aktueller Knoten Senke
            if name.index("_") != 1:
                is_senke = True
            for stelle in range(1, name.index("_")):
                try:
                    int(name[stelle])
                except:
                    is_senke = False
                    break
        if is_quelle:
            quellen_nummer += 1
            graph['s{}'.format(quellen_nummer)] = graph.pop(v)  # Quellen werden umbenannt zu "s<quellen_nummer>"
            for start in graph.keys():
                if v in graph[start].keys():
                    # Umbenennung auch überall, wo eine Kante zur aktuellen Quelle existiert
                    graph[start]['s{}'.format(quellen_nummer)] = graph[start].pop(v)
        if is_senke:
            senken_nummer += 1
            graph['t{}'.format(senken_nummer)] = graph.pop(v)  # Senken werden umbenannt zu "t<senken_nummer>"
            for start in graph.keys():
                if v in graph[start].keys():
                    # Umbenennung auch überall, wo eine Kante zur aktuellen Senke existiert
                    graph[start]['t{}'.format(senken_nummer)] = graph[start].pop(v)

    return graph
