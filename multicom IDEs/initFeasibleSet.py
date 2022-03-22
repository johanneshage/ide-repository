import polytope as pc

# Angenommen die Einflussraten u_i sind gegeben als Listen von Tupeln ("Startzeit", "Einflussrate") und somit stückweise
# linear.


# Funktion zur Berechnung von delta_v^+
def delta_plus(v):
    return list(self.G[v].keys())


def b_minus(i, v, theta):
    outneighbors = delta_minus(v)
    b = 0
    for w in outneighbors:
        e_ind = self.E.index((v, w))
        b += fm[i][e_ind]
