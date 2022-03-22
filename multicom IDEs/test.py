# init K: zulässige Menge durch [A,b], mit Ax = b äq zu x in K
I = 2
n = 7
m = 8

A_1 = np.zeros((I * (n-1), I * m))
A_2 = np.zeros((I, I * m))

b_1 = np.zeros(I * (n-1))
b_2 = np.zeros(I)

for i in range(I):
    V_i = self.V - ['t{}'.format(i+1)]
    for v in V_i:
        v_ind = self.V_i.index(v)
        b_1[i * (n-1) + v_ind] = self.b[i][self.V.index(v)]
        delta_v = self.get_outgoing_edges(v)  # Liste der Kantenindizes
        for e_ind in delta_v:
            A_1[i * (n-1) + v_ind][i * m + e_ind] = 1

    for v in self.V:
        delta_v_inact = self.get_outgoing_edges(v) - self.get_outgoing_active_edges(i, v)
        for e_ind in delta_v_inact:
            A_2[i][i * m + e_ind] = 1

