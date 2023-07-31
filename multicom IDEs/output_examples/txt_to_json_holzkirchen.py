import os
import pickle
import sys
sys.path.append("..")
import data

graph = data.graph
u = data.u
coords = data.v_coords
items = graph.items()
keys = graph.keys()

path = "new_holzkirchen-8.json"
#assert os.path.isfile(path)

V = list(graph.keys())
E = []  # Kanten
nu = []  # Kapazitaeten
r = []
I = len(u)
for delta in items:
    for w in list(delta[1].keys()):
        E.append((delta[0], w))  # speichere alle Kanten in E
        r.append(list(items)[list(keys).index(delta[0])][1][w][1])  # speichere Reisezeiten in r
        nu.append(list(items)[list(keys).index(delta[0])][1][w][0])  # Kapazitäten in nu
m = len(E)
n = len(V)


def loadall(filename):
    with open(filename, "rb") as f:
        while True:
            try:
                yield pickle.load(f)
            except EOFError:
                break


items = loadall('new_holzkirchen_komplett-8.txt')

for (no, item) in enumerate(items):
    if no == 0:
        fp = item.copy()
    elif no == 1:
        fm = item.copy()
    elif no == 2:
        q_global = item.copy()
    elif no == 3:
        q_ind = item.copy()
    else:
        global_phase = item.copy()

output_json = open(path, "w")
output_json.write('{"network": {\n "nodes": [')
output_json.close()
output_json = open(path, "a")
for v_ind in range(1928):
    v = V[v_ind]
    output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(v_ind, coords[v][0], coords[v][1]))
output_json.write(' {{"id": {0}, "x": {1}, "y": {2}, "label": "{3}"}},'.format(1928, coords['t2'][0], coords['t2'][1], "t2"))

for v_ind in range(1929, 2169, 1):
    v = V[v_ind]
    output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(v_ind, coords[v][0], coords[v][1]))
output_json.write(' {{"id": {0}, "x": {1}, "y": {2}, "label": "{3}"}},'.format(2169, coords['t1'][0], coords['t1'][1], "t1"))

for v_ind in range(2169, n-1, 1):
    v = V[v_ind]
    output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(v_ind, coords[v][0], coords[v][1]))
output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}}'.format(n-1, coords[V[n-1]][0], coords[V[n-1]][1]))

output_json.write('], \n "edges": [')
for e_ind in range(m-1):
    v_ind = V.index(E[e_ind][0])
    w_ind = V.index(E[e_ind][1])
    output_json.write(' {{"id": {0}, "from": {1}, "to": {2}, "capacity": {3}, "transitTime": {4} }},'.format(e_ind, v_ind, w_ind, nu[e_ind], r[e_ind]))

v_ind = V.index(E[m-1][0])
w_ind = V.index(E[m-1][1])
output_json.write(' {{"id": {0}, "from": {1}, "to": {2}, "capacity": {3}, "transitTime": {4} }}'.format(m-1, v_ind, w_ind, nu[m-1], r[m-1]))

output_json.write('], \n "commodities": [')
colors = ["red", "blue", "green"]
for i in range(I-1):
    output_json.write(' {{ "id": {0}, "color": "{1}" }},'.format(i, colors[i]))
output_json.write(' {{ "id": {0}, "color": "{1}" }}'.format(I-1, colors[I-1]))

output_json.write('] }, \n "flow": { \n "inflow": [')
fpes = {}
for v_ind in range(n):
    fpes[v_ind] = {}
for e_ind in range(m):
    v_ind = V.index(E[e_ind][0])
    output_json.write('{')
    for i in range(I):
        output_json.write('"{0}": {{ \n "times": ['.format(i))
        for val in fp[i][e_ind][:-1]:
            if 29.42 < val[0] < 29.425:
                print()
            if val[1] > 0:
                if val[0] not in fpes[v_ind]:
                    fpes[v_ind][val[0]] = 1
                else:
                    fpes[v_ind][val[0]] += 1
            output_json.write(' {},'.format(val[0]))
        output_json.write(' {}'.format(fp[i][e_ind][-1][0]))

        output_json.write('], \n "values": [')
        for val in fp[i][e_ind][:-1]:
            output_json.write(' {},'.format(val[1]))
        output_json.write(' {}'.format(fp[i][e_ind][-1][1]))

        if i < I - 1:
            output_json.write('] },')
        else:
            output_json.write('] }')
    if e_ind < m - 1:
        output_json.write('},')
    else:
        output_json.write('}')

output_json.write('], \n "outflow": [')
for e_ind in range(m):
    output_json.write('{')
    for i in range(I):
        output_json.write('"{0}": {{ \n "times": ['.format(i))
        for val in fm[i][e_ind][:-1]:
            output_json.write(' {},'.format(val[0]))
        output_json.write(' {}'.format(fm[i][e_ind][-1][0]))

        output_json.write('], \n "values": [')
        for val in fm[i][e_ind][:-1]:
            output_json.write(' {},'.format(val[1]))
        output_json.write(' {}'.format(fm[i][e_ind][-1][1]))
        if i < I - 1:
            output_json.write('] },')
        else:
            output_json.write('] }')
    if e_ind < m - 1:
        output_json.write('},')
    else:
        output_json.write('}')

q_times = {}
q_vals = {}
for e_ind in range(m):
    q_times[e_ind] = []
    q_vals[e_ind] = []
for e_ind in range(m):
    len_qe = len(q_global[e_ind])
    for t in range(len_qe):
        q_times[e_ind].append(global_phase[q_ind[e_ind][t]])
        q_vals[e_ind].append(q_global[e_ind][t])
output_json.write('], \n "queues": [')
for e_ind in range(m):
    output_json.write('{ "times": [')
    for t in q_times[e_ind][:-1]:
        if 29.42 < t < 29.425:
            print()
        output_json.write(' {},'.format(t))
    output_json.write(' {}'.format(q_times[e_ind][-1]))

    output_json.write('], "values": [')
    for val in q_vals[e_ind][:-1]:
        output_json.write(' {},'.format(val))
    output_json.write(' {}'.format(q_vals[e_ind][-1]))

    if e_ind < m - 1:
        output_json.write('], \n "domain": ["-Infinity", "Infinity"], "lastSlope": 0.0, "firstSlope": 0.0},')
    else:
        output_json.write('], \n "domain": ["-Infinity", "Infinity"], "lastSlope": 0.0, "firstSlope": 0.0}')

output_json.write('] } }')
output_json.close()
