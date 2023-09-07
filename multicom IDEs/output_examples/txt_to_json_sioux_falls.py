import pickle
import sys
sys.path.append("..")
import data

graph = data.graph
u = data.u
items = graph.items()
keys = graph.keys()

path = "sioux_falls-5.json"

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

coords = {}
coords[0] = [-96.77041974, -43.61282792]
coords[1] = [-96.71125063, -43.60581298]
coords[2] = [-96.77430341, -43.5729616]
coords[3] = [-96.74716843, -43.56365362]
coords[4] = [-96.73156909, -43.56403357]
coords[5] = [-96.71164389, -43.58758553]
coords[6] = [-96.69342281, -43.5638436]
coords[7] = [-96.71138171, -43.56232379]
coords[8] = [-96.73124137, -43.54859634]
coords[9] = [-96.73143801, -43.54527088]
coords[10] = [-96.74684071, -43.54413068]
coords[11] = [-96.78013678, -43.54394065]
coords[12] = [-96.79337655, -43.49070718]
coords[13] = [-96.75103549, -43.52930613]
coords[14] = [-96.73150355, -43.52940117]
coords[15] = [-96.71138171, -43.54674361]
coords[16] = [-96.71138171, -43.54128009]
coords[17] = [-96.69407825, -43.54674361]
coords[18] = [-96.71131617, -43.52959125]
coords[19] = [-96.71118508, -43.5153335]
coords[20] = [-96.73097920, -43.51048509]
coords[21] = [-96.73124137, -43.51485818]
coords[22] = [-96.75090441, -43.51485818]
coords[23] = [-96.74920028, -43.50316422]


def loadall(filename):
    with open(filename, "rb") as f:
        while True:
            try:
                yield pickle.load(f)
            except EOFError:
                break


items = loadall('sioux_falls-5.txt')

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
for v_ind in range(23):
    v = V[v_ind]
    if v.startswith("t"):
        output_json.write(' {{"id": {0}, "x": {1}, "y": {2}, "label": "{3}"}},'.format(v_ind, coords[v_ind][0], coords[v_ind][1], v))
    else:
        output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(v_ind, coords[v_ind][0], coords[v_ind][1]))
output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}}'.format(23, coords[23][0], coords[23][1]))

output_json.write('], \n "edges": [')
for e_ind in range(m-1):
    v_ind = V.index(E[e_ind][0])
    w_ind = V.index(E[e_ind][1])
    output_json.write(' {{"id": {0}, "from": {1}, "to": {2}, "capacity": {3}, "transitTime": {4} }},'.format(e_ind, v_ind, w_ind, nu[e_ind], r[e_ind]))

v_ind = V.index(E[m-1][0])
w_ind = V.index(E[m-1][1])
output_json.write(' {{"id": {0}, "from": {1}, "to": {2}, "capacity": {3}, "transitTime": {4} }}'.format(m-1, v_ind, w_ind, nu[m-1], r[m-1]))

output_json.write('], \n "commodities": [')
colors = ["red", "blue", "green", "purple", "orange"]
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
