import os
import data
import pickle

graph = data.graph
u = data.u
items = graph.items()
keys = graph.keys()

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


items = loadall('non_term-6.txt')

for (no, item) in enumerate(items):
    if no == 0:
        fp = item.copy()
    elif no == 1:
        fm = item.copy()
    elif no == 2:
        q_global = item.copy()
    else:
        global_phase = item.copy()

knoten_no = 12 * 7 * 5 + 6 * 7 * 5 + 9 * 7 * 5
nB = list(range(knoten_no)) + list(range(960, 960 + knoten_no)) + [1920, 1921]
kanten_no = 27 * 12 * 5
esB = list(range(kanten_no)) + list(range(1635, 1635 + kanten_no))
for i in range(I):
    # B2
    for bno in range(5):
        ano_start = 12 * bno + 135*i
        for ano in range(ano_start, ano_start + 11):
            esB.remove(4 + ano * 12 + 15 * i)
            esB.remove(8 + ano * 12 + 15 * i)
            esB.remove(11 + ano * 12 + 15 * i)
        esB.remove(3 + 11 * 12 + 12 * 12 * bno + 135 * i * 12 + 15 * i)
        # esB.remove(4 + 11 * 12 + 12 * 12 * bno + 135 * i * 12 + 15 * i)
        esB.remove(8 + 11 * 12 + 12 * 12 * bno + 135 * i * 12 + 15 * i)
        esB.remove(11 + 11 * 12 + 12 * 12 * bno + 135 * i * 12 + 15 * i)
    # B5
    for bno in range(5):
        ano_start = 12*5 + 6 * bno + 135*i
        for ano in range(ano_start, ano_start + 5):
            esB.remove(4 + ano * 12 + 15 * i)
            esB.remove(8 + ano * 12 + 15 * i)
            esB.remove(11 + ano * 12 + 15 * i)
        esB.remove(3 + 5 * 12 + 12 * 12 * 5 + 12 * 6 * bno + 135 * i * 12 + 15 * i)
        # esB.remove(4 + 5 * 12 + 12 * 12 * 5 + 12 * 6 * bno + 135 * i * 12 + 15 * i)
        esB.remove(8 + 5 * 12 + 12 * 12 * 5 + 12 * 6 * bno + 135 * i * 12 + 15 * i)
        esB.remove(11 + 5 * 12 + 12 * 12 * 5 + 12 * 6 * bno + 135 * i * 12 + 15 * i)
    # B7
    for bno in range(5):
        ano_start = 12*5 + 6 * 5 + 9 * bno + 135*i
        for ano in range(ano_start, ano_start + 8):
            esB.remove(4 + ano * 12 + 15 * i)
            esB.remove(8 + ano * 12 + 15 * i)
            esB.remove(11 + ano * 12 + 15 * i)
        esB.remove(3 + 8 * 12 + 12 * 12 * 5 + 12 * 6 * 5 + 12 * 9 * bno + 135 * i * 12 + 15 * i)
        # esB.remove(4 + 8 * 12 + 12 * 12 * 5 + 12 * 6 * 5 + 12 * 9 * bno + 135 * i * 12 + 15 * i)
        esB.remove(8 + 8 * 12 + 12 * 12 * 5 + 12 * 6 * 5 + 12 * 9 * bno + 135 * i * 12 + 15 * i)
        esB.remove(11 + 8 * 12 + 12 * 12 * 5 + 12 * 6 * 5 + 12 * 9 * bno + 135 * i * 12 + 15 * i)

output_json = open("output-flow.json", "w")
output_json.write('{"network": {\n "nodes": [')
output_json.close()
output_json = open("output-flow.json", "a")
for i in range(I):
    # B2
    for bno in range(5):
            for ano in range(12):
                output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(0 + ano*7 + bno*12*7 + 945 * i, 0.0 + bno * 5 + i * 87, 0.0 - 3 * ano))
                output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(1 + ano*7 + bno*12*7 + 945 * i, 0.0 + bno * 5 + i * 87, -2.0 - 3*ano))
                output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(2 + ano*7 + bno*12*7 + 945 * i, -1.0 + bno * 5 + i * 87, -2.0 - 3*ano))
                output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(3 + ano*7 + bno*12*7 + 945 * i, -2.0 + bno * 5 + i * 87, -1.0 - 3*ano))
                output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(4 + ano*7 + bno*12*7 + 945 * i, -1.0 + bno * 5 + i * 87, 0.0 - 3*ano))
                output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(5 + ano*7 + bno*12*7 + 945 * i, 1.0 + bno * 5 + i * 87, -2.0 - 3*ano))
                output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(6 + ano*7 + bno*12*7 + 945 * i, 1.0 + bno * 5 + i * 87, 0.0 - 3*ano))
    # B5
    for bno in range(5):
        for ano in range(6):
            output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(420 + ano*7 + bno*6*7 + 945 * i, 25.0 + bno * 5 + i * 87, 0.0 - 3 * ano))
            output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(421 + ano*7 + bno*6*7 + 945 * i, 25.0 + bno * 5 + i * 87, -2.0 - 3*ano))
            output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(422 + ano*7 + bno*6*7 + 945 * i, 24.0 + bno * 5 + i * 87, -2.0 - 3*ano))
            output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(423 + ano*7 + bno*6*7 + 945 * i, 23.0 + bno * 5 + i * 87, -1.0 - 3*ano))
            output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(424 + ano*7 + bno*6*7 + 945 * i, 24.0 + bno * 5 + i * 87, 0.0 - 3*ano))
            output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(425 + ano*7 + bno*6*7 + 945 * i, 26.0 + bno * 5 + i * 87, -2.0 - 3*ano))
            output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(426 + ano*7 + bno*6*7 + 945 * i, 26.0 + bno * 5 + i * 87, 0.0 - 3*ano))
    # B7
    for bno in range(5):
        for ano in range(9):
            output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(630 + ano*7 + bno*9*7 + 945 * i, 50.0 + bno * 5 + i * 87, 0.0 - 3 * ano))
            output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(631 + ano*7 + bno*9*7 + 945 * i, 50.0 + bno * 5 + i * 87, -2.0 - 3*ano))
            output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(632 + ano*7 + bno*9*7 + 945 * i, 49.0 + bno * 5 + i * 87, -2.0 - 3*ano))
            output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(633 + ano*7 + bno*9*7 + 945 * i, 48.0 + bno * 5 + i * 87, -1.0 - 3*ano))
            output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(634 + ano*7 + bno*9*7 + 945 * i, 49.0 + bno * 5 + i * 87, 0.0 - 3*ano))
            output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(635 + ano*7 + bno*9*7 + 945 * i, 51.0 + bno * 5 + i * 87, -2.0 - 3*ano))
            output_json.write(' {{"id": {0}, "x": {1}, "y": {2}}},'.format(636 + ano*7 + bno*9*7 + 945 * i, 51.0 + bno * 5 + i * 87, 0.0 - 3*ano))
output_json.write(' {{"id": {0}, "x": {1}, "y": {2}, "label": "{3}"}},'.format(1890, 119.0, -42.0, "red"))
output_json.write(' {{"id": {0}, "x": {1}, "y": {2}, "label": "{3}"}},'.format(1891, 39.0, -42.0, "blue"))
output_json.close()
with open("output-flow.json", 'rb+') as oj:
    oj.seek(-1, os.SEEK_END)
    oj.truncate()
    oj.close()

output_json = open("output-flow.json", "a")
output_json.write('], \n "edges": [')
for e_ind in esB:
    v_ind = V.index(E[e_ind][0])
    w_ind = V.index(E[e_ind][1])
    output_json.write(' {{"id": {0}, "from": {1}, "to": {2}, "capacity": {3}, "transitTime": {4} }},'.format(esB.index(e_ind), nB.index(v_ind), nB.index(w_ind), nu[e_ind], r[e_ind]))
output_json.close()
with open("output-flow.json", 'rb+') as oj:
    oj.seek(-1, os.SEEK_END)
    oj.truncate()
    oj.close()

output_json = open("output-flow.json", "a")
output_json.write('], \n "commodities": [')
colors = ["red", "blue", "green"]
for i in range(I):
    output_json.write(' {{ "id": {0}, "color": "{1}" }},'.format(i, colors[i]))
output_json.close()
with open("output-flow.json", 'rb+') as oj:
    oj.seek(-1, os.SEEK_END)
    oj.truncate()
    oj.close()

output_json = open("output-flow.json", "a")
output_json.write('] }, \n "flow": { \n "inflow": [')
for e_ind in esB:
    output_json.write('{')
    for i in range(I):
        output_json.write('"{0}": {{ \n "times": ['.format(i))
        for val in fp[i][e_ind]:
            output_json.write(' {},'.format(val[0]))
        output_json.close()
        with open("output-flow.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()
        output_json = open("output-flow.json", "a")
        output_json.write('], \n "values": [')
        for val in fp[i][e_ind]:
            output_json.write(' {},'.format(val[1]))
        output_json.close()
        with open("output-flow.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()
        output_json = open("output-flow.json", "a")
        output_json.write('] },')
    output_json.close()
    with open("output-flow.json", 'rb+') as oj:
        oj.seek(-1, os.SEEK_END)
        oj.truncate()
        oj.close()
    output_json = open("output-flow.json", "a")
    output_json.write('},')
output_json.close()
with open("output-flow.json", 'rb+') as oj:
    oj.seek(-1, os.SEEK_END)
    oj.truncate()
    oj.close()

output_json = open("output-flow.json", "a")
output_json.write('], \n "outflow": [')
for e_ind in esB:
    output_json.write('{')
    for i in range(I):
        output_json.write('"{0}": {{ \n "times": ['.format(i))
        for val in fm[i][e_ind]:
            output_json.write(' {},'.format(val[0]))
        output_json.close()
        with open("output-flow.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()
        output_json = open("output-flow.json", "a")
        output_json.write('], \n "values": [')
        for val in fm[i][e_ind]:
            output_json.write(' {},'.format(val[1]))
        output_json.close()
        with open("output-flow.json", 'rb+') as oj:
            oj.seek(-1, os.SEEK_END)
            oj.truncate()
            oj.close()
        output_json = open("output-flow.json", "a")
        output_json.write('] },')
    output_json.close()
    with open("output-flow.json", 'rb+') as oj:
        oj.seek(-1, os.SEEK_END)
        oj.truncate()
        oj.close()
    output_json = open("output-flow.json", "a")
    output_json.write('},')
output_json.close()
with open("output-flow.json", 'rb+') as oj:
    oj.seek(-1, os.SEEK_END)
    oj.truncate()
    oj.close()

output_json = open("output-flow.json", "a")
q_times = {}
q_vals = {}
for e_ind in esB:
    q_times[e_ind] = [0]
    q_vals[e_ind] = [0]
len_q_global = len(q_global)
for t in range(len_q_global):
    for e_ind in esB:
        if q_global[t][0, e_ind] != q_vals[e_ind][-1]:
            q_times[e_ind].append(global_phase[t])
            q_vals[e_ind].append(q_global[t][0, e_ind])
output_json.write('], \n "queues": [')
for e_ind in esB:
    output_json.write('{ "times": [')
    for t in q_times[e_ind]:
        output_json.write(' {},'.format(t))
    output_json.close()
    with open("output-flow.json", 'rb+') as oj:
        oj.seek(-1, os.SEEK_END)
        oj.truncate()
        oj.close()

    output_json = open("output-flow.json", "a")
    output_json.write('], "values": [')
    for val in q_vals[e_ind]:
        output_json.write(' {},'.format(val))
    output_json.close()
    with open("output-flow.json", 'rb+') as oj:
        oj.seek(-1, os.SEEK_END)
        oj.truncate()
        oj.close()

    output_json = open("output-flow.json", "a")
    output_json.write('], \n "domain": ["-Infinity", "Infinity"], "lastSlope": 0.0, "firstSlope": 0.0},')

output_json.close()
with open("output-flow.json", 'rb+') as oj:
    oj.seek(-1, os.SEEK_END)
    oj.truncate()
    oj.close()
output_json = open("output-flow.json", "a")
output_json.write('] } }')
output_json.close()
