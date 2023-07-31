import pickle
import matplotlib.pyplot  as  plt
from matplotlib import ticker
import numpy as np


def loadall(filename):
    with open(filename, "rb") as f:
        while True:
            try:
                yield pickle.load(f)
            except EOFError:
                break


s_err = [[], []]
s_err_rel = [[], []]
items = loadall('output_examples/holzkirchen_komplett-8.txt')
for (no, item) in enumerate(items):
    if no == 0:
        fp = item.copy()

eps = 10**(-12)
rel_fps1 = {}

rel_fps1[0] = [fp[0][5595]]
rel_fps1[0].append(fp[0][5596])
rel_fps1[0].append(fp[0][5597])

rel_fps1[1] = [fp[1][5595]]
rel_fps1[1].append(fp[1][5596])
rel_fps1[1].append(fp[1][5597])



formatter = ticker.ScalarFormatter(useMathText=False)
formatter.set_scientific(True)
formatter.set_powerlimits((-1,1))

s_times095 = [t for (t, v) in rel_fps1[0][0]]
s_end_times095 = [t for (t, v) in rel_fps1[0][0][1:]] + [2]
nz_times095 = []
nz_times0952 = []
nz_vals095 = []
nz_vals0952 = []
for tvi in range(len(rel_fps1[0][0])):
    if rel_fps1[0][0][tvi][1] > eps:
        nz_times095.append(rel_fps1[0][0][tvi][0])
        nz_times095.append(rel_fps1[0][0][tvi+1][0])
        nz_times0952.append(rel_fps1[0][0][tvi][0])
        nz_vals095.append(rel_fps1[0][0][tvi][1])
        nz_vals095.append(rel_fps1[0][0][tvi][1])
        nz_vals0952.append(rel_fps1[0][0][tvi][1])
s_vals095 = [v for (t, v) in rel_fps1[0][0]]

s_times195 = [t for (t, v) in rel_fps1[1][0]]
s_end_times195 = [t for (t, v) in rel_fps1[1][0][1:]] + [2]
nz_times195 = []
nz_times1952 = []
nz_vals195 = []
nz_vals1952 = []
for tvi in range(len(rel_fps1[1][0])):
    if rel_fps1[1][0][tvi][1] > eps:
        nz_times195.append(rel_fps1[1][0][tvi][0])
        nz_times195.append(rel_fps1[1][0][tvi+1][0])
        nz_times1952.append(rel_fps1[1][0][tvi][0])
        nz_vals195.append(rel_fps1[1][0][tvi][1])
        nz_vals195.append(rel_fps1[1][0][tvi][1])
        nz_vals1952.append(rel_fps1[1][0][tvi][1])
s_vals195 = [v for (t, v) in rel_fps1[1][0]]

s_times096 = [t for (t, v) in rel_fps1[0][1]]
s_end_times096 = [t for (t, v) in rel_fps1[0][1][1:]] + [2]
nz_times096 = []
nz_times0962 = []
nz_vals096 = []
nz_vals0962 = []
for tvi in range(len(rel_fps1[0][1])):
    if rel_fps1[0][1][tvi][1] > eps:
        nz_times096.append(rel_fps1[0][1][tvi][0])
        nz_times096.append(rel_fps1[0][1][tvi+1][0])
        nz_times0962.append(rel_fps1[0][1][tvi][0])
        nz_vals096.append(rel_fps1[0][1][tvi][1])
        nz_vals096.append(rel_fps1[0][1][tvi][1])
        nz_vals0962.append(rel_fps1[0][1][tvi][1])
s_vals096 = [v for (t, v) in rel_fps1[0][1]]

s_times196 = [t for (t, v) in rel_fps1[1][1]]
s_end_times196 = [t for (t, v) in rel_fps1[1][1][1:]] + [2]
nz_times196 = []
nz_times1962 = []
nz_vals196 = []
nz_vals1962 = []
for tvi in range(len(rel_fps1[1][1])):
    if rel_fps1[1][1][tvi][1] > eps:
        nz_times196.append(rel_fps1[1][1][tvi][0])
        nz_times196.append(rel_fps1[1][1][tvi+1][0])
        nz_times1962.append(rel_fps1[1][1][tvi][0])
        nz_vals196.append(rel_fps1[1][1][tvi][1])
        nz_vals196.append(rel_fps1[1][1][tvi][1])
        nz_vals1962.append(rel_fps1[1][1][tvi][1])
s_vals196 = [v for (t, v) in rel_fps1[1][1]]

s_times097 = [t for (t, v) in rel_fps1[0][2]]
s_end_times097 = [t for (t, v) in rel_fps1[0][2][1:]] + [2]
nz_times097 = []
nz_times0972 = []
nz_vals097 = []
nz_vals0972 = []
for tvi in range(len(rel_fps1[0][2])):
    if rel_fps1[0][2][tvi][1] > eps:
        nz_times097.append(rel_fps1[0][2][tvi][0])
        nz_times097.append(rel_fps1[0][2][tvi+1][0])
        nz_times0972.append(rel_fps1[0][2][tvi][0])
        nz_vals097.append(rel_fps1[0][2][tvi][1])
        nz_vals097.append(rel_fps1[0][2][tvi][1])
        nz_vals0972.append(rel_fps1[0][2][tvi][1])
s_vals097 = [v for (t, v) in rel_fps1[0][2]]

s_times197 = [t for (t, v) in rel_fps1[1][2]]
s_end_times197 = [t for (t, v) in rel_fps1[1][2][1:]] + [2]
nz_times197 = []
nz_times1972 = []
nz_vals197 = []
nz_vals1972 = []
for tvi in range(len(rel_fps1[1][2])):
    if rel_fps1[1][2][tvi][1] > eps:
        nz_times197.append(rel_fps1[1][2][tvi][0])
        nz_times197.append(rel_fps1[1][2][tvi+1][0])
        nz_times1972.append(rel_fps1[1][2][tvi][0])
        nz_vals197.append(rel_fps1[1][2][tvi][1])
        nz_vals197.append(rel_fps1[1][2][tvi][1])
        nz_vals1972.append(rel_fps1[1][2][tvi][1])
s_vals197 = [v for (t, v) in rel_fps1[1][2]]

#plt.plot(s_times0, total_err, 'k')
fig, axs = plt.subplots(3)
fig.set_figheight(6)
# fig.suptitle('Vertically stacked subplots')

#axs[0].yaxis.set_major_formatter(formatter)
#axs[1].yaxis.set_major_formatter(formatter)
#axs[0].plot(s_times095, s_vals095, 'r')
#axs[0].plot(s_times195, s_vals195, 'b')

axs[0].hlines(s_vals095, s_times095, s_end_times095, colors='r')
axs[0].vlines(nz_times095, [0] + nz_vals095, nz_vals095[1:] + [0], colors='r', linestyles='dotted')
axs[0].hlines(s_vals195, s_times195, s_end_times195, colors='b')
axs[0].vlines(nz_times195, [0] + nz_vals195, nz_vals195[1:] + [0], colors='b', linestyles='dotted')
axs[0].scatter(nz_times0952 + [nz_times095[-1]], nz_vals0952 + [0], color='r', marker='.')
axs[0].scatter(nz_times1952[1:] + [nz_times195[-1]], nz_vals1952[1:] + [0], color='b', marker='.')

axs[1].hlines(s_vals096, s_times096, s_end_times096, colors='r')
axs[1].vlines(nz_times096, [0] + nz_vals096, nz_vals096[1:] + [0], colors='r', linestyles='dotted')
axs[1].scatter(nz_times0962 + [nz_times096[-1]], nz_vals0962 + [0], color='r', marker='.')
#axs[1].hlines(s_vals196, s_times196, s_end_times196, colors='b')
#axs[1].vlines(nz_times196, np.zeros(len(nz_vals196)), nz_vals196, colors='b', linestyles='dotted')

axs[2].hlines(s_vals097, s_times097, s_end_times097, colors='r')
axs[2].vlines(nz_times097, nz_vals097, nz_vals097[1:] + [0], colors='r', linestyles='dotted')
axs[2].hlines(s_vals197, s_times197, s_end_times197, colors='b')
axs[2].vlines(nz_times197, [0] + nz_vals197, nz_vals197[1:] + [0], colors='b', linestyles='dotted')
axs[2].scatter(nz_times0972[1:] + [nz_times097[-1]], nz_vals0972[1:] + [0], color='r', marker='.')
axs[2].scatter(nz_times1972 + [nz_times197[-1]], nz_vals1972 + [0], color='b', marker='.')

#axs[0].legend(['Gut 1', '', 'Gut 2'])
axs[0].yaxis.set_label_position("right")
axs[1].yaxis.set_label_position("right")
axs[2].yaxis.set_label_position("right")
axs[0].set(ylabel='Einfluss Kante 5595', ylim=(0, 15.5), xlim=0, yticks=np.arange(16, step=2), xticks=[0.135, 0.48433, 0.54526, 0.647, 0.731, 1.433, 1.723, 1.99333])
plt.setp(axs[0].get_xticklabels(), rotation=40, horizontalalignment='right', fontsize='x-small')
axs[1].set(ylabel='Einfluss Kante 5596', ylim=(0, 15.5), xlim=0, yticks=np.arange(16, step=2), xticks=[0.135, 0.48433, 0.54526, 0.647, 0.731, 1.433, 1.723, 1.99333])
plt.setp(axs[1].get_xticklabels(), rotation=40, horizontalalignment='right', fontsize='x-small')
axs[2].set(xlabel='Zeit', ylabel='Einfluss Kante 5597', ylim=(0, 15.5), xlim=0, yticks=np.arange(16, step=2), xticks=[0.135, 0.48433, 0.54526, 0.647, 0.731, 1.433, 1.723, 1.99333])
plt.setp(axs[2].get_xticklabels(), rotation=40, horizontalalignment='right', fontsize='x-small')
fig.tight_layout()
plt.show()
