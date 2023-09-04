import pickle
import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np


def loadall(filename):
    with open(filename, "rb") as f:
        while True:
            try:
                yield pickle.load(f)
            except EOFError:
                break


s_err = [[], [], []]
s_err_rel = [[], [], []]
l_err_max = [[(0,0)], [(0,0)], [(0,0)]]
l_err_min = [[(0,0)], [(0,0)], [(0,0)]]
items = loadall('output_examples/bsp1_errors-5.txt')
for (no, item) in enumerate(items):
    it = item.copy()
    if no % 4 == 1:
        for (t, v) in it[0]:
            #if len(s_err_rel[0]) > 1 and v == 0 and s_err_rel[0][-1][1] == 0 and s_err_rel[0][-2][1] == 0:
            #    continue
            s_err_rel[0].append((t, v))
        for (t, v) in it[1]:
            #if len(s_err_rel[1]) > 1 and v == 0 and s_err_rel[1][-1][1] == 0 and s_err_rel[1][-2][1] == 0:
            #    continue
            s_err_rel[1].append((t, v))
        for (t, v) in it[2]:
            #if len(s_err_rel[2]) > 1 and v == 0 and s_err_rel[2][-1][1] == 0 and s_err_rel[2][-2][1] == 0:
            #    continue
            s_err_rel[2].append((t, v))
    elif no % 4 == 2:
        for (t, v) in it[0]:
            l_err_max[0].append((t, v))
        for (t, v) in it[1]:
            l_err_max[1].append((t, v))
        for (t, v) in it[2]:
            l_err_max[2].append((t, v))
    elif no % 4 == 3:
        for (t, v) in it[0]:
            l_err_min[0].append((t, v))
        for (t, v) in it[1]:
            l_err_min[1].append((t, v))
        for (t, v) in it[2]:
            l_err_min[2].append((t, v))
    else:
        for (t, v) in it[0]:
            #if len(s_err[0]) > 1 and v == 0 and s_err[0][-1][1] == 0 and s_err[0][-2][1] == 0:
            #    continue
            s_err[0].append((t, v))
        for (t, v) in it[1]:
            #if len(s_err[1]) > 1 and v == 0 and s_err[1][-1][1] == 0 and s_err[1][-2][1] == 0:
            #    continue
            s_err[1].append((t, v))
        for (t, v) in it[2]:
            #if len(s_err[2]) > 1 and v == 0 and s_err[2][-1][1] == 0 and s_err[2][-2][1] == 0:
            #    continue
            s_err[2].append((t, v))


formatter = ticker.ScalarFormatter(useMathText=False)
formatter.set_scientific(True)
formatter.set_powerlimits((-1,1))
s_times0 = [t for (t, v) in s_err[0]]
s_vals0 = [v for (t, v) in s_err[0]]
len_s_err0 = len(s_err[0])
len_s_err1 = len(s_err[1])
len_s_err2 = len(s_err[2])
#plt.plot(s_times0, s_vals0, 'r')
s_times1 = [t for (t, v) in s_err[1]]
s_vals1 = [v for (t, v) in s_err[1]]
s_times2 = [t for (t, v) in s_err[2]]
s_vals2 = [v for (t, v) in s_err[2]]
s_vals_rel0 = [v for (t, v) in s_err_rel[0]]
s_vals_rel1 = [v for (t, v) in s_err_rel[1]]
s_vals_rel2 = [v for (t, v) in s_err_rel[2]]
l_times_max0 = [t for (t, v) in l_err_max[0]]
l_times_min0 = [t for (t, v) in l_err_min[0]]
l_max_vals0 = [v for (t, v) in l_err_max[0]]
l_min_vals0 = [v for (t, v) in l_err_min[0]]
l_max_vals1 = [v for (t, v) in l_err_max[1]]
l_min_vals1 = [v for (t, v) in l_err_min[1]]
l_max_vals2 = [v for (t, v) in l_err_max[2]]
l_min_vals2 = [v for (t, v) in l_err_min[2]]

total_err = []
rel_err = []
times = list(set(s_times0 + s_times1 + s_times2))
times.sort()
all_times = []
t_jump = [[], [], []]
for t in times:
    total_v = 0
    total_v_rel = 0
    if t in s_times0:
        t0_ind = s_times0.index(t)
        total_v += s_vals0[t0_ind]
        total_v_rel += s_vals_rel0[t0_ind]
        if len(s_times0) > t0_ind + 1 and s_times0[t0_ind + 1] == t:
            t_jump[0].append(t)
    if t in s_times1:
        t1_ind = s_times1.index(t)
        total_v += s_vals1[t1_ind]
        total_v_rel += s_vals_rel1[t1_ind]
        if len(s_times1) > t1_ind + 1 and s_times1[t1_ind + 1] == t:
            t_jump[1].append(t)
    if t in s_times2:
        t2_ind = s_times2.index(t)
        total_v += s_vals2[t2_ind]
        total_v_rel += s_vals_rel2[t2_ind]
        if len(s_times2) > t2_ind + 1 and s_times2[t2_ind + 1] == t:
            t_jump[2].append(t)
    total_err.append(total_v)
    rel_err.append(total_v_rel)
    all_times.append(t)
    if np.any([t in t_jump[i] for i in range(3)]):
        if t in t_jump[0]:
            total_v += s_vals0[t0_ind + 1] - s_vals0[t0_ind]
            total_v_rel += s_vals_rel0[t0_ind + 1] - s_vals_rel0[t0_ind]
        if t in t_jump[1]:
            total_v += s_vals1[t1_ind + 1] - s_vals1[t1_ind]
            total_v_rel += s_vals_rel1[t1_ind + 1] - s_vals_rel1[t1_ind]
        if t in t_jump[2]:
            total_v += s_vals2[t2_ind + 1] - s_vals2[t2_ind]
            total_v_rel += s_vals_rel2[t2_ind + 1] - s_vals_rel2[t2_ind]
        total_err.append(total_v)
        rel_err.append(total_v_rel)
        all_times.append(t)

#plt.plot(s_times0, total_err, 'k')
plt.rcParams['svg.fonttype'] = 'none'
fig, axs = plt.subplots(3)
fig.set_figheight(8)
fig.set_figwidth(8)
# fig.suptitle('Vertically stacked subplots')

axs[0].yaxis.set_major_formatter(formatter)
axs[1].yaxis.set_major_formatter(formatter)
axs[2].yaxis.set_major_formatter(formatter)
# axs[0].plot(s_times0, s_vals0, 'r')
axs[0].plot([0, s_times0[0]], [0, s_vals0[0]], 'r')
axs[0].plot([0, s_times1[0]], [0, s_vals1[0]], 'b')
axs[0].plot([0, s_times2[0]], [0, s_vals2[0]], 'g')
for ind in range(0, len_s_err0, 2):
    axs[0].plot(s_times0[ind+1:ind+3], s_vals0[ind+1:ind+3], 'r')
    axs[0].plot(s_times0[ind:ind+2], s_vals0[ind:ind+2], 'r:')
    #axs[0].scatter(s_times0[ind+1], s_vals0[ind+1], color='r')
for ind in range(0, len_s_err1, 2):
    axs[0].plot(s_times1[ind+1:ind+3], s_vals1[ind+1:ind+3], 'b')
    axs[0].plot(s_times1[ind:ind+2], s_vals1[ind:ind+2], 'b:')
for ind in range(0, len_s_err2, 2):
    axs[0].plot(s_times2[ind+1:ind+3], s_vals2[ind+1:ind+3], 'g')
    axs[0].plot(s_times2[ind:ind+2], s_vals2[ind:ind+2], 'g:')

len_times = len(all_times)
if len_times % 2:
    len_times -= 1
axs[1].plot([0, all_times[0]], [0, total_err[0]], 'k')
axs[1].plot([0, all_times[0]], [0, rel_err[0]], 'c')
for ind in range(0, len_times, 2):
    #if np.any([all_times[ind] in t_jump[i] for i in range(2)]):
    #    axs[1].plot(all_times[ind:ind+2], total_err[ind:ind+2], 'k:')
    #else:
    #    axs[1].plot(all_times[ind:ind+2], total_err[ind:ind+2], 'k')
    axs[1].plot(all_times[ind+1:ind+3], total_err[ind+1:ind+3], 'k')
    axs[1].plot(all_times[ind+1:ind+3], rel_err[ind+1:ind+3], 'c')
    axs[1].plot(all_times[ind:ind+2], total_err[ind:ind+2], 'k:')
    axs[1].plot(all_times[ind:ind+2], rel_err[ind:ind+2], 'c:')

l_times = l_times_max0.copy()
len_l_times_max = len(l_times_max0)
len_l_times_min = len(l_times_min0)
axs[2].plot(l_times[0:2], l_max_vals0[0:2], 'darkred')
axs[2].plot(l_times[0:2], l_min_vals0[0:2], 'tomato')
axs[2].plot(l_times[0:2], l_max_vals2[0:2], 'green')
axs[2].plot(l_times[0:2], l_min_vals2[0:2], 'limegreen')
axs[2].plot(l_times[1:3], l_max_vals0[1:3], color='darkred', linestyle='dotted')
axs[2].plot(l_times[1:3], l_min_vals0[1:3], color='tomato', linestyle='dotted')
axs[2].plot(l_times[1:3], l_max_vals2[1:3], color='green', linestyle='dotted')
axs[2].plot(l_times[1:3], l_min_vals2[1:3], color='limegreen', linestyle='dotted')
for ind in range(2, len_l_times_max, 2):
    axs[2].plot(l_times[ind:ind+2], l_max_vals0[ind:ind+2], 'darkred')
    axs[2].plot(l_times[ind+1:ind+3], l_max_vals0[ind+1:ind+3], color='darkred', linestyle='dotted')
    axs[2].plot(l_times[ind:ind+2], l_max_vals2[ind:ind+2], 'green')
    axs[2].plot(l_times[ind+1:ind+3], l_max_vals2[ind+1:ind+3], color='green', linestyle='dotted')
for ind in range(2, len_l_times_min, 2):
    axs[2].plot(l_times[ind:ind+2], l_min_vals0[ind:ind+2], 'tomato')
    axs[2].plot(l_times[ind+1:ind+3], l_min_vals0[ind+1:ind+3], color='tomato', linestyle='dotted')
    axs[2].plot(l_times[ind:ind+2], l_min_vals2[ind:ind+2], 'limegreen')
    axs[2].plot(l_times[ind+1:ind+3], l_min_vals2[ind+1:ind+3], color='limegreen', linestyle='dotted')

axs[0].legend(['Gut 1', 'Gut 2', 'Gut 3'])
axs[1].legend(['absolut', 'relativ'])
axs[2].legend(['Gut $1^+$', 'Gut $1^-$', 'Gut $3^+$', 'Gut $3^-$'])
#axs[2].legend(['$\max \;\lambda_{i,v} - \lambda^{*}_{i, v}$', '$\min \;\lambda_{i,v} - \lambda^{*}_{i, v}$'])
axs[0].set(ylabel='maximaler IDE-Fehler', xlim=(0,5.5))#, ylim=(0, 1 * 10**(-8)))# , ylim=(0, 1.2 * 10**(-8)), xlim=(0, 2.2), xticks=[0.135, 0.389, 0.48433, 0.5925, 0.731, 1.433, 1.56233, 1.723, 1.99333])
#plt.setp(axs[0].get_xticklabels(), rotation=40, horizontalalignment='right', fontsize='x-small')
axs[1].set(xlim=(0,5.5), ylabel='IDE-Gesamtfehler')#, ylim=(0, 1.2 * 10**(-8)), xlim=(0, 2.2), xticks=[0.135, 0.389, 0.48433, 0.54526, 0.647, 0.7474, 1.433, 1.56233, 1.723, 1.99333])
#plt.setp(axs[1].get_xticklabels(), rotation=40, horizontalalignment='right', fontsize='x-small')
axs[2].set(xlabel='Zeit', xlim=(0, 5.5), ylabel='Fehler im Knotenlabel')
fig.tight_layout()
#plt.savefig("tast.svg")
plt.show()
