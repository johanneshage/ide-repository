import pickle
import matplotlib.pyplot  as  plt
from matplotlib import ticker


def loadall(filename):
    with open(filename, "rb") as f:
        while True:
            try:
                yield pickle.load(f)
            except EOFError:
                break


s_err = [[], []]
s_err_rel = [[], []]
items = loadall('output_examples/holzkirchen_final_err5-8.txt')
for (no, item) in enumerate(items):
    it = item.copy()
    if no % 2:
        for (t, v) in it[0]:
            s_err_rel[0].append((t, v))
        for (t, v) in it[1]:
            s_err_rel[1].append((t, v))
    else:
        for (t, v) in it[0]:
            s_err[0].append((t, v))
        for (t, v) in it[1]:
            s_err[1].append((t, v))


formatter = ticker.ScalarFormatter(useMathText=False)
formatter.set_scientific(True)
formatter.set_powerlimits((-1,1))
s_times0 = [t for (t, v) in s_err[0]]
s_vals0 = [v for (t, v) in s_err[0]]
len_s_err = len(s_err[0])
#plt.plot(s_times0, s_vals0, 'r')
s_times1 = [t for (t, v) in s_err[1]]
s_vals1 = [v for (t, v) in s_err[1]]
#plt.plot(s_times1, s_vals1, 'b')
s_vals_rel0 = [v for (t, v) in s_err_rel[0]]
s_vals_rel1 = [v for (t, v) in s_err_rel[1]]

total_err = []
rel_err = []
for t in range(len(s_times0)):
    total_err.append(s_vals0[t] + s_vals1[t])
    rel_err.append(s_vals_rel0[t] + s_vals_rel1[t])
#plt.plot(s_times0, total_err, 'k')
fig, axs = plt.subplots(2)
#fig.set_figheight(6)
#fig.set_figwidth(8)

axs[0].yaxis.set_major_formatter(formatter)
axs[1].yaxis.set_major_formatter(formatter)
# axs[0].plot(s_times0, s_vals0, 'r')
for ind in range(0, len_s_err, 2):
    axs[0].plot(s_times0[ind+1:ind+3], s_vals0[ind+1:ind+3], 'r')
    axs[0].plot(s_times0[ind:ind+2], s_vals0[ind:ind+2], 'r:')
    #axs[0].scatter(s_times0[ind+1], s_vals0[ind+1], color='r')
    axs[0].plot(s_times0[ind+1:ind+3], s_vals1[ind+1:ind+3], 'b')
    axs[0].plot(s_times0[ind:ind+2], s_vals1[ind:ind+2], 'b:')
    axs[1].plot(s_times0[ind+1:ind+3], total_err[ind+1:ind+3], 'k')
    axs[1].plot(s_times0[ind+1:ind+3], rel_err[ind+1:ind+3], 'c')
    axs[1].plot(s_times0[ind:ind+2], total_err[ind:ind+2], 'k:')
    axs[1].plot(s_times0[ind:ind+2], rel_err[ind:ind+2], 'c:')

axs[1].legend(['absolut', 'relativ'])
axs[0].set(ylabel='maximaler IDE-Fehler', ylim=(0, 1.2 * 10**(-8)), xlim=(0, 2.2), xticks=[0.135, 0.389, 0.48433, 0.5925, 0.731, 1.433, 1.56233, 1.723, 1.99333])
plt.setp(axs[0].get_xticklabels(), rotation=40, horizontalalignment='right', fontsize='x-small')
axs[1].set(xlabel='Zeit', ylabel='IDE-Gesamtfehler', ylim=(0, 1.7 * 10**(-8)), xlim=(0, 2.2), xticks=[0.135, 0.389, 0.48433, 0.54526, 0.647, 0.7474, 1.433, 1.56233, 1.723, 1.99333])
plt.setp(axs[1].get_xticklabels(), rotation=40, horizontalalignment='right', fontsize='x-small')
fig.tight_layout()
plt.show()
