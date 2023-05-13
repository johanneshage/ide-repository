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
items = loadall('output_examples/no_term_errors250-6.txt')
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
# fig.suptitle('Vertically stacked subplots')

axs[0].yaxis.set_major_formatter(formatter)
axs[1].yaxis.set_major_formatter(formatter)
axs[0].plot(s_times0, s_vals0, 'r')
axs[0].plot(s_times0, s_vals1, 'b')
axs[1].plot(s_times0, rel_err, 'c')
axs[1].plot(s_times0, total_err, 'k')
axs[1].legend(['relativ', 'absolut'])
axs[0].set(ylabel='Fehler', ylim=0, xlim=0)
axs[1].set(xlabel='Zeit', ylabel='Fehler', ylim=0, xlim=0)
plt.show()
