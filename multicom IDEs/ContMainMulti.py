import data
from ContAppMultiOT import ContAppMultiOT
from ContAppMulti import ContAppMulti


class ContMainMulti:
    def __init__(self, graph, u):
        self.graph = graph
        self.u = u

        # self.app = ContAppMultiOT(self.graph, self.u)
        self.app = ContAppMulti(self.graph, self.u)


def main():
    try:
        graph = data.graph
        u = data.u
    except AttributeError:
        raise AttributeError('Daten unvollständig!')
    ContMainMulti(graph, u)
    return 0


main()
