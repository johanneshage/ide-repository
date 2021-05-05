import Math.cont_data as data
from Math.ContApp import ContApp

class ContMain:
    def __init__(self, graph, u):
        self.graph = graph
        self.u = u

        self.app = ContApp(self.graph, self.u)

def main():
    try:
        graph = data.graph
        u = data.u
    except AttributeError:
        raise AttributeError('Daten unvollständig!')
    ContMain(graph, u)
    return 0

main()
