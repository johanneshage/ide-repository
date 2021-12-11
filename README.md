# Instantaneous Dynamic Equilibria for dynamic Flows

This repository contains a calculation tool for Instantaneous Dynamic Equilibrium Flows (IDE-flows) for both the discrete and the continuous case.

For the definition and more theoretical information about IDE-flows, see [this](https://arxiv.org/abs/2007.07808).

## Running the program
### Discrete Case
In [Math/data.py](Math/data.py) the input data has to be specified. This includes:
* `R`: a list containing for each player the exact point in time where the player enters the network
* `ST`: a list containing tuples of length 2 according to the following example: the tuple `(i,j)` with index `k` in `ST` indicates that player `k+1` wants to travel from source <code>s<sub>i</sub></code> to sink <code>t<sub>j</sub></code>
* `alpha`: value in `[0,1]`. Describes the significance of the traveltime and the waiting time for the cost of an edge. For values towards `0` this is weighted towards the traveltime, for values closer to `1` the waiting time takes more impact. Default value: `0.5`, where both quantities take the same impact
* `kanten_queue`: a list containing all edges, where virtual players shall be added
* `start_queue`: for each edge in `kanten_queue` this list contains a point in time, where the virtual inflow into the edge starts
* `ende_queue`: for each edge in `kanten_queue` this list contains a point in time, where the virtual inflow into the edge stops
* `y0_queue`: for each edge in `kanten_queue` this list contains the amount of virtual inflow that will be added
* `graph`: directed graph as dictonary, if this is not specified, the program will read the graph from the file [GraphenGEXF/myGraph.gexf](GraphenGEXF/myGraph.gexf)

Once [Math/data.py](Math/data.py) contains all information, run [Math/Main.py](Math/Main.py), which will create a [Math/Application.py](Math/Application.py).

### Continuous Case
In [Math/cont_data.py](Math/cont_data.py) the input data has to be specified. This includes:
* `u`: a list of inflow rates. For each node this contains a list of tuples of length 2, where each tuple contains two values `<starting time>` and `<flow value>`, indicating an inflow of `<flow value>` flow volume, starting at `<starting time>` until the `<starting time>` of the next tuple in the list
* `graph`: directed graph as dictonary, if this is not specified, the program will read the graph from the file [GraphenGEXF/myGraph.gexf](GraphenGEXF/myGraph.gexf)
* `table_output`: Boolean value. If `True`, an `OutputTable` will be created, otherwise the program will return the raw data (see [Math/cont_data.py](Math/cont_data.py))

Once [Math/cont_data.py](Math/cont_data.py) contains all information, run [Math/ContMain.py](Math/ContMain.py), which will create a [Math/ContApp.py](Math/ContApp.py).

### Graph Editor
The [GraphEditor](sigma.js-1.2.1/GraphEditor) can be used to create `.gexf`-files for graphs which can then be read into the programs described above. It is based on `Sigma`, a JavaScript library dedicated to graph drawing, mainly developed by [@jacomyal](https://github.com/jacomyal) and [@Yomguithereal](https://github.com/Yomguithereal). See also [sigma.js-1.2.1/README.md](sigma.js-1.2.1/README.md).

To create graphs, open [sigma.js-1.2.1/GraphEditor/meinGraph.html](sigma.js-1.2.1/GraphEditor/meinGraph.html). Once the graph is created, its `.gexf`-file can be downloaded via the `Download`-Button. Saving it in [GraphenGEXF/myGraph.gexf](GraphenGEXF/myGraph.gexf) will result in the graph being read into the program.
