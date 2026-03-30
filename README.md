# BDR-Solver: Bitset-Accelerated Bipartite Degree Realization

This repository provides an efficient Python implementation of the algorithmic framework for the **Bipartite Degree Realization (BDR)** problem, as described in the paper:

> **Emrah Akyar and Handan Akyar**, *"A Bitset-Accelerated Algorithmic Framework for Bipartite Degree Realization"*, 2026.

## Overview

The BDR-Solver determines whether a given unpartitioned degree sequence can be realized as a simple bipartite graph. The framework uses a three-phase approach:
1. **Phase 1:** Bitset-accelerated partition search with suffix reachability pruning.
2. **Phase 2:** Efficient Gale-Ryser validation and early-exit structural checks.
3. **Phase 3:** Explicit graph construction using a greedy bipartite Havel-Hakimi variant.

## Key Features

- **High Performance:** Utilizes word-level bitset parallelism for state-space exploration.
- **Smart Pruning:** Incorporates suffix-based reachability to bypass infeasible search branches.
- **Visualization:** Integrates with `NetworkX` and `Matplotlib` to visualize the resulting bipartite graphs.

## Requirements

To run the solver, you need Python 3.x and the following libraries:
- `networkx`
- `matplotlib`

You can install them via pip:
```bash
pip install networkx matplotlib
```

## Usage

```bash
from bdr_solver import BDRSolver

#Define a sequence
d = [5, 5, 3, 2, 2, 2, 1, 1, 1]

#Initialize and solve
solver = BDRSolver(d)
result = solver.solve()

if result:
    X, Y = result
    print(f"Partition found: X={X}, Y={Y}")
    G = solver.construct_graph(X, Y)
else:
    print("No realization exists.")
```

## Contact

**Emrah Akyar** Department of Mathematics, Eskişehir Technical University  
Email: [eakyar@eskisehir.edu.tr](mailto:eakyar@eskisehir.edu.tr)
