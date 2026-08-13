# BDR-Solver: Bitset-Accelerated Bipartite Degree Realization

This repository provides a Python implementation of the algorithmic framework
for the Bipartite Degree Realization (BDR) problem described in:

> Emrah Akyar and Handan Akyar,
> "A Bitset-Accelerated Algorithmic Framework for Bipartite Degree Realization",
> 2026.

## Overview

The Bipartite Degree Realization (BDR) problem asks whether a given
unpartitioned sequence of nonnegative integers can be realized as the degree
sequence of a simple bipartite graph.

For a prescribed bipartition, bipartite realizability is characterized by the
classical Gale--Ryser theorem. In the unpartitioned setting, however, a valid
bipartition must first be identified.

The implementation combines established graph-theoretic and algorithmic
ingredients in a unified framework for this unpartitioned problem:

1. **Bitset-accelerated partition search.**
   A balanced partition must satisfy
   `sum(X) = sum(Y) = sum(d)/2`. Bitset-based subset-sum computation is used to
   test this necessary condition efficiently.

2. **Suffix-reachability pruning and candidate validation.**
   Precomputed suffix bitsets are used during recursive search to discard
   branches that cannot reach the target sum. Candidate partitions are first
   checked against the partition-dependent degree bounds

   `max(X) <= |Y|` and `max(Y) <= |X|`,

   and are then validated using the Gale--Ryser criterion.

3. **Explicit graph construction.**
   Once a valid partition is found, a greedy bipartite realization procedure
   using a max-priority queue constructs an explicit simple bipartite graph.

The worst-case partition search remains exponential. The bitset and
suffix-reachability mechanisms are intended to reduce the number of search
states explored on instances for which these pruning conditions are effective.

## Key Features

- **Bitset-based subset-sum feasibility:** Python arbitrary-precision integers
  are used as compact bitsets for word-level subset-sum operations.
- **Suffix-reachability pruning:** branches that cannot reach the required
  balanced sum are discarded during recursive search.
- **Partition-dependent structural checks:** candidate partitions violating
  `max(X) <= |Y|` or `max(Y) <= |X|` are rejected before the full
  Gale--Ryser test.
- **Gale--Ryser validation:** candidate bipartitions are formally checked for
  bipartite realizability.
- **Explicit realization:** a heap-based greedy construction produces a
  realizing bipartite graph when a valid partition is found.
- **Visualization:** `NetworkX` and `Matplotlib` can be used to visualize the
  resulting graph.

## Complexity

Let `n` denote the length of the input sequence and let

`T = sum(d)/2`.

The suffix bitsets are computed in `O(nT/w)` word operations, where `w` is the
machine word size in a conventional bitset model. Storing all suffix bitsets
requires `O(nT)` bits, equivalently `O(nT/w)` machine words.

The recursive partition search has exponential worst-case complexity. Suffix
reachability and structural checks may substantially reduce the number of
visited states for particular instances, but they do not change the
worst-case exponential nature of the problem.

For a fixed candidate partition, the Gale--Ryser validation is implemented in
linear time after the inherited ordering of the subsequences is taken into
account. The explicit graph-construction phase uses a priority queue and
requires `O(m log q)` time, where `m` is the number of edges and `q = |Y|`.

## Requirements

Python 3.x and the following libraries are required:

- `networkx`
- `matplotlib`

Install them with:

```bash
pip install networkx matplotlib
```

The Python standard-library module `heapq` is used for the priority-queue
implementation and requires no additional installation.

## Usage

```python
from bdr_solver import BDRSolver

# Degree sequence from the illustrative case study
d = [5, 5, 3, 2, 2, 2, 1, 1, 1]

solver = BDRSolver(d)
result = solver.solve()

if result:
    X, Y = result

    print("Bipartite realization found!")
    print(f"Set X: {X}")
    print(f"Set Y: {Y}")

    G = solver.construct_graph(X, Y)
else:
    print("The sequence is not bipartite-realizable.")
```

For the sequence above, the solver may return the balanced bipartition

```text
X = [5, 5, 1]
Y = [3, 2, 2, 2, 1, 1]
```

which satisfies the partition-dependent degree bounds and the Gale--Ryser
criterion.

A solver is not required to return this particular partition if another valid
bipartition exists.

## Scope

The present implementation contains the bitset-based balanced-partition
search, suffix-reachability pruning, partition-dependent structural checks,
Gale--Ryser validation, and explicit graph construction.

Density-based tractability conditions, connectivity requirements, and
cut-constrained realization criteria discussed as possible extensions in the
paper are **not part of the current implementation**.

## Reproducibility

The source code in this repository accompanies the paper cited above and is
provided to facilitate reproducibility and experimentation with the proposed
algorithmic framework.

## Contact

**Emrah Akyar**  
Department of Mathematics  
Eskişehir Technical University  
Eskişehir, Türkiye

Email: eakyar@eskisehir.edu.tr

## License

This repository is distributed under the MIT License.
