import heapq
import networkx as nx
import matplotlib.pyplot as plt

class BDRSolver:
    """
    Implements the Bitset-Accelerated Algorithmic Framework for 
    Bipartite Degree Realization (BDR) as described in the paper.
    """
    def __init__(self, degree_sequence):
        # Ensure the sequence is sorted in nonincreasing order
        self.d = sorted(degree_sequence, reverse=True)
        self.n = len(self.d)
        self.S = sum(self.d)
        self.T = self.S // 2
        self.suffix_bitsets = []

    def precompute_suffix_bitsets(self):
        """
        Phase 1: Precomputes suffix bitsets to enable constant-time 
        reachability pruning during the recursive search.
        """
        if self.S % 2 != 0:
            return False
        
        # Initialize bitsets B_1 to B_{n+1}
        # Python's arbitrary-precision integers act as efficient bitsets
        self.suffix_bitsets = [0] * (self.n + 1)
        self.suffix_bitsets[self.n] = 1  # Base case: sum 0 is always reachable
        
        for i in range(self.n - 1, -1, -1):
            # Transition: B_i = B_{i+1} OR (B_{i+1} << d_i)
            self.suffix_bitsets[i] = self.suffix_bitsets[i+1] | (self.suffix_bitsets[i+1] << self.d[i])
            # Masking bits beyond target T to optimize memory
            mask = (1 << (self.T + 1)) - 1
            self.suffix_bitsets[i] &= mask
            
        # Check if the total target T is reachable at all
        return (self.suffix_bitsets[0] >> self.T) & 1

    def structural_check(self, X, Y):
        """
        Phase 2: O(1) partition-dependent degree filtering.

        The candidate subsequences X and Y inherit the nonincreasing order
        of self.d, so their maximum entries are X[0] and Y[0]. A simple
        bipartite realization requires max(X) <= |Y| and max(Y) <= |X|.
        """
        if not X or not Y:
            return False
        return X[0] <= len(Y) and Y[0] <= len(X)

    def gale_ryser_check(self, X, Y):
        """
        Phase 2: Formal verification using the Gale-Ryser majorization criteria.
        Checks if X is majorized by the conjugate sequence Y*.
        """
        X = sorted(X, reverse=True)
        Y = sorted(Y, reverse=True)
        
        if sum(X) != sum(Y):
            return False
            
        # Compute the conjugate sequence Y* in linear time.
        # Because Y is nonincreasing, a single moving boundary gives
        # y_k^* = |{j : y_j >= k}| for k = 1, ..., y_1.
        y_max = Y[0] if Y else 0
        y_star = []
        count = len(Y)
        for k in range(1, y_max + 1):
            while count > 0 and Y[count - 1] < k:
                count -= 1
            y_star.append(count)
        
        # Majorization check: Prefix sums of X <= Prefix sums of Y*
        prefix_sum_x = 0
        prefix_sum_ystar = 0
        
        # Extend sequences for comparison
        max_len = max(len(X), len(y_star))
        x_ext = X + [0] * (max_len - len(X))
        y_s_ext = y_star + [0] * (max_len - len(y_star))
        
        for x_val, ys_val in zip(x_ext, y_s_ext):
            prefix_sum_x += x_val
            prefix_sum_ystar += ys_val
            if prefix_sum_x > prefix_sum_ystar:
                return False
        return True

    def solve(self):
        """
        Main entry point for the unpartitioned search.
        Includes symmetry breaking by fixing d_1 in X.
        """
        if not self.precompute_suffix_bitsets():
            return None
        
        # Symmetry breaking: Fix the first element d[0] in set X
        return self._search(1, self.d[0], [0])

    def _search(self, idx, current_sum, x_indices):
        """
        Recursive backtracking with Bitset-Based Suffix Reachability Pruning.
        """
        # Base case: Target sum reached
        if current_sum == self.T:
            X = [self.d[i] for i in x_indices]
            y_indices = [i for i in range(self.n) if i not in x_indices]
            Y = [self.d[i] for i in y_indices]
            
            # Phase 2: Validation
            if self.structural_check(X, Y) and self.gale_ryser_check(X, Y):
                return X, Y
            return None

        # Pruning: Bounds and Reachability Check
        if idx >= self.n or current_sum > self.T:
            return None

        needed = self.T - current_sum
        # Direct reachability bit test using the precomputed suffix bitset
        if not ((self.suffix_bitsets[idx] >> needed) & 1):
            return None

        # Branch 1: Include d[idx] in X
        if current_sum + self.d[idx] <= self.T:
            res = self._search(idx + 1, current_sum + self.d[idx], x_indices + [idx])
            if res: return res

        # Branch 2: Include d[idx] in Y (Exclude from X)
        return self._search(idx + 1, current_sum, x_indices)

    def construct_graph(self, X, Y):
        """
        Phase 3: Explicit graph construction using a max-priority queue.

        Since Python's heapq is a min-heap, negative residual degrees are
        stored so that the vertex of largest remaining degree is popped first.
        """
        G = nx.Graph()
        U_nodes = [f"u{i}" for i in range(len(X))]
        V_nodes = [f"v{j}" for j in range(len(Y))]
        G.add_nodes_from(U_nodes, bipartite=0)
        G.add_nodes_from(V_nodes, bipartite=1)

        # Max-heap on positive residual degrees in V.
        heap = [(-Y[j], V_nodes[j]) for j in range(len(Y)) if Y[j] > 0]
        heapq.heapify(heap)

        for i, deg_u in enumerate(X):
            if deg_u > len(heap):
                raise ValueError("Construction failed: insufficient positive residual degrees in V.")

            selected = []
            for _ in range(deg_u):
                neg_deg_v, v_node = heapq.heappop(heap)
                deg_v = -neg_deg_v
                G.add_edge(U_nodes[i], v_node)
                if deg_v - 1 > 0:
                    selected.append((-(deg_v - 1), v_node))

            # Reinsert only after all neighbors of u_i have been selected,
            # preventing multiple edges from u_i to the same vertex.
            for item in selected:
                heapq.heappush(heap, item)

        if heap:
            raise ValueError("Construction failed: positive residual degrees remain in V.")

        return G

# --- Execution Example ---
if __name__ == "__main__":
    # Test sequence from the paper's Illustrative Case Study
    test_sequence = [5, 5, 3, 2, 2, 2, 1, 1, 1]
    
    solver = BDRSolver(test_sequence)
    result = solver.solve()

    if result:
        partition_X, partition_Y = result
        print("Bipartite realization found!")
        print(f"Set X: {partition_X}")
        print(f"Set Y: {partition_Y}")
        
        # Generate the graph
        realization = solver.construct_graph(partition_X, partition_Y)
        
        # Visualization using NetworkX and Matplotlib
        u_side = [n for n, d in realization.nodes(data=True) if d['bipartite']==0]
        pos = nx.bipartite_layout(realization, u_side)
        
        plt.figure(figsize=(10, 7))
        nx.draw(realization, pos, with_labels=True, 
                node_color=['skyblue' if n.startswith('u') else 'lightgreen' for n in realization.nodes()],
                node_size=700, font_size=10)
        plt.title(f"Bipartite Realization of sequence {test_sequence}")
        plt.show()
    else:
        print("The sequence is not bipartite-realizable.")
