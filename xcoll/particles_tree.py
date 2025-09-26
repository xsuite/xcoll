# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from __future__ import annotations

import numpy as np
from typing import Iterable, Callable, Optional, Dict, List, Tuple

from xtrack.particles import LAST_INVALID_STATE


class ParticlesTree:
    def __init__(self, particle_ids, parent_ids=None, *, mask=None):
        """
        Usage:
          - ParticlesTree(particle_ids, parent_ids)
          - ParticlesTree(part)                         # uses all particles
          - ParticlesTree(part, mask=...)               # custom boolean mask
        """

        # ---- Accept either arrays or an xtrack.Particles-like object ----
        if parent_ids is None and hasattr(particle_ids, "particle_id") and hasattr(particle_ids, "parent_particle_id"):
            part = particle_ids
            if mask is None:
                if hasattr(part, "state"):
                    mask = np.asarray(part.state) > LAST_INVALID_STATE
                else:
                    mask = slice(None)  # use all rows
            ids = np.asarray(part.particle_id)[mask]
            par = np.asarray(part.parent_particle_id)[mask]
        else:
            ids = np.asarray(particle_ids)
            par = np.asarray(parent_ids)
            if par is None:
                raise ValueError("When passing arrays, both particle_ids and parent_ids are required.")

        # ---- Basic validation / normalisation ----
        if ids.ndim != 1 or par.ndim != 1 or ids.shape != par.shape:
            raise ValueError("ids and parent_ids must be 1D arrays of equal length")
        ids = ids.astype(np.int64, copy=False)
        par = par.astype(np.int64, copy=False)
        n = ids.size

        if n == 0:
            self.ids = ids
            self.parent_idx = np.zeros(0, dtype=np.int64)
            self.root_idx = np.zeros(0, dtype=np.int64)
            self.root_ids = np.zeros(0, dtype=np.int64)
            self.children_ptr = np.zeros(1, dtype=np.int64)
            self.children_idx = np.zeros(0, dtype=np.int64)
            self._order = None
            self._ids_sorted = None
            return

        # ---- id -> index (one sort) ----
        order = np.argsort(ids)
        ids_sorted = ids[order]
        parent_idx = order[np.searchsorted(ids_sorted, par)]

        # ---- pointer jumping to roots ----
        P = parent_idx.copy()
        max_iters = int(np.ceil(np.log2(max(1, n)))) + 3
        for _ in range(max_iters):
            newP = P[P]
            if np.array_equal(newP, P):
                break
            P = newP
        else:
            raise ValueError("Parent chains did not converge (cycle suspected).")

        root_idx = P
        root_ids = ids[root_idx]

        # ---- children CSR (exclude self-loops) ----
        self_idx = np.arange(n, dtype=np.int64)
        mask_child = parent_idx != self_idx
        if np.any(mask_child):
            parents = parent_idx[mask_child]
            children = np.nonzero(mask_child)[0]
            k = np.argsort(parents)
            parents = parents[k]
            children = children[k]

            counts = np.zeros(n, dtype=np.int64)
            uniq, cnts = np.unique(parents, return_counts=True)
            counts[uniq] = cnts

            ptr = np.empty(n + 1, dtype=np.int64)
            np.cumsum(counts, out=ptr[1:])
            ptr[0] = 0
            csr_idx = children
        else:
            ptr = np.zeros(n + 1, dtype=np.int64)
            csr_idx = np.zeros(0, dtype=np.int64)

        # ---- store ----
        self.ids = ids
        self.parent_idx = parent_idx
        self.root_idx = root_idx
        self.root_ids = root_ids
        self.children_ptr = ptr
        self.children_idx = csr_idx

        # lazy lookup cache (always 1-D)
        self._order = None
        self._ids_sorted = None

    # keep the rest of the class as before …
    def index_of(self, particle_id: int) -> int:
        if self._order is None:
            self._order = np.argsort(self.ids)
            self._ids_sorted = self.ids[self._order].ravel()  # ensure 1-D
        j = int(np.searchsorted(self._ids_sorted, int(particle_id)))
        if j >= self._ids_sorted.size or self._ids_sorted[j] != particle_id:
            raise KeyError(f"particle_id {particle_id} not found")
        return int(self._order[j])

    def indices_of(self, particle_ids: Iterable[int]) -> np.ndarray:
        x = np.asarray(list(particle_ids))
        if self._order is None:
            self._order = np.argsort(self.ids)
            self._ids_sorted = self.ids[self._order]
        j = np.searchsorted(self._ids_sorted, x)
        ok = (j < self._ids_sorted.size) & (self._ids_sorted[j] == x)
        if not np.all(ok):
            missing = x[~ok]
            raise KeyError(f"particle_id(s) not found: {missing.tolist()}")
        return self._order[j]

    def children_of_index(self, i: int) -> np.ndarray:
        a, b = int(self.children_ptr[i]), int(self.children_ptr[i + 1])
        return self.children_idx[a:b]

    def descendants_indices(self, i: int) -> np.ndarray:
        out = []
        stack = [int(i)]
        while stack:
            u = stack.pop()
            ks = self.children_of_index(u)
            if ks.size:
                out.extend(ks.tolist())
                stack.extend(ks.tolist())
        return np.asarray(out, dtype=np.int64)

    def descendants_ids(self, particle_id: int) -> np.ndarray:
        i = self.index_of(particle_id)
        return self.ids[self.descendants_indices(i)]

    def chain_of_id(self, particle_id: int, include_self: bool = True) -> List[int]:
        i = self.index_of(particle_id)
        chain: List[int] = []
        if include_self:
            chain.append(int(self.ids[i]))
        while True:
            p = int(self.parent_idx[i])
            if p == i:
                break
            chain.append(int(self.ids[p]))
            i = p
        return chain

    # ---------- ASCII tree (sorted) ----------
    def print_tree(self,
                   root_particle_id: Optional[int] = None,
                   label: Optional[Callable[[int], str]] = None,
                   max_depth: Optional[int] = None,
                   sort_children: bool = True,
                   sort_roots: bool = True) -> None:
        """
        Print an ASCII tree/forest. Children and roots can be sorted by particle_id.
        """
        if label is None:
            label = lambda pid: str(pid)

        def kids(i: int) -> np.ndarray:
            k = self.children_of_index(i)
            if sort_children and k.size:
                k = k[np.argsort(self.ids[k])]
            return k

        def emit(i: int, prefix: str, depth: int, is_last: bool):
            connector = "" if depth == 0 else ("└─ " if is_last else "├─ ")
            print(prefix + connector + label(int(self.ids[i])))
            if max_depth is not None and depth >= max_depth:
                return
            ks = kids(i)
            if ks.size == 0:
                return
            new_prefix = prefix + ("" if depth == 0 else ("   " if is_last else "│  "))
            for j, c in enumerate(ks):
                emit(int(c), new_prefix, depth + 1, j == ks.size - 1)

        if root_particle_id is None:
            roots = self.originals_indices()
            if sort_roots and roots.size:
                roots = roots[np.argsort(self.ids[roots])]
            for idx_r, r in enumerate(roots):
                emit(int(r), "", 0, True)
                if idx_r != len(roots) - 1:
                    print()  # blank line between trees
        else:
            emit(self.index_of(root_particle_id), "", 0, True)

    # ---------- layout (pure NumPy; sorted) ----------
    def layout_positions(self,
                         root_particle_id: Optional[int] = None,
                         max_depth: Optional[int] = None,
                         sort_children: bool = True,
                         sort_roots: bool = True,
                         x_gap: float = 1.0,
                         y_gap: float = 1.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Compute a tidy top-down layout. Returns (node_idx, x, y, depth) for the drawn nodes.
        """
        n = self.ids.size
        # roots
        if root_particle_id is None:
            roots = self.originals_indices()
            if sort_roots and roots.size:
                roots = roots[np.argsort(self.ids[roots])]
            roots = roots.tolist()
        else:
            roots = [self.index_of(int(root_particle_id))]

        def _kids(i: int) -> np.ndarray:
            k = self.children_of_index(i)
            if sort_children and k.size:
                k = k[np.argsort(self.ids[k])]
            return k

        depth = np.zeros(n, dtype=np.int64)
        seen = np.zeros(n, dtype=bool)
        post: List[int] = []

        # DFS to set depths and produce postorder
        stack: List[Tuple[int, int, int]] = []
        for r in roots:
            stack.append((int(r), 0, 0))
            while stack:
                u, st, d = stack.pop()
                if st == 0:
                    if seen[u]:
                        continue
                    seen[u] = True
                    depth[u] = d
                    stack.append((u, 1, d))
                    ks = _kids(u)
                    if max_depth is not None and d >= max_depth:
                        ks = ks[:0]
                    for c in ks[::-1]:
                        stack.append((int(c), 0, d + 1))
                else:
                    post.append(u)

        # (leaf widths computed implicitly by placing children first)
        x = np.zeros(n, dtype=float)
        y = depth.astype(float) * y_gap
        ordered_nodes: List[int] = []
        cur_x = 0.0

        def place(u: int):
            nonlocal cur_x
            ks = _kids(u)
            if max_depth is not None and depth[u] >= max_depth:
                ks = ks[:0]
            if ks.size == 0:
                x[u] = cur_x
                cur_x += x_gap
            else:
                left_edge = cur_x
                for c in ks:
                    place(int(c))
                right_edge = cur_x - x_gap
                x[u] = (left_edge + right_edge) / 2.0
            ordered_nodes.append(u)

        for r in roots:
            place(int(r))
            cur_x += x_gap  # gap between trees

        node_idx = np.array(ordered_nodes, dtype=np.int64)
        return node_idx, x[node_idx], y[node_idx], depth[node_idx]

    # ---------- plotting (Matplotlib only; sorted) ----------
    def plot_tree(self,
                  root_particle_id: Optional[int] = None,
                  max_depth: Optional[int] = None,
                  sort_children: bool = True,
                  sort_roots: bool = True,
                  with_labels: bool = True,
                  label: Optional[Callable[[int], str]] = None,
                  node_marker: str = "o",
                  node_size: float = 50.0,
                  linewidth: float = 1.0,
                  ax=None):
        """
        Draw a subtree or whole forest using Matplotlib.
        """
        import matplotlib.pyplot as plt

        if label is None:
            label = lambda pid: str(pid)

        idx, xs, ys, depth = self.layout_positions(root_particle_id=root_particle_id,
                                                   max_depth=max_depth,
                                                   sort_children=sort_children,
                                                   sort_roots=sort_roots)
        slot = {int(i): k for k, i in enumerate(idx)}

        # edges
        for i in idx:
            p = int(self.parent_idx[i])
            if p == i:
                continue
            if int(i) not in slot or int(p) not in slot:
                continue
            x0, y0 = xs[slot[p]], ys[slot[p]]
            x1, y1 = xs[slot[i]], ys[slot[i]]
            ax = ax or plt.gca()
            ax.plot([x0, x1], [y0, y1], linewidth=linewidth)

        # nodes
        ax = ax or plt.gca()
        ax.scatter(xs, ys, s=node_size, marker=node_marker)

        if with_labels:
            for i, x0, y0 in zip(idx, xs, ys):
                ax.text(x0, y0, label(int(self.ids[int(i)])), ha="center", va="bottom")

        ax.invert_yaxis()
        ax.set_aspect("equal", adjustable="box")
        ax.axis("off")
        return ax
