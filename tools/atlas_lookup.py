#!/usr/bin/env python3
"""Look up the exact minimal threshold circuit for ANY 4-input Boolean
function (all 65,536 truth tables), from the n=4 atlas.

Pure standard library — no solver required. The atlas stores one verified
minimal circuit per NPN equivalence class; this tool maps your truth table
to its class, transforms the stored circuit back through the NPN
transformation, and re-verifies the result against your table before
printing it.

Usage:
  python tools/atlas_lookup.py 0x6996                 # free-weight (1,1) optimum
  python tools/atlas_lookup.py 0x6996 --regime w2     # |w|<=2 Pareto frontier
  python tools/atlas_lookup.py 0x6996 --regime w1 --metric wire_primary
  python tools/atlas_lookup.py 0x1ee --constructive   # readable form, if one exists

Truth-table convention: bit m of the 16-bit table is f(x) for the
assignment with x_j = (m >> j) & 1, j = 0..3.

Circuit convention: a circuit is a list of layers; each layer is a list of
gates [w_0, ..., w_{k-1}, bias]; a gate fires (outputs 1) iff
sum(w_i * input_i) + bias >= 0. Layer-1 gates read the 4 inputs; each later
layer reads the previous layer's gate outputs; the last layer is a single
output gate. Cost = wires (nonzero weights) + gates (with >=1 wire).
"""
import argparse
import json
import sys
from itertools import permutations
from pathlib import Path

N = 4
DATA = Path(__file__).resolve().parent.parent / "data"


# ---------------------------------------------------------------- evaluation

def eval_circuit(ckt, x):
    """x = list of 4 bits; returns the output bit."""
    prev = x
    for layer in ckt:
        prev = [1 if sum(w * v for w, v in zip(g[:-1], prev)) + g[-1] >= 0 else 0
                for g in layer]
    return prev[0]


def table_of(ckt):
    U = 0
    for m in range(16):
        if eval_circuit(ckt, [(m >> j) & 1 for j in range(N)]):
            U |= 1 << m
    return U


def cost_of(ckt):
    wires = gates = 0
    for layer in ckt:
        for g in layer:
            nz = sum(1 for w in g[:-1] if w != 0)
            if nz:
                wires += nz
                gates += 1
    return dict(cost=wires + gates, wires=wires, gates=gates)


# ---------------------------------------------------- NPN group and transforms

def transform_table(T, perm, neg, flip):
    """g(x) = f(y) ^ flip with y_j = x_{perm[j]} ^ neg_j."""
    U = 0
    for m in range(16):
        y = 0
        for j in range(N):
            bit = ((m >> perm[j]) & 1) ^ ((neg >> j) & 1)
            y |= bit << j
        v = ((T >> y) & 1) ^ flip
        if v:
            U |= 1 << m
    return U


def npn_canon(T):
    return min(transform_table(T, p, s, f)
               for p in permutations(range(N)) for s in range(16) for f in (0, 1))


def transform_circuit(ckt, perm, neg, flip):
    """Apply an input permutation/negation to layer 1 and an output negation
    to the final gate. Correctness is not derived from convention: the
    caller tries group elements until the transformed circuit VERIFIES."""
    l1 = []
    for g in ckt[0]:
        w, b = list(g[:-1]), g[-1]
        nw = [0] * N
        for j in range(N):
            wj = w[j]
            if (neg >> j) & 1:
                b += wj
                wj = -wj
            nw[perm[j]] = wj
        l1.append(nw + [b])
    new = [l1] + [[list(g) for g in layer] for layer in ckt[1:]]
    if flip:
        g = new[-1][0]
        new[-1][0] = [-w for w in g[:-1]] + [-g[-1] - 1]
    return new


def realize(ckt, T):
    """Return a transform of the stored canonical circuit that computes T
    exactly (verified). Guaranteed to exist since T is in the class orbit."""
    for p in permutations(range(N)):
        for s in range(16):
            for f in (0, 1):
                cand = transform_circuit(ckt, p, s, f)
                if table_of(cand) == T:
                    return cand, (p, s, f)
    raise RuntimeError("no NPN transform realized the table — data corrupt?")


# ------------------------------------------------------------------- lookup

def load_atlas():
    recs = {}
    with open(DATA / "n4_atlas.jsonl") as f:
        for line in f:
            r = json.loads(line)
            recs[r["canon"]] = r
    return recs


def load_constructive():
    recs = {}
    with open(DATA / "n4_constructive_optima.jsonl") as f:
        for line in f:
            r = json.loads(line)
            recs[r["canon"]] = r
    return recs


def pick_circuit(rec, regime, metric):
    if regime == "free":
        e = rec["regimes"]["free"]["balanced_11"]
        return e["ckt"], dict(proven=e["proven"], cost=e["cost"],
                              source="exact solver (free weights, cost = wires + gates)")
    reg = rec["regimes"][regime]
    idx = reg["metrics"][metric]["idx"]
    pt = reg["frontier"][idx]
    return pt["ckt"], dict(proven=pt["proven"], cost=reg["metrics"][metric]["cost"],
                           source=f"exact solver (|w|<={regime[1]}, metric {metric}, "
                                  f"frontier point {idx})")


def fmt(ckt):
    lines = []
    for L, layer in enumerate(ckt):
        src = [f"x{j}" for j in range(N)] if L == 0 else \
              [f"h{L-1}_{j}" for j in range(len(ckt[L-1]))]
        for gi, g in enumerate(layer):
            terms = " + ".join(f"{w}*{s}" for w, s in zip(g[:-1], src) if w != 0)
            name = "out" if L == len(ckt) - 1 else f"h{L}_{gi}"
            lines.append(f"  {name} = [ {terms or '0'} + ({g[-1]}) >= 0 ]")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table", help="16-bit truth table, e.g. 0x6996 or 27030")
    ap.add_argument("--regime", choices=["free", "w2", "w1"], default="free")
    ap.add_argument("--metric", default="balanced_11",
                    choices=["balanced_11", "node_primary", "wire_primary",
                             "wire10", "gate10"])
    ap.add_argument("--constructive", action="store_true",
                    help="prefer the constructive (single-gate / shell) optimum "
                         "when one is stored (free regime only)")
    ap.add_argument("--json", action="store_true", help="print circuit as JSON only")
    args = ap.parse_args()

    T = int(args.table, 0) & 0xFFFF
    if T in (0, 0xFFFF):
        print(f"constant function {T & 1 if T else 0}"
              if T == 0 else "constant function 1")
        print("cost 0 (no gates)")
        return

    canon = npn_canon(T)
    rec = load_atlas()[canon]

    info = None
    if args.constructive and args.regime == "free":
        c = load_constructive().get(canon)
        if c is not None:
            ckt, info = c["ckt"], dict(proven=True, cost=c["cost"],
                                       source=f"constructive ({c['form']}); optimality "
                                              "inherited from the exact-solver optimum")
    if info is None:
        ckt, info = pick_circuit(rec, args.regime, args.metric)

    out, (p, s, f) = realize(ckt, T)
    assert table_of(out) == T
    c = cost_of(out)

    if args.json:
        print(json.dumps(out))
        return
    print(f"truth table 0x{T:04x}  (NPN class 0x{canon:04x}, "
          f"transform perm={p} neg=0b{s:04b} outflip={f})")
    print(f"source: {info['source']}")
    print(f"cost {c['cost']} = {c['wires']} wires + {c['gates']} gates"
          f"   optimality proven: {info['proven']}")
    print(fmt(out))
    print("verified: circuit reproduces the requested table on all 16 inputs")


if __name__ == "__main__":
    main()
