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

def gate_sources(L, w, acts, x):
    """The activation vector a gate's weights line up against.

    Layer 1 reads the inputs. A later gate reads EVERY earlier layer's
    outputs, in layer order, and then MAY additionally read the raw inputs
    (a skip connection), in which case its weight vector ends with N
    raw-input weights. Three widths, three models, one file format:

        len(prev)            strict layered   (n4_atlas.jsonl)
        len(all earlier) + N full DAG         (n4_skip.jsonl)
        len(prev) + N        one-step skip    (accepted, never emitted)

    The last two coincide only at L = 1, where all three models coincide
    anyway, so the width still decides without an extra tag. The one-step
    case is kept because it is a real circuit — a full-DAG circuit with
    zeros in the slots it skips over — and files published before
    2026-08-06 stored five values that were minimal only under it.
    """
    if L == 0:
        if len(w) != len(x):
            raise ValueError(f"layer-1 gate has {len(w)} weights, expected {len(x)}")
        return x
    prev, earlier = acts[-1], [v for layer in acts for v in layer]
    if len(w) == len(prev):
        return prev                                    # strict layered
    if len(w) == len(earlier) + len(x):
        return earlier + list(x)                       # full DAG
    if len(w) == len(prev) + len(x):
        return list(prev) + list(x)                    # one-step skip
    raise ValueError(f"gate has {len(w)} weights; expected {len(prev)} "
                     f"(layered), {len(earlier) + len(x)} (full DAG) or "
                     f"{len(prev) + len(x)} (one-step skip)")


def eval_circuit(ckt, x):
    """x = list of 4 bits; returns the output bit."""
    acts = []
    for L, layer in enumerate(ckt):
        acts.append([1 if sum(w * v for w, v in
                              zip(g[:-1], gate_sources(L, g[:-1], acts, x))) + g[-1] >= 0
                     else 0
                     for g in layer])
    return acts[-1][0]


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
    """Apply an input permutation/negation to every gate that reads the raw
    inputs, and an output negation to the final gate.

    In a strict-layered circuit only layer 1 reads the inputs. With skip
    connections a later gate reads them too, in the tail N slots of its
    weight vector, and those slots need exactly the same permutation and
    negation as layer 1 — transforming layer 1 alone would silently produce
    a circuit for a different function. Correctness is not derived from
    convention: the caller tries group elements until the transformed
    circuit VERIFIES against the requested table.

    Only the tail matters here, and the head is carried through untouched,
    so the layer-source convention (previous layer only, or every earlier
    layer) does not affect this function — but the tail must be located by
    the same widths gate_sources uses.
    """
    new, sizes = [], []
    for L, layer in enumerate(ckt):
        out = []
        nprev = sizes[-1] if sizes else 0
        nall = sum(sizes)
        for g in layer:
            w, b = list(g[:-1]), g[-1]
            if L == 0:
                head, tail = [], w                    # all slots are inputs
            elif len(w) == nall + N:
                head, tail = w[:nall], w[nall:]       # full DAG: tail is inputs
            elif len(w) == nprev + N:
                head, tail = w[:nprev], w[nprev:]     # one-step: tail is inputs
            else:
                head, tail = w, []                    # layered: no input slots
            if tail:
                nt = [0] * N
                for j in range(N):
                    wj = tail[j]
                    if (neg >> j) & 1:
                        b += wj
                        wj = -wj
                    nt[perm[j]] = wj
            else:
                nt = []
            out.append(list(head) + nt + [b])
        new.append(out)
        sizes.append(len(layer))
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


def load_skip():
    recs = {}
    with open(DATA / "n4_skip.jsonl") as f:
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
        if L == 0:
            src = [f"x{j}" for j in range(N)]
        else:
            prev = [f"h{L-1}_{j}" for j in range(len(ckt[L - 1]))]
            earlier = [f"h{k}_{j}" for k in range(L) for j in range(len(ckt[k]))]
            ins = [f"x{j}" for j in range(N)]
            # the width says which sources the slots line up against; same
            # three cases as gate_sources, and the same order
            width = len(layer[0]) - 1 if layer else len(prev)
            if width == len(earlier) + N:
                src = earlier + ins            # full DAG
            elif width == len(prev) + N:
                src = prev + ins               # one-step skip
            else:
                src = prev                     # strict layered
        for gi, g in enumerate(layer):
            terms = " + ".join(f"{w}*{s}" for w, s in zip(g[:-1], src) if w != 0)
            name = "out" if L == len(ckt) - 1 else f"h{L}_{gi}"
            lines.append(f"  {name} = [ {terms or '0'} + ({g[-1]}) >= 0 ]")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table", help="16-bit truth table, e.g. 0x6996 or 27030")
    ap.add_argument("--model", choices=["skip", "layered"], default="skip",
                    help="skip (default) = the standard circuit model, where a "
                         "gate may also read the raw inputs; layered = every "
                         "gate reads only the previous layer. Free weights "
                         "only for skip; the capped regimes are layered-only.")
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

    if args.model == "skip":
        if args.regime != "free":
            sys.exit("--model skip has free weights only; the |w|<=2 and |w|<=1 "
                     "Pareto frontiers are tabulated in the layered model. Use "
                     "--model layered with --regime, or drop --regime.")
        if args.constructive:
            sys.exit("--constructive is a layered-model table; use --model layered.")
        s = load_skip()[canon]
        ckt = s["ckt"]
        info = dict(proven=s["proven"], cost=s["cost"],
                    source=f"exact solver (skip model, free weights, W={s['W']}); "
                           f"layered cost for the same class is {s['layered_cost']}"
                           + (f", so the layering restriction costs {s['tax']}"
                              if s["tax"] else " (unaffected: depth 1)"))
    else:
        info = None
    if info is None and args.constructive and args.regime == "free":
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
