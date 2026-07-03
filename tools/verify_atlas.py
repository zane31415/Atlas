#!/usr/bin/env python3
"""Re-verify every circuit in the atlas data files. Pure standard library.

Checks, per data file:
  n4_atlas.jsonl              every free circuit and every capped frontier
                              point computes its class's truth table; capped
                              points respect their weight bound; stored
                              gate/wire counts match the circuit.
  n4_constructive_optima.jsonl  every circuit computes its table AND its cost
                              equals the free optimum stored in n4_atlas.
  n4_fold_price.jsonl         every fold circuit computes its table, its
                              layer-1 gates respect the stated bipartition,
                              and premium = fold_cost - free_opt.

Exit code 0 iff everything passes. Run from the repository root:
  python tools/verify_atlas.py
"""
import json
import sys
from pathlib import Path

N = 4
DATA = Path(__file__).resolve().parent.parent / "data"


def eval_circuit(ckt, x):
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
    wires = gates = maxw = 0
    for layer in ckt:
        for g in layer:
            nz = sum(1 for w in g[:-1] if w != 0)
            if nz:
                wires += nz
                gates += 1
                maxw = max(maxw, max(abs(w) for w in g[:-1]))
    return wires, gates, maxw


def check(cond, msg, fails):
    if not cond:
        fails.append(msg)


def main():
    fails = []
    atlas = {}
    with open(DATA / "n4_atlas.jsonl") as f:
        for line in f:
            r = json.loads(line)
            atlas[r["canon"]] = r

    n_free = n_pts = unproven = 0
    for T, r in atlas.items():
        fr = r["regimes"]["free"]["balanced_11"]
        if fr["ckt"]:
            check(table_of(fr["ckt"]) == T, f"atlas free ckt wrong: 0x{T:04x}", fails)
            w, g, _ = cost_of(fr["ckt"])
            check(w + g == fr["cost"], f"atlas free cost mismatch: 0x{T:04x}", fails)
            n_free += 1
        for regime, bound in (("w2", 2), ("w1", 1)):
            for i, pt in enumerate(r["regimes"][regime]["frontier"]):
                if not pt["ckt"]:
                    continue
                U = table_of(pt["ckt"])
                # constant classes store a bias-only gate; table may be the
                # complement realization of the same class member
                check(U == T or (T in (0,) and U in (0, 0xFFFF)),
                      f"atlas {regime} pt wrong: 0x{T:04x}[{i}]", fails)
                w, g, mw = cost_of(pt["ckt"])
                check(mw <= bound, f"atlas {regime} weight bound: 0x{T:04x}[{i}]", fails)
                check(w == pt["w"] and g == pt["g"],
                      f"atlas {regime} g/w mismatch: 0x{T:04x}[{i}]", fails)
                n_pts += 1
                if not pt["proven"]:
                    unproven += 1
    print(f"n4_atlas.jsonl: {len(atlas)} classes; {n_free} free circuits and "
          f"{n_pts} frontier points verified; {unproven} frontier points "
          f"flagged proven=false (documented)")

    n_c = 0
    with open(DATA / "n4_constructive_optima.jsonl") as f:
        for line in f:
            r = json.loads(line)
            if not r["ckt"]:
                continue
            T = r["canon"]
            check(table_of(r["ckt"]) == T, f"constructive wrong: 0x{T:04x}", fails)
            w, g, _ = cost_of(r["ckt"])
            free = atlas[T]["regimes"]["free"]["balanced_11"]["cost"]
            check(w + g == r["cost"] == free,
                  f"constructive not at optimum: 0x{T:04x}", fails)
            n_c += 1
    print(f"n4_constructive_optima.jsonl: {n_c} circuits verified at the free optimum")

    n_f = fold_unproven = 0
    with open(DATA / "n4_fold_price.jsonl") as f:
        for line in f:
            r = json.loads(line)
            T = r["canon"]
            check(table_of(r["ckt"]) == T, f"fold ckt wrong: 0x{T:04x}", fails)
            w, g, _ = cost_of(r["ckt"])
            check(w + g == r["fold_cost"], f"fold cost mismatch: 0x{T:04x}", fails)
            check(r["premium"] == r["fold_cost"] - r["free_opt"],
                  f"fold premium arithmetic: 0x{T:04x}", fails)
            A, B = (set(int(v) for v in part.strip("()").split(",") if v.strip())
                    for part in r["bipartition"].split("|"))
            for gate in r["ckt"][0]:
                sup = {j for j in range(N) if gate[j] != 0}
                check(sup <= A or sup <= B,
                      f"fold mask violated: 0x{T:04x}", fails)
            n_f += 1
            if not r["proven"]:
                fold_unproven += 1
    print(f"n4_fold_price.jsonl: {n_f} fold circuits verified "
          f"({fold_unproven} flagged proven=false: premiums are upper bounds)")

    if fails:
        print(f"\nFAILURES ({len(fails)}):")
        for m in fails:
            print(" ", m)
        sys.exit(1)
    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()
