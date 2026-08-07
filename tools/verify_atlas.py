#!/usr/bin/env python3
"""Re-verify every circuit in the atlas data files. Pure standard library.

Checks, per data file:
  n4_atlas.jsonl              every free circuit and every capped frontier
                              point computes its class's truth table; capped
                              points respect their weight bound; stored
                              gate/wire counts match the circuit; and the
                              stored free circuit is the MAGNITUDE-MINIMAL
                              optimum — max|w| <= 3, and equal to the k*
                              implied by that class's own w1/w2 costs, so a
                              regenerated file cannot quietly revert to raw
                              solver output sitting on the search bound.
  n4_constructive_optima.jsonl  every circuit computes its table AND its cost
                              equals the free optimum stored in n4_atlas.
  n4_fold_price.jsonl         every fold circuit computes its table, its
                              layer-1 gates respect the stated bipartition,
                              and premium = fold_cost - free_opt.
  n4_skip.jsonl               every skip-model circuit computes its table and
                              respects its weight bound; its layered_cost
                              agrees with n4_atlas; and tax = layered - skip
                              is arithmetic AND non-negative (the layered
                              model is the skip model with the input wires
                              forbidden, so it can never be cheaper).

Exit code 0 iff everything passes. Run from the repository root:
  python tools/verify_atlas.py
"""
import json
import sys
from pathlib import Path

N = 4
DATA = Path(__file__).resolve().parent.parent / "data"


def gate_sources(L, w, acts, x):
    """Layer 1 reads the inputs. A later gate reads every earlier layer in
    order, and with a skip connection also the raw inputs in the tail N
    slots. Three widths, three models — len(prev) is strict layered,
    len(all earlier)+N is the full DAG, len(prev)+N is a one-step skip —
    and they are distinct except at L=1, where the models coincide. One
    evaluator covers every data file."""
    if L == 0:
        if len(w) != len(x):
            raise ValueError(f"layer-1 gate has {len(w)} weights, expected {len(x)}")
        return x
    prev, earlier = acts[-1], [v for layer in acts for v in layer]
    if len(w) == len(prev):
        return prev
    if len(w) == len(earlier) + len(x):
        return earlier + list(x)
    if len(w) == len(prev) + len(x):
        return list(prev) + list(x)
    raise ValueError(f"gate has {len(w)} weights; expected {len(prev)} "
                     f"(layered), {len(earlier) + len(x)} (full DAG) or "
                     f"{len(prev) + len(x)} (one-step skip)")


def eval_circuit(ckt, x):
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
            w, g, mw = cost_of(fr["ckt"])
            check(w + g == fr["cost"], f"atlas free cost mismatch: 0x{T:04x}", fails)
            check(mw == fr["mw"], f"atlas free mw mismatch: 0x{T:04x}", fails)
            # the free regime stores the MAGNITUDE-MINIMAL optimum, and no n=4
            # class needs more than 3. Raw solver output would sit at the W=7
            # search bound instead, so this is the regression guard for a
            # regenerated file quietly losing the property.
            check(mw <= 3, f"atlas free weight not magnitude-minimal: "
                           f"0x{T:04x} has max|w|={mw}", fails)
            # k* is pinned by this class's OWN certified capped costs
            w1c = r["regimes"]["w1"]["metrics"]["balanced_11"]["cost"]
            w2c = r["regimes"]["w2"]["metrics"]["balanced_11"]["cost"]
            kstar = 1 if w1c == fr["cost"] else (2 if w2c == fr["cost"] else 3)
            check(mw == kstar, f"atlas free mw {mw} != k*={kstar} implied by "
                               f"the capped costs: 0x{T:04x}", fails)
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
    flag_note = (f"({fold_unproven} flagged proven=false: premiums are upper bounds)"
                 if fold_unproven else "(all proven)")
    print(f"n4_fold_price.jsonl: {n_f} fold circuits verified {flag_note}")

    n_s = skip_unproven = taxed = 0
    with open(DATA / "n4_skip.jsonl") as f:
        for line in f:
            r = json.loads(line)
            T = r["canon"]
            check(table_of(r["ckt"]) == T, f"skip ckt wrong: 0x{T:04x}", fails)
            w, g, mw = cost_of(r["ckt"])
            check(w + g == r["cost"] and w == r["w"] and g == r["g"],
                  f"skip cost mismatch: 0x{T:04x}", fails)
            check(mw == r["mw"] and mw <= r["W"],
                  f"skip weight bound: 0x{T:04x}", fails)
            # cross-file: the layered model is the skip model with the input
            # wires forbidden, so the layered optimum can never be cheaper
            lay = atlas[T]["regimes"]["free"]["balanced_11"]["cost"]
            check(r["layered_cost"] == lay,
                  f"skip/atlas layered_cost disagree: 0x{T:04x}", fails)
            check(r["cost"] <= lay, f"skip dearer than layered: 0x{T:04x}", fails)
            check(r["tax"] == lay - r["cost"], f"skip tax arithmetic: 0x{T:04x}", fails)
            n_s += 1
            taxed += r["tax"] > 0
            if not r["proven"]:
                skip_unproven += 1
    print(f"n4_skip.jsonl: {n_s} skip circuits verified; {taxed} classes cost "
          f"strictly less than in the layered model"
          + (f"; {skip_unproven} flagged proven=false" if skip_unproven
             else "; all proven"))

    if fails:
        print(f"\nFAILURES ({len(fails)}):")
        for m in fails:
            print(" ", m)
        sys.exit(1)
    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()
