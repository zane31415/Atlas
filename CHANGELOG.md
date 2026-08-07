# Changelog

Dates are the dates of the work, not of any release.

## Unreleased (since v1.0.1)

Two changes to `data/`, one a **correction** and one an **improvement**.
Anything citing `v1.0.1` for a skip-model cost should move to this revision.

### Corrected — the skip model was too narrow at depth 3

`data/n4_skip.jsonl` was solved in a *one-step* skip model, where a gate
reads the previous layer and the raw inputs but **not** the hidden layers
before it. The standard model has no such restriction. The two agree at
depth ≤ 2 and differ at depth 3, where the one-step model cannot express

```
u1  reads u0
out reads u0 AND u1
```

for want of a weight slot. **Five of the 222 costs were too high by 1:**

| class | v1.0.1 | now |
|---|---|---|
| `0x17ac`, `0x07b4`, `0x179a` | 14 | **13** |
| `0x1be4` | 15 | **14** |
| `0x6996` (parity-4) | 17 | **16** |

The table is regenerated from a full-DAG sweep, 222/222 CP-SAT proven, over
an architecture family derived per class from a lower bound instead of a
fixed list (parity-4 searched to depth 6; depth 4 proven never to pay).
Records now carry `model: "full"`, and a gate's weight vector past layer 1
is `[every earlier layer's outputs..., x_0..x_3]`.

Knock-on numbers: maximum cost 17 → 16, depth-3 optima 7 → 10 classes,
layering tax mean 2.86 → 2.88 and maximum 8 → 9. Unchanged: 174/222 optima
at `arch=[1]`, 207 taxed classes, the 43/57 gate/wire split of the tax,
and every value in `data/n4_atlas.jsonl`'s cost columns.

### Improved — the layered atlas stores magnitude-minimal circuits

The `wires + gates` cost never charged for weight magnitude, so solver
witnesses had no reason to be small and expanded to fill the `W = 7` search
bound; 148 of 222 sat on it. `data/n4_atlas.jsonl`'s free-regime circuits
are now the magnitude-minimal representatives of the same optima —
**`max|w| ≤ 3` for every class** (57 ternary, 120 at 2, 44 at 3). NOR-4 was
`[-7,-7,-7,-7,0]`, is now `[-1,-1,-1,-1,0]`.

No cost, gate count, wire count, architecture or `proven` flag changes. Only
`ckt` and `mw`, on the 198 classes where a smaller witness exists. The
minimality bound is the atlas's own certified capped costs, so it carries
their architecture-family scope and nothing more; `tools/verify_atlas.py`
now enforces it.

`data/n4_skip.jsonl` has had **no** such sweep — its `mw` column is still
solver slack and should not be read as a precision requirement.

### Tools

- Both tools' `gate_sources()` decode three widths (layered, full-DAG,
  one-step) instead of two. Circuits published before this revision still
  read correctly.
- `tools/verify_atlas.py` additionally checks that stored `mw` matches the
  circuit and that free-regime circuits are magnitude-minimal.
- CI gains two steps that exercise a depth-3 skip circuit, including one
  reached through an NPN transform. The previous seven never touched one,
  which is why a decoder bug there could have shipped green.
- `mm_oracle.py` updated to the version the current tables were built with
  (it takes `skips=False|True|'full'`; the shipped copy previously supported
  the layered model only). `arch_family.py` added — the derived architecture
  family and the lower bound behind it.

## v1.0.1

Fold-price certification recheck: all 48 records proven (v1.0.0 flagged 12;
no premium changed). Skip-model table published for the first time.

## v1.0.0

Initial public release: layered atlas with capped Pareto frontiers,
constructive optima, fold prices, categories, lookup and verification tools.
