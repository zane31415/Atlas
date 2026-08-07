# Atlas — exact minimal threshold circuits for all 4-input Boolean functions

[![verify](https://github.com/zane31415/Atlas/actions/workflows/verify.yml/badge.svg)](https://github.com/zane31415/Atlas/actions/workflows/verify.yml)

Exact minimum-cost realizations of **every 4-input Boolean function**
(65,536 truth tables; 222 NPN equivalence classes) as feedforward
integer-weight **threshold circuits**, with machine-checkable verification
of every stored circuit and explicit flags on every value whose optimality
proof is incomplete.

The cost model is **evaluation cost**: `cost = wires + gates`, where wires
= nonzero weights (the interconnect currency of Kane–Williams, STOC 2016)
and gates = threshold nodes.

## Two circuit models, both tabulated

Whether a gate may read the **raw inputs** as well as the layer below it
turns out to change every number in this atlas, so both models are stored
and the difference is a first-class object.

| model | who can read the inputs | file |
|---|---|---|
| **skip** (default) | any gate — every earlier layer *and* the raw inputs, the standard model in circuit complexity and the one Kane–Williams use | `data/n4_skip.jsonl` |
| **layered** | only the first layer | `data/n4_atlas.jsonl` |

Measured exhaustively at n=4 at the same weight bound: **207 of 222 classes
are strictly cheaper in the skip model** — median 2, mean 2.88, maximum 9.
The 15 unaffected classes are the depth-1 ones, where the two models
coincide by definition. **Parity-4 costs 25 layered and 16 skip, both
CP-SAT–proven.**

The skip search is not over a fixed list of architectures. It runs over
every shape, of any depth and any width, that a support-and-wiring lower
bound cannot retire against the incumbent — a finite family, derived per
class rather than chosen. Parity-4 was searched to depth 6.

The restriction is *not* a harmless rescaling: **4.8% of class pairs swap
cost order** between the models, so comparative claims do not transfer.

The clearest single case is XOR of two variables (`0x0ff0`). Layered, it
needs 3 gates and 6 wires; with one skip wire it needs 2 gates and 5:

```
h0_0 = [ -1*x2 + 1*x3 + (0) >= 0 ]          cost 7  (skip)
out  = [ -2*h0_0 + -1*x2 + 1*x3 + (1) >= 0 ]   vs   cost 9 (layered)
```

**Why the skip table is the default.** It is the standard model, it is the
one the cost model's own citation uses, and the structural fact below is
invisible without it: 174 of the 222 optima use an architecture that
*cannot exist* in the layered model.

### What is stored in each

The **skip** table has one regime — free weights — with the proven optimal
circuit, its architecture, and its layered counterpart for comparison.

The **layered** table is the older and deeper one, with three weight
regimes:

| regime | what is stored |
|---|---|
| **free** (unbounded integer weights) | the (1,1)-optimal circuit per class: full weights, biases, architecture |
| **\|w\| ≤ 2** and **\|w\| ≤ 1** | the complete **Pareto frontier** over (gates, wires) — every non-dominated point with a verified circuit, so the optimum under *any* cost `a·wires + b·gates` (a,b > 0) is a lookup, no solver needed |

Plus two derived tables, both layered-model: **constructive optima** (85
classes where a human-readable construction — a single gate, or a
shell/decision-list circuit — provably matches the exact optimum) and the
**price of decomposability** (for the 48 disjoint-decomposable classes, the
exact minimum cost over circuits that respect the decomposition, versus the
unrestricted optimum).

So the capped regimes, the constructive forms and the fold prices are
currently available **only** in the layered model. Re-deriving them under
skips is open work, not an oversight.

## Quick start (no dependencies — standard library only)

Look up the minimal circuit for any truth table:

```
$ python tools/atlas_lookup.py 0x6996
truth table 0x6996  (NPN class 0x6996, ...)
source: exact solver (skip model, free weights, W=7); layered cost for the
        same class is 25, so the layering restriction costs 9
cost 16 = 13 wires + 3 gates   optimality proven: True
  h0_0 = [ 6*x1 + -1*x2 + (-6) >= 0 ]
  h1_0 = [ -3*h0_0 + -2*x0 + 1*x1 + -2*x2 + 2*x3 + (1) >= 0 ]
  out  = [ -6*h0_0 + -6*h1_0 + -5*x0 + 3*x1 + -3*x2 + 3*x3 + (5) >= 0 ]
verified: circuit reproduces the requested table on all 16 inputs
```

`0x6996` is 4-bit parity; note both upper gates reading the raw inputs
directly, and `out` reading **both** hidden gates — the skip to `h0_0`
across `h1_0` is what makes this 16 rather than 17.

Add `--model layered` for the strict-layered optimum, and then
`--constructive` for the readable form (popcount shells), `--regime w2|w1`
for capped weights, `--metric node_primary|wire_primary|wire10|gate10` for
other cost ratios. `--json` gives machine output. The tool maps your table
to its NPN class, transforms the stored circuit back, and **re-verifies on
all 16 inputs before printing** — you never have to trust the transform.

Re-verify the entire dataset from scratch:

```
$ python tools/verify_atlas.py
n4_atlas.jsonl: 222 classes; 221 free circuits and 484 frontier points verified; ...
ALL CHECKS PASSED
```

## Reference facts readable off the tables

### Skip model

- **174 of 222 optima use `arch=[1]`** — a single hidden gate feeding an
  output that also reads the inputs. That architecture *degenerates* in the
  layered model (an output gate reading one bit can only compute that bit,
  its negation, or a constant), which is why layered gate counts jump
  0, 1, 3, 4, 5 with **no circuit ever using exactly 2 gates**. The gap is
  not a curiosity about a missing rung: it is where 78% of n=4 optima live
  once the model permits them. Skip gate counts run 0, 1, 2, 3 with no gap,
  and the maximum cost falls from 25 to 16.
- **Depth 3 pays at n=4 — but only in this model.** Ten classes have a
  depth-3 optimum (`arch=[1,1]`: two hidden layers of one gate each), all
  proven, parity-4 among them. In the layered model the count is **0 of
  222**: depth 3 never pays there at free weights. Every one of the ten
  uses the shape layering cannot express, so this is the same `arch=[1]`
  degeneracy one layer deeper. A depth claim is only ever a claim about a
  model.
- **Depth 4 never pays at n=4**, in either model — proven, not assumed. The
  search reaches depth 6 where the bound allows it, and no class improves
  past depth 3.

| class | skip cost | arch | layered cost |
|---|---|---|---|
| `0x6996` (parity-4) | 16 | `[1,1]` | 25 |
| `0x1698`, `0x169a`, `0x19e3` | 14 | `[1,1]` | 17 |
| `0x1be4` | 14 | `[1,1]` | 18 |
| `0x17ac`, `0x179a`, `0x07b4` | 13 | `[1,1]` | 16 |
| `0x0672`, `0x0776` | 13 | `[1,1]` | 15 |

- **The layering tax is architecture-quantized**, not smooth. Its
  distribution {2: 130, 3: 28, 4: 5, 5: 37, 6: 4, 7: 1, 8: 1, 9: 1} has a
  secondary mode at 5 which is a single collapse: 37 classes whose layered
  optimum needs three hidden gates but whose skip optimum needs one, each
  dropping exactly 2 gates and 3 wires.
- Every one of the 207 taxed classes saves at least one **gate**; the bill
  splits 43% gates / 57% wires.

### Layered model

- **1,882** of 65,536 functions are threshold functions (single gate);
  their minimum-wire gates are stored.
- **Weights never need to exceed 3.** Every class attains its free-regime
  optimal cost with `|w| ≤ 3`, and the stored circuit is the
  magnitude-minimal one: 1 constant, 57 classes ternary (`|w| ≤ 1`), 120 at
  2, 44 at 3. This is a property of the *representative*, not of the cost —
  the (1,1) metric never charged for magnitude, so the circuits stored
  before 2026-08-06 were raw solver output sitting on the `W = 7` search
  bound (148 of 222 did). Same optima, same costs, same architectures,
  smaller numbers. The skip table has had no such sweep, so its `mw` column
  is still solver slack.
- **Depth never pays at free weights at n=4** *in this model*: every
  free-regime optimum is depth ≤ 2, 0 of 222. (Contrast the skip model
  above, where 10 classes have a proven depth-3 optimum — the statement is
  model-dependent, not a fact about threshold circuits.) Depth-3
  realizations become competitive here only under weight caps: at |w| ≤ 2
  they appear only as wire-saving frontier points, but at |w| ≤ 1 a depth-3
  circuit is the strict `wires+gates` optimum for 13 of the 222 classes
  (all proven; for three of them no depth-2 point makes the frontier at
  all).
- **33 of 222 classes have a genuine gate↔wire tradeoff under weight caps**
  (a multi-point Pareto frontier). 30 of the 33 are binate in all four
  variables; XOR-decomposable classes hit the tradeoff at 54% vs a 15%
  base rate. The lone fully-unate exception (0x011f) is a threshold
  function that needs a weight of magnitude 2 — cap starvation, not
  structure.
- **Decomposable ≠ decomposed.** 46 of the 48 disjoint-decomposable
  classes pay a strictly positive premium for *any* circuit that respects
  their block structure (modal premium: +2, on 25 of 48 classes; category
  means 2.3–3.3; two classes fold at no cost; all 48 records proven). For XOR- and MUX-type decompositions the
  premium provably includes a depth-3 requirement: a depth-2 circuit whose
  first layer is confined to the blocks computes `[F_A(x_A) + F_B(x_B) ≥ 0]`,
  and a 2×2 exchange argument shows no such split realizes an XOR/MUX
  block table at any width. The per-bipartition feasibility verdicts (a
  small LP over block potentials) are stored with each record.
- Constructive (readable) forms achieve the exact optimum for **85/222
  classes**: all 15 threshold classes, all 5 symmetric classes, and 65
  others via shell/decision-list circuits. For the remaining 137 classes
  the stored optima are solver witnesses with no known readable form at
  equal cost.

## Honesty flags — read before citing numbers

- **Optimality proofs.** Every circuit is *verified* (it computes its
  table — you can re-check this yourself, see above). Optimality is
  *proven* per point by an exact CP-SAT solve except where flagged
  `proven: false`:
  - `n4_atlas.jsonl`: **3 frontier points**, all on class 0x6996 (parity),
    each a verified upper bound carrying a wire lower bound (`wlb`).
  - `n4_fold_price.jsonl`: **none** — as of v1.0.1 all 48 records are
    proven. (v1.0.0 flagged 12; the certification recheck closed all 12
    without changing any premium — every flagged upper bound was already
    tight.)
- **Witness non-uniqueness.** A stored circuit is *one* minimum-cost
  circuit; minimum-cost circuits are generally not unique, and
  structural statements about "the" optimum should be phrased as
  statements about the stored witness unless the table proves otherwise
  (the fold-price table is such a proof: a positive proven premium
  certifies that *no* block-respecting circuit matches the optimum).
- **Scope.** n = 4, exact and exhaustive. Nothing here is an asymptotic
  claim; several of the regularities above are known to be
  small-n-specific (e.g., unbounded weights are worth roughly one layer of
  depth in general, but which functions exploit the exchange changes
  with n).

## Files

```
data/n4_skip.jsonl                SKIP model: per-class minimal circuits (free
                                  weights), with the layered cost and the tax
data/n4_atlas.jsonl               LAYERED model: per-class minimal circuits +
                                  capped Pareto frontiers
data/n4_constructive_optima.jsonl readable circuits matching the exact optimum (85)
data/n4_fold_price.jsonl          price of decomposability per decomposable class (48)
data/n4_categories.jsonl          structural category per class
data/SCHEMAS.md                   precise schemas and encoding conventions
n4_summary.csv                    one-row-per-class browsable summary, both
                                  models side by side (layered_* and skip_*)
tools/atlas_lookup.py             truth table -> verified minimal circuit (stdlib)
tools/verify_atlas.py             re-verify every stored circuit (stdlib)
mm_oracle.py                      the exact-synthesis library used to build the
                                  tables (CP-SAT / LP; needs numpy, scipy, ortools)
```

Model conventions (bit order, circuit encoding, NPN group) are specified at
the top of [data/SCHEMAS.md](data/SCHEMAS.md).

## Reproducing / extending

`mm_oracle.py` is the solver library the tables were built with: exact
minimum-wire realization per architecture (CP-SAT), feasibility at a cost
budget, circuit verification, and NPN canonicalization. `python
mm_oracle.py` runs its self-checks. Dependencies for solving only:
`numpy`, `scipy`, `ortools` (the lookup and verification tools need
nothing).

## Relation to prior work

- **Muroga's threshold-function enumerations** (S. Muroga, *Threshold Logic
  and its Applications*, Wiley, 1971) are the classical tables for single
  threshold gates. Our single-gate census reproduces the classical count —
  1,882 of the 65,536 4-input functions are threshold functions — as an
  independent cross-check; the atlas extends the enumeration from single
  gates to minimum-cost multi-gate circuits.
- **The cost model** (wires = nonzero weights as the complexity currency)
  follows D. M. Kane and R. Williams, "Super-linear gate and
  super-quadratic wire lower bounds for depth-two and depth-three threshold
  circuits," *STOC 2016* (arXiv:1511.07860). Their LTF∘LTF gates may read
  input variables as well as previous gate outputs — i.e. their model
  *admits* skips — so `data/n4_skip.jsonl` is the table comparable to their
  quantity, and the layered values in `data/n4_atlas.jsonl` are a strictly
  larger one. Earlier revisions of this README and of `mm_oracle.py`
  attributed the cost model to them without that qualification; that
  attribution was wrong and is corrected here.
- **The layered restriction itself** is a studied class: A. Gál and
  J.-T. K. Jang, "The size and depth of layered Boolean circuits,"
  *Information Processing Letters* 111(5):213–217, 2011. In the neural-net
  setting the same asymmetry appears as H. Lin and S. Jegelka, "ResNet with
  one-neuron hidden layers is a universal approximator," *NeurIPS 2018*
  (arXiv:1806.10909): skip connections make width-1 layers universal while
  plain narrow nets are not.
- **The nearest relative in spirit** is Knuth's exhaustive small-n
  optimal-circuit computation for 4- and 5-variable functions over
  two-input Boolean gates (*TAOCP* Vol. 4A, §7.1.2, "Boolean evaluation").
  This atlas plays the same role for integer-weight threshold gates under a
  wire+gate cost, adds complete Pareto frontiers under weight caps, and
  ships per-point optimality certificates/flags and re-verification tools.

Exact synthesis of threshold networks is an active EDA topic (e.g.
SAT/CP-based exact synthesis, threshold-logic decomposition); this
repository is a reference *table with certificates* for the complete
4-input space, not a synthesis tool for larger n.

## License

Code (tools, solver library): **MIT** — see [LICENSE](LICENSE).
Data (`data/*.jsonl`, `n4_summary.csv`): **CC0 1.0** (public domain
dedication) — see [data/LICENSE](data/LICENSE).

## Citing

If you use these tables or tools, please cite this repository (tagged
releases are immutable reference points):

```
Atlas: exact minimal threshold circuits for all 4-input Boolean functions.
https://github.com/zane31415/Atlas, v1.0.1, 2026.
```
