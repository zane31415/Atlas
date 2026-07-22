# Atlas — exact minimal threshold circuits for all 4-input Boolean functions

[![verify](https://github.com/zane31415/Atlas/actions/workflows/verify.yml/badge.svg)](https://github.com/zane31415/Atlas/actions/workflows/verify.yml)

Exact minimum-cost realizations of **every 4-input Boolean function**
(65,536 truth tables; 222 NPN equivalence classes) as feedforward
integer-weight **threshold circuits**, with machine-checkable verification
of every stored circuit and explicit flags on every value whose optimality
proof is incomplete.

The cost model is **evaluation cost**: `cost = wires + gates`, where wires
= nonzero weights (the interconnect currency of Kane–Williams, STOC 2016 —
but see the model note below) and gates = threshold nodes. Three weight
regimes are tabulated:

> ### ⚠️ Model scope: every cost here is a **strict-layered-model** cost
>
> In this atlas each gate reads **only the previous layer**, so the output
> gate cannot read the raw inputs — no skip connections. The standard model
> in circuit complexity *does* allow them, and the restriction is not free.
> Measured exhaustively at n=4: **all 207 non-LTF classes get cheaper once
> skips are allowed**, by a median of 2 and a mean of 2.86, and **parity-4
> drops from 25 to 17** (both values CP-SAT–proven). Only the 15 depth-1
> classes are unaffected, trivially. The tax is *non-uniform*: 4.8% of class
> pairs swap cost order, so it is not a harmless rescaling.
>
> This also qualifies the Kane–Williams reference above. Their LTF∘LTF model
> lets the output gate read input variables as well as previous gate outputs
> — i.e. it *admits* skips — so the costs tabulated here are a strictly
> larger quantity than the one they bound. Earlier revisions of this README
> and of `mm_oracle.py` attributed the cost model to them without that
> qualification; that attribution was wrong and is corrected here.
>
> Everything below is correct and proven **within the layered model**. Read
> "minimum cost" as "minimum cost among strict-layered circuits" throughout.
> The skip-model measurements are not yet published in this repository; this
> repository ships proven values only, and that table will land when it is
> packaged to the same standard as the data here.

| regime | what is stored |
|---|---|
| **free** (unbounded integer weights) | the (1,1)-optimal circuit per class: full weights, biases, architecture |
| **\|w\| ≤ 2** and **\|w\| ≤ 1** | the complete **Pareto frontier** over (gates, wires) — every non-dominated point with a verified circuit, so the optimum under *any* cost `a·wires + b·gates` (a,b > 0) is a lookup, no solver needed |

Plus two derived tables: **constructive optima** (85 classes where a
human-readable construction — a single gate, or a shell/decision-list
circuit — provably matches the exact optimum) and the **price of
decomposability** (for the 48 disjoint-decomposable classes, the exact
minimum cost over circuits that respect the decomposition, versus the
unrestricted optimum).

## Quick start (no dependencies — standard library only)

Look up the minimal circuit for any truth table:

```
$ python tools/atlas_lookup.py 0x6996
truth table 0x6996  (NPN class 0x6996, ...)
cost 25 = 20 wires + 5 gates   optimality proven: True
  h0_0 = [ -2*x0 + -2*x1 + 1*x2 + 2*x3 + (1) >= 0 ]
  ...
verified: circuit reproduces the requested table on all 16 inputs
```

`0x6996` is 4-bit parity. Add `--constructive` for the readable form
(popcount shells), `--regime w2|w1` for capped weights, `--metric
node_primary|wire_primary|wire10|gate10` for other cost ratios, `--json`
for machine output. The tool maps your table to its NPN class, transforms
the stored circuit back, and **re-verifies on all 16 inputs before
printing** — you never have to trust the transform.

Re-verify the entire dataset from scratch:

```
$ python tools/verify_atlas.py
n4_atlas.jsonl: 222 classes; 221 free circuits and 484 frontier points verified; ...
ALL CHECKS PASSED
```

## Reference facts readable off the tables

- **1,882** of 65,536 functions are threshold functions (single gate);
  their minimum-wire gates are stored.
- **Depth never pays at free weights at n=4**: every free-regime optimum
  is depth ≤ 2. Depth-3 realizations only become competitive under weight
  caps: at |w| ≤ 2 they appear only as wire-saving frontier points, but at
  |w| ≤ 1 a depth-3 circuit is the strict `wires+gates` optimum for 13 of
  the 222 classes (all proven; for three of them no depth-2 point makes
  the frontier at all).
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
data/n4_atlas.jsonl               per-class minimal circuits + capped Pareto frontiers
data/n4_constructive_optima.jsonl readable circuits matching the exact optimum (85)
data/n4_fold_price.jsonl          price of decomposability per decomposable class (48)
data/n4_categories.jsonl          structural category per class
data/SCHEMAS.md                   precise schemas and encoding conventions
n4_summary.csv                    one-row-per-class browsable summary
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
  circuits," *STOC 2016* (arXiv:1511.07860) — in the *currency*, not in the
  *circuit model*: their LTF∘LTF gates may read input variables as well as
  previous gate outputs, whereas this atlas is strictly layered (see the
  model note at the top). The values here are therefore an upper bound on
  the Kane–Williams quantity, not equal to it.
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
