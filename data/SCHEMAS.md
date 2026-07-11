# Data schemas

## Conventions (all files)

**Truth table**: a 4-input Boolean function is a 16-bit integer `T`; bit `m`
of `T` is `f(x)` for the assignment with `x_j = (m >> j) & 1`, `j = 0..3`.

**NPN class**: functions are grouped under input permutation, input
negation, and output negation (the NPN group, 768 elements; 222 classes at
n=4). Both wire and gate counts and any `|w|` bound are invariant under this
group, so one stored circuit per class suffices; `tools/atlas_lookup.py`
maps any table to its class and transforms the stored circuit back
(re-verifying the result).

**Circuit encoding** (`ckt` fields): a list of layers; each layer a list of
gates `[w_0, ..., w_{k-1}, bias]`. A gate outputs 1 iff
`sum(w_i * input_i) + bias >= 0`. Layer 1 reads the four inputs; each later
layer reads the previous layer's gate outputs (strict layering, no skip
connections); the final layer is a single output gate. Empty `ckt` = a
constant function (no gates).

**Cost**: `cost = wires + gates`, where wires = number of nonzero weights
across all gates and gates = number of gates with at least one nonzero
input weight. Weight magnitudes do not enter the cost; they are constrained
only in the capped regimes.

**`proven`**: `true` = the exact solver certified optimality (or, in
`n4_constructive_optima`, equality with a certified optimum); `false` = the
value is a verified upper bound whose optimality proof timed out (each such
point carries a wire lower bound `wlb` where applicable).

## n4_atlas.jsonl — one line per NPN class (222 lines)

```
canon   : int   NPN-canonical truth table (class id)
orbit   : int   class orbit size (number of distinct tables in the class)
regimes :
  free  : unbounded integer weights
    balanced_11 : {g, w, mw, arch, ckt, proven, cost}
                  the (1,1)-optimal circuit: g gates, w wires, mw = max |weight|,
                  arch = hidden-layer sizes, cost = w + g
  w2, w1 : |w| <= 2 and |w| <= 1 (bounds on weights AND biases)
    frontier : [{g, w, mw, arch, ckt, proven, wlb}, ...]
               the full Pareto frontier over (gates, wires) — every
               non-dominated point, each with a verified circuit. The
               optimum under ANY cost a*wires + b*gates with a,b > 0 is a
               frontier point, so all such metrics are lookups.
    metrics  : {balanced_11, node_primary, wire_primary, wire10, gate10}
               each -> {idx, cost}: the frontier index optimal under that
               metric (balanced = w+g; node_primary = min gates, wires
               tiebreak; wire_primary = min wires, gates tiebreak;
               wire10/gate10 = 10*w+g / w+10*g)
    feasible, timeouts : bookkeeping
```

## n4_constructive_optima.jsonl — one line per class with a constructive optimum (85 lines)

Classes where a solver-free construction provably achieves the exact free
(1,1) optimum. `form` = `single-gate` (the function is a threshold
function; one minimum-wire gate), `shell` (a decision-list / shell circuit:
one hidden layer of chunk detectors plus a priority-weighted output
threshold), or `constant`. `shells` = number of shells peeled. Optimality
is inherited by cost-equality with `n4_atlas`'s certified optimum.

```
canon, category, cost, form, shells, ckt
```

## n4_fold_price.jsonl — one line per disjoint-decomposable class (48 lines)

The **price of decomposability**: exact minimum (1,1) cost over
*fold-respecting* circuits — circuits whose layer-1 gates each read
variables from only one block of some bipartition of the function's support
(later layers unconstrained) — versus the unrestricted free optimum.

```
canon, category      : class id and its structural category
support              : variables the function depends on
free_opt             : unrestricted (1,1) optimum (from n4_atlas)
lp_depth2            : per-bipartition verdict of the depth-2 potential-split
                       LP (false = provably NO depth-2 fold-respecting
                       circuit exists at any width for that bipartition)
fold_cost, premium   : best fold-respecting cost; premium = fold_cost - free_opt
bipartition, arch,
split, depth, ckt    : the best fold circuit found (verified)
proven               : true = no cheaper fold-respecting circuit exists
                       (exhaustive within the search cap); false = premium
                       is an upper bound (solver time-limit)
inf_ub               : OPTIONAL (present on records that went through the
                       certification recheck): per-candidate map
                       "bipartition#arch#split" -> largest cost budget at
                       which that candidate was proven infeasible
                       (solver bookkeeping; consumers may ignore it)
```

As of v1.0.1 **all 48 records are proven** (the earlier 12 `proven: false`
flags were closed by a certification recheck; no premium changed).

## n4_categories.jsonl — one line per class with a cached free optimum (221 lines)

Structural category per class (`LTF`, `symmetric`, `dec:AND/OR`, `dec:XOR`,
`dec:MUX`, `prime:tangle`) plus the gap of a greedy shell-peeling
construction against the free optimum (`opt`, `peel`, `gap`, `shells`,
`frontdoor_gap`). `dec:*` = disjoint-support decomposable, split by
combining operator; `prime:tangle` = not disjointly decomposable.
