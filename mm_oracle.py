"""Exact minimum-wire oracle for integer-weight threshold circuits.

Minimum-WIRE (nonzero-weight) realization of a Boolean function as a
feedforward integer-weight THRESHOLD circuit of a given architecture,
via CP-SAT (primary) or scipy.optimize.milp/HiGHS (cross-check). Gates are
GENERAL threshold functions at every layer including the output, which is
required for the answer to be an NPN class invariant: restricting the
output gate (to AND/OR, say) gives a gate set that is not closed under
input negation, so two members of the same NPN class could then get
different costs. The threshold/bias is free; wires = total nonzero input
weights.

MODEL (load-bearing — there are THREE, pass `skips=` explicitly):

  skips=False   strict layered: each layer reads ONLY the previous layer,
                so the output gate cannot even see the raw inputs.
  skips=True    one-step skip: a gate reads the previous layer AND the raw
                inputs — but NOT any earlier hidden layer.
  skips='full'  full DAG: a gate reads every earlier layer AND the inputs.
                This is the standard circuit-complexity model.

`True` and `'full'` COINCIDE for depth <= 2 (with one hidden layer, the
previous layer is the only earlier layer), so they differ only on depth-3
and deeper architectures. MODEL CORRECTION (2026-07-26): `skips=True` was
described as "the standard model"; that is right at depth <= 2 and WRONG at
depth 3, where the standard model lets the output gate read the first
hidden layer too. Costs measured with skips=True on a depth>=3 architecture
are therefore UPPER bounds relative to the standard model, not minima in it.

CITATION CORRECTION (2026-07-22): this docstring used to call the wire
count "the Kane-Williams quantity". Kane-Williams (STOC 2016,
arXiv:1511.07860 §1) define LTF-of-LTF with an output gate that may take
input variables as well as previous LTF outputs — i.e. WITH skips. That is
a depth-2 statement, so skips=True and skips='full' both match it there.
The default layered path computes a strictly larger quantity.

Architecture is a list of hidden-layer sizes, of ANY length:
  []        -> depth-1 (single output gate reading inputs)
  [k]       -> depth-2 (k hidden gates -> output)
  [a, b]    -> depth-3 (a -> b -> output)
  [a, b, c] -> depth-4, and so on
Which architectures a certified sweep must actually search is a separate
question with its own answer; see `arch_family.live_archs`. Under
skips=False the output gate reads the last layer only (the inputs when
there is no hidden layer); under skips=True it also reads the inputs, and
under skips='full' every earlier layer as well. Activations of input
"gates" are the constant data bits.

Returns the minimum wire count and a circuit, or None if INFEASIBLE
for the given architecture and weight bound (an infeasibility is the
wire-cost lower-bound certificate, the analog of the LP refusal).
"""
import itertools
import os
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import csr_matrix
from ortools.sat.python import cp_model


# ---------------------------------------------------------------------------
# Truth-table helpers (assignment index a: bit i = (a >> i) & 1)
# ---------------------------------------------------------------------------

def truth_bits(T, n):
    """List of f(a) for a in 0..2^n-1 from packed integer T."""
    return [(T >> a) & 1 for a in range(1 << n)]


def assignments(n):
    """All input rows as tuples (x_0..x_{n-1})."""
    return [tuple((a >> i) & 1 for i in range(n)) for a in range(1 << n)]


# ---------------------------------------------------------------------------
# Variable pool for the MIP
# ---------------------------------------------------------------------------

class _Pool:
    def __init__(self):
        self.lb, self.ub, self.integ, self.cost = [], [], [], []
        self.n = 0

    def add(self, lb, ub, integer=True, cost=0.0):
        i = self.n
        self.lb.append(lb); self.ub.append(ub)
        self.integ.append(1 if integer else 0); self.cost.append(cost)
        self.n += 1
        return i


def _solve_arch(rows, labels, hidden, W, time_limit=None, mip_gap=0.0,
                wire_budget=None, feasibility_only=False):
    """Min-wire realization for one architecture. Returns dict with
    'wires', 'status', 'circuit' or None if proven infeasible.
    'status' in {'optimal','infeasible','timeout','error'}."""
    n = len(rows[0])
    P = len(rows)
    layers = list(hidden) + [1]          # hidden layers then the output gate
    prev_sizes = [n] + list(hidden)      # fan-in source size per layer
    pool = _Pool()

    # --- per-gate variables ---
    # gate[L][g] = dict(w=[idx...], b=idx, z=[idx...], y=[idx per point] or None)
    gate = []
    for L, size in enumerate(layers):
        psz = prev_sizes[L]
        Bmax = W * psz + 1
        glist = []
        is_output = (L == len(layers) - 1)
        for g in range(size):
            w = [pool.add(-W, W) for _ in range(psz)]
            z = [pool.add(0, 1, cost=1.0) for _ in range(psz)]   # objective
            b = pool.add(-Bmax, Bmax)
            y = None if is_output else [pool.add(0, 1) for _ in range(P)]
            glist.append(dict(w=w, b=b, z=z, y=y, psz=psz, Bmax=Bmax))
        gate.append(glist)

    rowsmat, lbs, ubs = [], [], []
    def add_row(coeffs, lb, ub):
        rowsmat.append(coeffs); lbs.append(lb); ubs.append(ub)

    INF = np.inf

    # --- indicator linking: |w| <= W*z  ->  w - W z <= 0 ; w + W z >= 0 ---
    for L in range(len(layers)):
        for gd in gate[L]:
            for wj, zj in zip(gd['w'], gd['z']):
                add_row({wj: 1.0, zj: -float(W)}, -INF, 0.0)
                add_row({wj: 1.0, zj:  float(W)}, 0.0, INF)

    # --- McCormick products u = w * a, where a is a previous-layer
    #     activation (binary var). Only needed when the source layer is
    #     hidden (L >= 1). For L == 0 the source is the constant inputs. ---
    def gate_pre(L, g, p):
        """Return (list of (coeff_or_varidx) terms) representing the linear
        pre-activation s = sum_j w_j * a_{j,p} + b for gate (L,g) at point p,
        as a dict {var: coeff} plus constant, adding McCormick vars/rows as
        needed."""
        gd = gate[L][g]
        terms = {gd['b']: 1.0}
        const = 0.0
        if L == 0:
            x = rows[p]
            for j, wj in enumerate(gd['w']):
                if x[j]:
                    terms[wj] = terms.get(wj, 0.0) + 1.0
        else:
            src = gate[L - 1]
            Wf = float(W)
            for j, wj in enumerate(gd['w']):
                a = src[j]['y'][p]                  # binary activation var
                u = pool.add(-W, W)
                # u <= W a ; u >= -W a ; u <= w + W(1-a) ; u >= w - W(1-a)
                add_row({u: 1.0, a: -Wf}, -INF, 0.0)
                add_row({u: 1.0, a:  Wf}, 0.0, INF)
                add_row({u: 1.0, wj: -1.0, a:  Wf}, -INF, Wf)
                add_row({u: 1.0, wj: -1.0, a: -Wf}, -Wf, INF)
                terms[u] = terms.get(u, 0.0) + 1.0
        return terms, const

    # --- activation constraints ---
    for L in range(len(layers)):
        is_output = (L == len(layers) - 1)
        for g in range(layers[L]):
            gd = gate[L][g]
            psz = gd['psz']; Bmax = gd['Bmax']
            M = W * psz + Bmax
            for p in range(P):
                terms, const = gate_pre(L, g, p)
                if is_output:
                    # s >= 0 if label 1 ; s <= -1 if label 0
                    if labels[p] == 1:
                        add_row(terms, -const, INF)        # s >= 0
                    else:
                        add_row(terms, -INF, -1.0 - const)  # s <= -1
                else:
                    # s - M*y in [-M, -1]
                    yv = gd['y'][p]
                    t = dict(terms); t[yv] = t.get(yv, 0.0) - float(M)
                    add_row(t, -float(M) - const, -1.0 - const)

    # --- symmetry breaking: within each hidden layer, wires
    #     non-increasing (sum z_g >= sum z_{g+1}) ---
    for L in range(len(hidden)):
        for g in range(layers[L] - 1):
            za = gate[L][g]['z']; zb = gate[L][g + 1]['z']
            add_row({**{i: 1.0 for i in za}, **{i: -1.0 for i in zb}}, 0.0, INF)

    # --- optional wire-budget constraint: total z <= budget ---
    all_z = [zi for L in range(len(layers)) for gd in gate[L] for zi in gd['z']]
    if wire_budget is not None:
        add_row({zi: 1.0 for zi in all_z}, 0.0, float(wire_budget))

    # --- assemble ---
    Nv = pool.n
    data, ri, ci = [], [], []
    for r, coeffs in enumerate(rowsmat):
        for v, c in coeffs.items():
            data.append(c); ri.append(r); ci.append(v)
    A = csr_matrix((data, (ri, ci)), shape=(len(rowsmat), Nv))
    cons = LinearConstraint(A, np.array(lbs), np.array(ubs))
    integrality = np.array(pool.integ)
    bounds = Bounds(np.array(pool.lb, float), np.array(pool.ub, float))
    options = {}
    if time_limit is not None:
        options['time_limit'] = time_limit
    if mip_gap is not None:
        options['mip_rel_gap'] = mip_gap

    obj = np.zeros(pool.n) if feasibility_only else np.array(pool.cost)
    res = milp(c=obj, constraints=cons, integrality=integrality,
               bounds=bounds, options=options)

    if res.status == 2:
        return dict(wires=None, status='infeasible', circuit=None)
    if res.status == 1 or res.x is None:
        # iteration/time limit without proven optimum
        return dict(wires=None, status='timeout', circuit=None)
    if res.status != 0:
        return dict(wires=None, status='error', circuit=None)

    x = res.x
    wires = int(round(res.fun))
    circuit = []
    for L in range(len(layers)):
        glist = []
        for gd in gate[L]:
            w = [int(round(x[i])) for i in gd['w']]
            b = int(round(x[gd['b']]))
            glist.append((w, b))
        circuit.append(glist)
    return dict(wires=wires, status='optimal', circuit=circuit)


def gate_sources(depth, L, skips=False):
    """Ordered source blocks read by a gate at layer L of a `depth`-layer
    circuit (depth = len(hidden) + 1). A block is an int (that layer's
    outputs) or the string 'x' (the raw inputs); a gate's weight vector is
    their concatenation IN THIS ORDER. Single source of truth shared by the
    solver, the evaluator and any external reader of a stored circuit.

      skips=False  -> [L-1]              strict layered
      skips=True   -> [L-1, 'x']         one-step skip (prev layer + inputs)
      skips='full' -> [0, 1, .., L-1, 'x']  full DAG (all earlier layers)

    The three coincide for L <= 1, so ONLY depth>=3 circuits distinguish
    True from 'full'."""
    if L == 0:
        return ['x']
    if skips == 'full':
        return list(range(L)) + ['x']
    return [L - 1] + (['x'] if skips else [])


def eval_circuit(circuit, x, skips=False):
    """Model-aware evaluation of a stored circuit on one input row."""
    n = len(x)
    acts = []
    for L, layer in enumerate(circuit):
        vals = []
        for (w, b) in layer:
            s, off = b, 0
            for blk in gate_sources(len(circuit), L, skips):
                src = list(x) if blk == 'x' else acts[blk]
                for j, v in enumerate(src):
                    if v:
                        s += w[off + j]
                off += len(src)
            vals.append(1 if s >= 0 else 0)
        acts.append(vals)
    return acts[-1][0]


def verify_circuit(circuit, T, n, skips=False):
    """Confirm a returned circuit computes truth table T exactly, IN THE
    STATED MODEL. Returns (ok, wires). Every oracle result should be
    verified, as the project verifies every IIS witness."""
    labels = truth_bits(T, n)
    rows = assignments(n)
    wires = sum(sum(1 for wj in w if wj != 0) for layer in circuit for (w, b) in layer)
    for p, x in enumerate(rows):
        if eval_circuit(circuit, x, skips=skips) != labels[p]:
            return False, wires
    return True, wires


def _solve_arch_cpsat(rows, labels, hidden, W, time_limit=None,
                      num_workers=None, wire_budget=None, feasibility_only=False,
                      cost_budget=None, l1_masks=None, require_active=False,
                      skips=False, support=None):
    """CP-SAT min-wire realization for one architecture (the primary
    backend; scipy/_solve_arch is kept only for cross-checking). Boolean
    activations make weight*activation products and reified thresholds
    native, no big-M. Returns dict(wires, status, circuit, wlb) where wlb is
    the integer wire LOWER bound (ceil of CP-SAT's best objective bound) when
    minimizing — so the caller can read the optimality gap (wires - wlb) on a
    timed-out incumbent, not just the upper bound. wlb is None when there is
    no wire objective (feasibility_only) or no model (infeasible).

    l1_masks (optional): list of allowed-input index sets, one per FIRST
    hidden-layer gate (len == hidden[0]); weights outside a gate's mask are
    fixed to 0. Used for fold-respecting (support-confined) synthesis. With
    masks, the wire-count symmetry-breaking order is applied only between
    adjacent SAME-mask gates (differently-masked gates are not
    interchangeable, so the global ordering would be unsound).

    support (optional): indices of the input variables the function actually
    depends on. Adds the IMPLIED constraint that each of them is read by at
    least one gate -- if no gate carries a nonzero weight on x_i then the
    output is independent of x_i, so any correct circuit must read it. The
    solver can only reach this through the correctness constraints, which is
    expensive; stating it directly is what makes the tight-budget
    infeasibility proofs tractable. At a wire budget of exactly s + G - 1
    (the floor: s input-sourced wires plus one out-wire per hidden gate) it
    forces the whole structure -- every support variable read exactly once,
    every hidden gate read exactly once, i.e. a read-once tree.
    Sound ONLY for the EXACT support. A superset is a real restriction and
    silently loses optima: it forces a wire to a variable the function does
    not depend on. (Caught by the self-check on 0x6666 at n=4, which is
    XOR(x0,x1) with support {0,1} -- declaring {1,2} raised its certified
    cost from 7 to 8.) Use arch_family.support_vars, never a guess.

    require_active (optional): force every gate to have >=1 in-wire and every
    hidden gate to be read by >=1 next-layer wire ('trimmed' circuits only).
    Sound ONLY when the caller enumerates all sub-architectures separately
    (the trimming lemma): any circuit trims to this
    form at <= cost by deleting unread gates and absorbing zero-wire
    (constant) gates into downstream biases."""
    n = len(rows[0]); P = len(rows)
    layers = list(hidden) + [1]
    prev_sizes = [n] + list(hidden)
    m = cp_model.CpModel()

    def src_blocks(L):
        """Ordered source blocks feeding any gate at layer L. A block is
        ('L', idx) for a hidden layer's outputs or ('x', None) for the raw
        inputs; the weight vector is their concatenation IN THIS ORDER.

        skips=False  : ['L', L-1]                  (strict layered)
        skips=True   : ['L', L-1] + ['x']          (one-step skip)
        skips='full' : ['L', 0..L-1] + ['x']       (full DAG)

        The three agree at L<=1, so depth<=2 results are model-independent
        once skips is on at all; only depth>=3 distinguishes True from
        'full'. See gate_sources()."""
        if L == 0:
            return [('x', None)]
        if skips == 'full':
            return [('L', l) for l in range(L)] + [('x', None)]
        return [('L', L - 1)] + ([('x', None)] if skips else [])

    def block_size(blk):
        return n if blk[0] == 'x' else layers[blk[1]]

    def src_pos(L_src, g, L_dst):
        """Index of layer-L_src gate g inside a layer-L_dst weight vector,
        or None when L_dst cannot read it in this model."""
        off = 0
        for blk in src_blocks(L_dst):
            if blk == ('L', L_src):
                return off + g
            off += block_size(blk)
        return None

    gate = []
    for L, size in enumerate(layers):
        # skip wires: layers past the first may also read the raw inputs
        # directly (the standard DAG model). Off by default so every
        # existing call keeps the strict layered semantics.
        psz = sum(block_size(b) for b in src_blocks(L)); Bmax = W * psz + 1
        is_out = (L == len(layers) - 1)
        glist = []
        for g in range(size):
            w = [m.NewIntVar(-W, W, f'w{L}_{g}_{j}') for j in range(psz)]
            z = [m.NewBoolVar(f'z{L}_{g}_{j}') for j in range(psz)]
            for wj, zj in zip(w, z):
                m.Add(wj == 0).OnlyEnforceIf(zj.Not())   # z=0 => w=0
            if L == 0 and l1_masks is not None:
                allowed = set(l1_masks[g])
                for j in range(psz):
                    if j not in allowed:
                        m.Add(w[j] == 0)
                        m.Add(z[j] == 0)
            b = m.NewIntVar(-Bmax, Bmax, f'b{L}_{g}')
            y = None if is_out else [m.NewBoolVar(f'y{L}_{g}_{p}')
                                     for p in range(P)]
            glist.append(dict(w=w, b=b, z=z, y=y))
        gate.append(glist)

    def preact_terms(L, g, p):
        gd = gate[L][g]
        terms = [gd['b']]
        x = rows[p]
        off = 0
        for blk in src_blocks(L):
            if blk[0] == 'x':                    # raw inputs: constant data
                for j in range(n):
                    if x[j]:
                        terms.append(gd['w'][off + j])
            else:
                src = gate[blk[1]]
                for j in range(len(src)):
                    a = src[j]['y'][p]
                    u = m.NewIntVar(-W, W, f'u{L}_{g}_{off + j}_{p}')
                    m.Add(u == gd['w'][off + j]).OnlyEnforceIf(a)
                    m.Add(u == 0).OnlyEnforceIf(a.Not())
                    terms.append(u)
            off += block_size(blk)
        return terms

    for L in range(len(layers)):
        is_out = (L == len(layers) - 1)
        for g in range(layers[L]):
            for p in range(P):
                s = sum(preact_terms(L, g, p))
                if is_out:
                    if labels[p] == 1:
                        m.Add(s >= 0)
                    else:
                        m.Add(s <= -1)
                else:
                    y = gate[L][g]['y'][p]
                    m.Add(s >= 0).OnlyEnforceIf(y)
                    m.Add(s <= -1).OnlyEnforceIf(y.Not())

    if support:
        # every relevant variable must be read SOMEWHERE (implied, not a
        # restriction -- see docstring). Also state the aggregate count,
        # which is what the LP relaxation can actually use against a budget.
        inslots = []
        for L in range(len(layers)):
            off = 0
            for blk in src_blocks(L):
                if blk[0] == 'x':
                    inslots.append((L, off))
                off += block_size(blk)
        for i in support:
            m.AddBoolOr([gate[L][g]['z'][off + i]
                         for (L, off) in inslots
                         for g in range(layers[L])])
        m.Add(sum(gate[L][g]['z'][off + i]
                  for (L, off) in inslots
                  for g in range(layers[L])
                  for i in support) >= len(support))

    if require_active:
        # trimmed-form constraints (see docstring for the soundness contract)
        for L in range(len(layers)):
            for g in range(layers[L]):
                m.AddBoolOr(gate[L][g]['z'])          # >=1 in-wire
        for L in range(len(layers) - 1):
            for g in range(layers[L]):
                # read SOMEWHERE downstream. Under skips='full' a hidden gate
                # may be read by any later layer, not only the next one, so
                # the next-layer-only version would be an unsound extra
                # constraint (it could cut off a legal trimmed circuit).
                lit = []
                for M in range(L + 1, len(layers)):
                    j = src_pos(L, g, M)
                    if j is not None:
                        lit += [gate[M][g2]['z'][j] for g2 in range(layers[M])]
                m.AddBoolOr(lit)

    # symmetry breaking: per hidden layer, gate wire-counts non-increasing
    # (with l1_masks, layer 0 orders only within same-mask runs — see docstring)
    for L in range(len(hidden)):
        for g in range(layers[L] - 1):
            if L == 0 and l1_masks is not None and \
                    set(l1_masks[g]) != set(l1_masks[g + 1]):
                continue
            m.Add(sum(gate[L][g]['z']) >= sum(gate[L][g + 1]['z']))

    allz = [zi for L in range(len(layers)) for gd in gate[L] for zi in gd['z']]
    if wire_budget is not None:
        m.Add(sum(allz) <= int(wire_budget))
    if cost_budget is not None:
        # cost = wires + active gates (a gate is active iff >=1 of its wires is
        # on) — the same quantity circuit_cost reports. Used as a feasibility
        # cap for arch pruning: "can this arch realize T at cost <= budget?"
        gacts = []
        for L in range(len(layers)):
            for gd in gate[L]:
                ga = m.NewBoolVar(f'ga{id(gd)}')
                m.AddMaxEquality(ga, gd['z'])        # ga = OR(wire indicators)
                gacts.append(ga)
        m.Add(sum(allz) + sum(gacts) <= int(cost_budget))
    if not feasibility_only:
        m.Minimize(sum(allz))

    solver = cp_model.CpSolver()
    if time_limit is not None:
        solver.parameters.max_time_in_seconds = float(time_limit)
    solver.parameters.num_search_workers = (
        num_workers or int(os.environ.get('MM_WORKERS', 0))
        or min(8, os.cpu_count() or 1))
    st = solver.Solve(m)

    def _extract():
        return [[([solver.Value(wj) for wj in gd['w']], solver.Value(gd['b']))
                 for gd in gate[L]] for L in range(len(layers))]
    def _wires(circ):
        return sum(1 for layer in circ for (w, b) in layer for wj in w if wj != 0)
    def _wlb():
        # integer wire lower bound from the objective (None if no objective)
        if feasibility_only:
            return None
        return int(-(-solver.BestObjectiveBound() // 1))      # ceil for a min

    if st == cp_model.OPTIMAL:
        circ = _extract()
        # with no objective (feasibility), OPTIMAL just means a model was found
        return dict(wires=_wires(circ),
                    status=('feasible' if feasibility_only else 'optimal'),
                    circuit=circ, wlb=None if feasibility_only else _wires(circ))
    if st == cp_model.INFEASIBLE:
        return dict(wires=None, status='infeasible', circuit=None, wlb=None)
    if st == cp_model.FEASIBLE:
        # incumbent found; for feasibility that IS success, else optimality unproven
        circ = _extract()
        return dict(wires=_wires(circ),
                    status=('feasible' if feasibility_only else 'timeout'),
                    circuit=circ, wlb=_wlb())
    return dict(wires=None, status='timeout', circuit=None, wlb=None)


def feasible_depth3_below(T, n, budget, W=2, archs=None, time_limit=None):
    """Win-check: is there a VERIFIED depth-3 circuit with <= budget wires?
    Sound for win-hunting — optimality is irrelevant, a verified circuit
    below W_2 is a win. Returns status in {'win','nowin','provisional'}:
      'win'         -> found+verified circuit (wires<=budget); fields set.
      'nowin'       -> every architecture PROVEN infeasible at the budget.
      'provisional' -> at least one architecture timed out (unknown)."""
    labels = truth_bits(T, n)
    rows = assignments(n)
    if archs is None:
        archs = [(2, 2), (2, 3), (3, 2), (3, 3), (2, 4), (4, 2)]
    any_timeout = False
    for (a, b) in archs:
        r = _solve_arch_cpsat(rows, labels, [a, b], W, time_limit=time_limit,
                              wire_budget=budget, feasibility_only=True)
        if r['status'] == 'feasible' and r['circuit'] is not None:
            ok, wires = verify_circuit(r['circuit'], T, n)
            if ok and wires <= budget:
                return dict(status='win', wires=wires, arch=(a, b),
                            circuit=r['circuit'])
        elif r['status'] == 'timeout':
            any_timeout = True
        # 'infeasible' -> this architecture cannot win; keep going
    return dict(status=('provisional' if any_timeout else 'nowin'),
                wires=None, arch=None, circuit=None)


def min_wires(T, n, hidden, W=3, time_limit=None, backend='cpsat',
              cost_budget=None, feasibility_only=False, l1_masks=None,
              skips=False, support=None):
    """Public: min wires for architecture `hidden` (list of hidden-layer
    sizes) realizing truth table T at weight bound W. cost_budget (wires +
    active gates) caps the realization for arch pruning — with
    feasibility_only it answers 'realizable at cost <= budget?' (status
    feasible / infeasible) without the optimize-to-proof cost. l1_masks:
    per-first-layer-gate allowed-input sets (fold-respecting synthesis);
    cpsat backend only."""
    labels = truth_bits(T, n)
    rows = assignments(n)
    if backend == 'cpsat':
        return _solve_arch_cpsat(rows, labels, hidden, W, time_limit=time_limit,
                                 cost_budget=cost_budget,
                                 feasibility_only=feasibility_only,
                                 l1_masks=l1_masks, skips=skips,
                                 support=support)
    assert support is None, "support requires the cpsat backend"
    assert l1_masks is None, "l1_masks requires the cpsat backend"
    assert not skips, "skips require the cpsat backend"
    return _solve_arch(rows, labels, hidden, W, time_limit=time_limit)


def best_shallow(T, n, W=3, kmax=None, time_limit=None):
    """C_2: min wires over depth-1 and depth-2 (k=1..kmax). Returns
    (wires, ('depth1'| ('depth2', k)), circuit). Stops when adding a gate
    no longer improves (objective is non-increasing in k)."""
    if kmax is None:
        kmax = n + 2
    best = None
    r1 = min_wires(T, n, [], W, time_limit)
    if r1['status'] == 'optimal':
        best = (r1['wires'], 'depth1', r1['circuit'])
    prev = None
    no_improve = 0
    for k in range(1, kmax + 1):
        r = min_wires(T, n, [k], W, time_limit)
        if r['status'] != 'optimal':
            continue
        cand = (r['wires'], ('depth2', k), r['circuit'])
        if best is None or cand[0] < best[0]:
            best = cand
        if prev is not None and r['wires'] >= prev:
            no_improve += 1
            if no_improve >= 1:        # one extra gate gave nothing more
                break
        else:
            no_improve = 0
        prev = r['wires']
    return best


def best_depth3(T, n, W=3, archs=None, cutoff=None, time_limit=None):
    """Min wires over a set of depth-3 (a,b) architectures. If cutoff is
    given, early-exit on the first architecture strictly below it (the
    'does a depth-3 win exist' question). Returns (wires, (a,b), circuit)
    or None if none feasible."""
    if archs is None:
        archs = [(2, 2), (2, 3), (3, 2), (3, 3), (2, 4), (4, 2)]
    best = None
    for (a, b) in archs:
        r = min_wires(T, n, [a, b], W, time_limit)
        if r['status'] != 'optimal':
            continue
        cand = (r['wires'], (a, b), r['circuit'])
        if best is None or cand[0] < best[0]:
            best = cand
        if cutoff is not None and cand[0] < cutoff:
            break
    return best


# ---------------------------------------------------------------------------
# NPN canonicalization
# ---------------------------------------------------------------------------

def npn_canonical(T, n, return_orbit_size=False):
    """Canonical NPN representative of truth table T: the minimum packed
    integer over all 2^n * n! * 2 transforms (input permutation, input
    negation, output negation). Sound as a cost invariant only because the
    gate set is general threshold functions, which IS closed under this
    group; see the module docstring."""
    N = 1 << n
    full = (1 << N) - 1
    orbit = set()
    best = None
    base_bits = [(T >> a) & 1 for a in range(N)]
    for perm in itertools.permutations(range(n)):
        for m in range(1 << n):                            # input-negation mask
            # build transformed (pre-output-negation) truth table
            val = 0
            for a in range(N):
                # old var j gets x^new_{perm[j]} XOR m_j
                b = 0
                for j in range(n):
                    bit = (a >> perm[j]) & 1
                    bit ^= (m >> j) & 1
                    b |= bit << j
                if base_bits[b]:
                    val |= (1 << a)
            for o in (0, 1):
                t = val ^ full if o else val
                orbit.add(t)
                if best is None or t < best:
                    best = t
    if return_orbit_size:
        return best, len(orbit)
    return best


def enumerate_npn_classes(n):
    """All NPN representatives at n inputs by bucketing every truth table.
    Feasible only for small n (n<=4: 2^16 tables)."""
    seen = {}
    N = 1 << n
    for T in range(1 << N):
        c = npn_canonical(T, n)
        seen.setdefault(c, 0)
        seen[c] += 1
    return seen


# ---------------------------------------------------------------------------
# XOR-fold structure of a win
# ---------------------------------------------------------------------------
# "Is a depth-3 win an XOR fold?" is a statement about the FUNCTION, not a
# gate: in a strict threshold circuit every gate is an LTF and every LTF is
# unate, so no gate is ever literally XOR. The decidable form is "is T the
# symmetric difference of two halfspaces, T = G1 XOR G2 for two LTFs". The
# fold may need halfspaces of higher weight than the per-gate bound Wmax, in
# which case the wire-optimal circuit synthesizes it from bounded-weight
# gates rather than showing an XOR at the top (e.g. 0x000707f8: fold weight
# 3, gate bound 2, NOR-topped circuit).

def _ltf_tables(n, W):
    """{truth_table_int: (w, b)} for every integer-weight LTF over n inputs
    with weights in [-W, W] (first representative per distinct table). An
    LTF outputs 1 at assignment a iff b + sum_j w_j * a_j >= 0."""
    N = 1 << n
    tables = {}
    Bmax = W * n
    for w in itertools.product(range(-W, W + 1), repeat=n):
        base = [sum(w[j] for j in range(n) if (a >> j) & 1) for a in range(N)]
        for b in range(-Bmax - 1, Bmax + 2):
            tt = 0
            for a in range(N):
                if base[a] + b >= 0:
                    tt |= (1 << a)
            tables.setdefault(tt, (w, b))
    return tables


def is_xor_of_two_ltfs(T, n, Wmax=4):
    """Decide whether T (n inputs, packed int) is the symmetric difference of
    two integer-weight threshold gates, T = G1 XOR G2, searching weight
    bounds W = 1..Wmax. Returns dict(found, weight, G1, G2):
      weight -> SMALLEST W at which both halfspaces are simultaneously
                realizable within [-W, W] (the fold's weight requirement);
      G1, G2 -> (w, b) representative witnesses at that W (not necessarily
                individually minimal-weight). found=False/weight=None if no
                fold exists up to Wmax.
    T = G1 XOR G2 because a XOR (a XOR T) = T; we look for a in the LTF set
    whose XOR-complement a^T is also an LTF."""
    for W in range(1, Wmax + 1):
        tables = _ltf_tables(n, W)
        keys = set(tables)
        for a in keys:
            if (a ^ T) in keys:
                return dict(found=True, weight=W,
                            G1=tables[a], G2=tables[a ^ T])
    return dict(found=False, weight=None, G1=None, G2=None)


_GATE2_NAMES = {
    (0, 0, 0, 0): 'FALSE', (1, 1, 1, 1): 'TRUE',
    (0, 0, 0, 1): 'AND',   (1, 1, 1, 0): 'NAND',
    (0, 1, 1, 1): 'OR',    (1, 0, 0, 0): 'NOR',
    (0, 1, 0, 1): 'x0',    (0, 0, 1, 1): 'x1',
    (1, 0, 1, 0): 'NOT x0', (1, 1, 0, 0): 'NOT x1',
    (0, 1, 0, 0): 'x0 AND NOT x1', (0, 0, 1, 0): 'NOT x0 AND x1',
    (1, 1, 0, 1): 'x0 OR NOT x1',  (1, 0, 1, 1): 'NOT x0 OR x1',
    (0, 1, 1, 0): 'XOR',   (1, 0, 0, 1): 'XNOR',  # not LTF-realizable
}


def decode_output_gate(circuit):
    """Classify the Boolean function the output (last-layer) gate computes
    over its k layer-inputs. Returns dict(k, weights, bias, truth, name,
    unate). `truth` is the 2^k-tuple of gate outputs in index order (input j
    = bit j of the index); `name` is a recognized gate name for k<=2, else
    'threshold-k'. unate is always True (XOR/XNOR are not threshold
    functions); a non-trivial XOR character lives in the realized FUNCTION,
    never in a single gate."""
    w, b = circuit[-1][0]
    k = len(w)
    truth = tuple(1 if b + sum(w[j] for j in range(k) if (m >> j) & 1) >= 0
                  else 0 for m in range(1 << k))
    if k == 1:
        name = {(0, 1): 'x0', (1, 0): 'NOT x0',
                (0, 0): 'FALSE', (1, 1): 'TRUE'}.get(truth, '?')
    elif k == 2:
        name = _GATE2_NAMES.get(truth, 'threshold-2')
    else:
        name = f'threshold-{k}'
    return dict(k=k, weights=list(w), bias=b, truth=truth, name=name,
                unate=True)


# ---------------------------------------------------------------------------
# Evaluation-cost objective: alpha*wires + beta*gates
# ---------------------------------------------------------------------------
# Cost model = inference cost of a Boolean threshold net. Each wire is one
# weighted input (a conditional add on binary activations); each gate is one
# accumulate+compare. So cost = alpha*(#wires) + beta*(#gates); (alpha,beta)
# is a deployment dial ((1,1) = binary net, (2,1) = real-multiply MAC). Unused
# zero-fanin gates do not count. NOTE: weights are NOT penalized here (a
# multiply is a multiply regardless of |w|), so this is run at FREE weights
# (a generous bound W); the min-wire weight cap Wmax=2 is a different,
# weight-starved object -- under a cap a function can be forced onto a wider
# or deeper architecture purely for want of magnitude, so capped and free
# costs are not comparable point by point.

def circuit_cost(circuit, alpha=1, beta=1):
    """Evaluation-cost score alpha*wires + beta*gates of a threshold circuit.
    wires = total nonzero input weights; gates = gates with >=1 wire (a
    zero-fanin gate does no work and is not counted). Returns
    dict(cost, wires, gates, maxw)."""
    wires = gates = maxw = 0
    for layer in circuit:
        for (w, b) in layer:
            nz = sum(1 for wj in w if wj != 0)
            if nz:
                wires += nz
                gates += 1
                maxw = max(maxw, max(abs(wj) for wj in w))
    return dict(cost=alpha * wires + beta * gates, wires=wires,
                gates=gates, maxw=maxw)


def best_cost(T, n, alpha=1, beta=1, W=4, archs=None, time_limit=None,
              skips=False):
    """Minimum alpha*wires + beta*gates realization of T over a set of
    architectures (each a list of hidden-layer sizes; [] = depth-1). Run at
    free-ish weights via bound W (default 4; record maxw and bump W if it is
    tight). Returns the best dict(cost, arch, wires, gates, maxw, circuit,
    per) or None. Soundness of min-over-archs: min_wires(arch) returns a
    realization with wires<=any same-arch optimum, and every smaller used
    architecture is itself enumerated, so the architecture-wise minimum is
    the true global min cost (the argument is the one just given)."""
    if archs is None:
        archs = ([[]] + [[k] for k in range(1, n + 3)]
                 + [[a, b] for a in range(2, n + 1) for b in (2, 3)])
    best = None
    per = {}
    for h in archs:
        r = min_wires(T, n, list(h), W, time_limit, skips=skips)
        if r['status'] != 'optimal' or r['circuit'] is None:
            continue
        c = circuit_cost(r['circuit'], alpha, beta)
        per[tuple(h)] = c['cost']
        if best is None or c['cost'] < best['cost']:
            best = dict(cost=c['cost'], arch=tuple(h), wires=c['wires'],
                        gates=c['gates'], maxw=c['maxw'], circuit=r['circuit'])
    if best is not None:
        best['per'] = per
    return best


if __name__ == '__main__':
    import sys, time
    # quick self-test
    print("NPN class counts (known: n2=4, n3=14, n4=222):")
    for n in (2, 3):
        print(f"  n={n}: {len(enumerate_npn_classes(n))}")

    # XOR-fold structure of the depth-3 wins
    print("\nXOR-fold test (is_xor_of_two_ltfs); expected fold weights "
          "constructed=1, 0x1eee0f=2, 0x707f8=3 (only 0x707f8 exceeds the "
          "Wmax=2 gate bound):")
    for name, T, exp in [('constructed 0xfffe0001', 0xfffe0001, 1),
                         ('sampled     0x1eee0f', 0x001eee0f, 2),
                         ('sampled     0x707f8', 0x000707f8, 3)]:
        r = is_xor_of_two_ltfs(T, 5, Wmax=4)
        ok = r['found'] and r['weight'] == exp
        print(f"  {name}: found={r['found']} weight={r['weight']} "
              f"{'OK' if ok else 'MISMATCH (exp %d)' % exp}")
        assert ok, (name, r)

    # decode_output_gate on the verified 0x707f8 win (NOR top) -- no solve
    win707f8 = [[([-1, -1, -2, 0, 0], 1), ([0, 0, 0, 1, 0], -1),
                 ([0, 0, 0, 0, -1], 0)],
                [([-2, 1, -2], 1), ([1, -2, 2], -3)],
                [([-1, -2], 0)]]
    vok, wires = verify_circuit(win707f8, 0x000707f8, 5)
    assert vok and wires == 13, (vok, wires)
    og = decode_output_gate(win707f8)
    print(f"  0x707f8 win verified (wires={wires}); output gate = "
          f"{og['name']} truth={og['truth']} (expected NOR (1, 0, 0, 0))")
    assert og['name'] == 'NOR' and og['truth'] == (1, 0, 0, 0), og

    # cost objective: the verified 0x707f8 win circuit has 13 wires, 6 gates
    cc = circuit_cost(win707f8, alpha=1, beta=1)
    print(f"  circuit_cost(0x707f8 win, (1,1)): {cc}")
    assert cc == dict(cost=19, wires=13, gates=6, maxw=2), cc
    # best_cost of AND3 (0x80, n=3) = one depth-1 gate: 3 wires, 1 gate
    bc = best_cost(0x80, 3, alpha=1, beta=1, W=3)
    print(f"  best_cost(AND3, (1,1)): cost={bc['cost']} arch={bc['arch']} "
          f"wires={bc['wires']} gates={bc['gates']}")
    assert bc['cost'] == 4 and bc['arch'] == () and bc['gates'] == 1, bc

    # ---- three circuit models (2026-07-26) ----
    # (a) the source layout is what the solver and the evaluator agree on
    assert gate_sources(3, 2, False) == [1]
    assert gate_sources(3, 2, True) == [1, 'x']
    assert gate_sources(3, 2, 'full') == [0, 1, 'x']
    assert gate_sources(2, 1, True) == gate_sources(2, 1, 'full') == [0, 'x']
    # (b) skips=True and skips='full' must AGREE at depth<=2, and each must
    #     reproduce its own known value. XOR(x0,x1) at n=4: 9 layered, 7 skip.
    for sk, exp in ((False, 9), (True, 7), ('full', 7)):
        best = None
        for h in ([], [1], [2], [3]):
            r = min_wires(0x6666, 4, h, W=7, time_limit=60, skips=sk)
            if r['circuit'] is None or r['status'] != 'optimal':
                continue
            ok, _ = verify_circuit(r['circuit'], 0x6666, 4, skips=sk)
            assert ok, ('model eval disagrees with solver', sk, h)
            c = circuit_cost(r['circuit'])['cost']
            best = c if best is None else min(best, c)
        print(f"  XOR2@n4 depth<=2, skips={sk!r}: cost {best} (expected {exp})")
        assert best == exp, (sk, best, exp)
    # (c) at depth 3 they need NOT agree, but 'full' can never be worse:
    #     it strictly extends the source set of the output gate.
    r1 = min_wires(0x6666, 4, [1, 1], W=7, time_limit=120, skips=True)
    r2 = min_wires(0x6666, 4, [1, 1], W=7, time_limit=120, skips='full')
    c1 = circuit_cost(r1['circuit'])['cost'] if r1['circuit'] else None
    c2 = circuit_cost(r2['circuit'])['cost'] if r2['circuit'] else None
    ok2 = r2['circuit'] is None or verify_circuit(r2['circuit'], 0x6666, 4,
                                                  skips='full')[0]
    print(f"  XOR2@n4 arch [1,1]: one-step {c1} ({r1['status']}) vs "
          f"full-DAG {c2} ({r2['status']})")
    assert ok2, 'full-DAG circuit failed its own verifier'
    assert not (c1 is not None and c2 is not None and r2['status'] == 'optimal'
                and c2 > c1), (c1, c2)

    # ---- support= is an IMPLIED constraint FOR THE TRUE SUPPORT: it may
    # change the TIME, never the answer. Checked on certified optima.
    def _support(T, n):
        return [i for i in range(n)
                if any(((T >> a) & 1) != ((T >> (a ^ (1 << i))) & 1)
                       for a in range(1 << n))]
    assert _support(0x6666, 4) == [0, 1], _support(0x6666, 4)
    for T, n, h, exp in ((0x6666, 4, [1], 7),                   # XOR2@n4
                         (0x6666, 4, [2], 7),
                         (0x80, 3, [], 4)):                     # AND3
        sup = _support(T, n)
        a = min_wires(T, n, h, W=7, time_limit=120, skips='full')
        b = min_wires(T, n, h, W=7, time_limit=120, skips='full', support=sup)
        ca = circuit_cost(a['circuit'])['cost'] if a['circuit'] else None
        cb = circuit_cost(b['circuit'])['cost'] if b['circuit'] else None
        print(f"  support-invariance 0x{T:x} n={n} {h}: {ca} vs {cb} "
              f"(expected {exp})")
        assert ca == cb == exp, (T, h, ca, cb, exp)
        assert a['status'] == b['status'], (T, h, a['status'], b['status'])
    print("self-check OK")
