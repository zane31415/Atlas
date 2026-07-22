"""M-M depth-3 parameter-cost oracle.

Minimum-WIRE (nonzero-weight) realization of a Boolean function as a
feedforward integer-weight THRESHOLD circuit of a given architecture,
via scipy.optimize.milp (HiGHS). General threshold gates at every layer
including the output (required for NPN soundness, MM_log DP-1). The
threshold/bias is free; cost = total nonzero input weights = total
wires.

MODEL (load-bearing): STRICT LAYERED. Each layer reads only the previous
layer, so the output gate cannot see the raw inputs — no skip
connections. The standard circuit-complexity model allows them and the
restriction is not free; see the model note in README.md. Earlier
revisions of this docstring called the wire count "the Kane-Williams
quantity"; that was wrong, because Kane-Williams (STOC 2016,
arXiv:1511.07860) let the output gate read input variables as well as
previous gate outputs. What this oracle computes is an upper bound on
their quantity, in the layered model.

Architecture is a list of hidden-layer sizes:
  []     -> depth-1 (single output gate reading inputs)
  [k]    -> depth-2 (k hidden gates -> output)
  [a, b] -> depth-3 (a -> b -> output)
The output gate reads the last layer (the inputs when there is no
hidden layer). Activations of input "gates" are the constant data bits.

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


def verify_circuit(circuit, T, n):
    """Confirm a returned circuit computes truth table T exactly. Returns
    (ok, wires). Every oracle result should be verified, as the project
    verifies every IIS witness."""
    labels = truth_bits(T, n)
    rows = assignments(n)
    wires = sum(sum(1 for wj in w if wj != 0) for layer in circuit for (w, b) in layer)
    for p, x in enumerate(rows):
        act = list(x)
        for layer in circuit:
            nxt = []
            for (w, b) in layer:
                s = b + sum(w[j] * act[j] for j in range(len(w)))
                nxt.append(1 if s >= 0 else 0)
            act = nxt
        if act[0] != labels[p]:
            return False, wires
    return True, wires


def _solve_arch_cpsat(rows, labels, hidden, W, time_limit=None,
                      num_workers=None, wire_budget=None, feasibility_only=False,
                      cost_budget=None, l1_masks=None):
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
    interchangeable, so the global ordering would be unsound)."""
    n = len(rows[0]); P = len(rows)
    layers = list(hidden) + [1]
    prev_sizes = [n] + list(hidden)
    m = cp_model.CpModel()

    gate = []
    for L, size in enumerate(layers):
        psz = prev_sizes[L]; Bmax = W * psz + 1
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
        if L == 0:
            x = rows[p]
            for j, wj in enumerate(gd['w']):
                if x[j]:
                    terms.append(wj)
        else:
            src = gate[L - 1]
            for j, wj in enumerate(gd['w']):
                a = src[j]['y'][p]
                u = m.NewIntVar(-W, W, f'u{L}_{g}_{j}_{p}')
                m.Add(u == wj).OnlyEnforceIf(a)
                m.Add(u == 0).OnlyEnforceIf(a.Not())
                terms.append(u)
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
              cost_budget=None, feasibility_only=False, l1_masks=None):
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
                                 l1_masks=l1_masks)
    assert l1_masks is None, "l1_masks requires the cpsat backend"
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
# NPN canonicalization (MM_log section 2)
# ---------------------------------------------------------------------------

def npn_canonical(T, n, return_orbit_size=False):
    """Canonical NPN representative of truth table T: the minimum packed
    integer over all 2^n * n! * 2 transforms (input permutation, input
    negation, output negation). General-threshold soundness in MM_log."""
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
# XOR-fold structure of a win (D5_log 2026-06-15)
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
# Evaluation-cost objective: alpha*wires + beta*gates (D4METRIC_log 2026-06-15)
# ---------------------------------------------------------------------------
# Cost model = inference cost of a Boolean threshold net. Each wire is one
# weighted input (a conditional add on binary activations); each gate is one
# accumulate+compare. So cost = alpha*(#wires) + beta*(#gates); (alpha,beta)
# is a deployment dial ((1,1) = binary net, (2,1) = real-multiply MAC). Unused
# zero-fanin gates do not count. NOTE: weights are NOT penalized here (a
# multiply is a multiply regardless of |w|), so this is run at FREE weights
# (a generous bound W); the min-wire weight cap Wmax=2 is a different,
# weight-starved object (see the cap-stability audit in D5_log).

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


def best_cost(T, n, alpha=1, beta=1, W=4, archs=None, time_limit=None):
    """Minimum alpha*wires + beta*gates realization of T over a set of
    architectures (each a list of hidden-layer sizes; [] = depth-1). Run at
    free-ish weights via bound W (default 4; record maxw and bump W if it is
    tight). Returns the best dict(cost, arch, wires, gates, maxw, circuit,
    per) or None. Soundness of min-over-archs: min_wires(arch) returns a
    realization with wires<=any same-arch optimum, and every smaller used
    architecture is itself enumerated, so the architecture-wise minimum is
    the true global min cost (proof in D4METRIC_log)."""
    if archs is None:
        archs = ([[]] + [[k] for k in range(1, n + 3)]
                 + [[a, b] for a in range(2, n + 1) for b in (2, 3)])
    best = None
    per = {}
    for h in archs:
        r = min_wires(T, n, list(h), W, time_limit)
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

    # XOR-fold structure of the D5 wins (re-derives D5_log 2026-06-15 table)
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
    print("self-check OK")
