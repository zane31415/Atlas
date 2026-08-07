"""ARCHITECTURE FAMILY — which shapes a certified sweep must actually search.

Written 2026-07-26 after the full-DAG correction found that BOTH of the
things we were treating as fixed were model-dependent:

  1. the wiring model (fixed: mm_oracle now has skips='full'), and
  2. the ARCHITECTURE FAMILY. `d4_metric`'s 32 shapes stop at depth 3, and
     `d5_metric`'s stop at a curated depth-3 set. Both cutoffs were argued
     in the LAYERED model, where a one-gate layer DEGENERATES (a gate that
     reads only one bit computes that bit or its negation). With skips it
     does not degenerate, and the measurements followed immediately:
     parity-4 is 17 over depth<=3 but **16 on [1,1,1]** (verified, and 15
     is proven infeasible there). A hand-picked family cannot be called
     exhaustive once its justification is gone.

So the family is not a list any more. It is DERIVED: every architecture
whose lower bound leaves it able to beat the incumbent, of any depth and
any width. That set is finite, small, and provably complete.

THE BOUNDS (all from the trimming lemma, FOLDPRICE_log 2026-07-11: any
circuit trims at no greater cost to one where every gate has >= 1 in-wire
and every hidden gate is read downstream; the trimmed circuit sits on a
narrower architecture, and every narrower architecture is enumerated
separately, so pruning on these is sound).

  (a) SUPPORT bound, new and the one that closes the depth question.
      Let s = |support(f)| and G = sum(h) + 1 = the gate count of a trimmed
      circuit on h. Every one of the s relevant inputs must be read by at
      least one gate, or the circuit ignores it: that is >= s wires whose
      SOURCE is an input. Every one of the G-1 hidden gates must be read
      downstream: that is >= G-1 wires whose SOURCE is a gate. Those two
      wire sets are disjoint, so
          wires >= s + G - 1     and     cost >= s + 2G - 1.
      This is what makes the family finite: G <= (incumbent + 1 - s) / 2,
      i.e. a bound on the TOTAL number of gates, hence on both depth and
      width at once. It is tight on depth-1 LTFs (s inputs, 1 gate, cost
      s+1).
  (b) DEPTH-2, h = [k]:   cost >= 3k + 1.
      k hidden gates need >= 1 in-wire each, and each must be read by the
      output, with k+1 gates in total.
  (c) DEPTH-3, h = [a, b]: cost >= 2a + 2b + max(a, b) + 1.
      Re-derived for the full DAG model (SKIP_log 2026-07-26) and unchanged
      from the one-step derivation: with W2, W3 the in-wire counts of layer
      2 and the output, W2 >= b, W3 >= b, and every layer-1 gate is read by
      layer 2 or by the output so W2 + W3 - b >= a; minimising W2 + W3
      gives b + max(a, b).

cost_lb takes the MAX of whichever apply. (a) and (b) cross at k = s - 4;
(a) and (c) cross at max(a, b) = s.

NOT used: the STRONG_L1 layer-1 bound, walked back 2026-07-24 as proven
only at free weights (SKIP_log). Everything here is weight-independent.
"""

MAX_SANITY_GATES = 12          # guard: refuse to enumerate absurd families


def support_vars(T, n):
    """Indices of the variables f actually depends on. Feed this to
    mm_oracle's `support=` so the solver is told outright that a correct
    circuit must read each of them -- it is the difference between a proof
    and a timeout on tight-budget architectures."""
    out = []
    for i in range(n):
        bit = 1 << i
        if any(((T >> a) & 1) != ((T >> (a ^ bit)) & 1) for a in range(1 << n)):
            out.append(i)
    return out


def support_size(T, n):
    """Number of variables f actually depends on."""
    return len(support_vars(T, n))


def cost_lb(h, s):
    """Sound lower bound on the (1,1) cost of any TRIMMED circuit on
    architecture h computing a function of support size s. See module
    docstring for the three arguments; valid in the layered, one-step and
    full-DAG models alike (none of them appeals to who may read whom)."""
    h = list(h)
    G = sum(h) + 1
    lb = s + 2 * G - 1
    if not h:
        # depth-1: one gate reading s inputs, cost s+1 -- and 0 for the
        # constant function, which trims to no gates at all. Return 0 so the
        # bound holds without a special case; [] is visited first anyway.
        return 0 if s == 0 else lb
    if len(h) == 1:
        lb = max(lb, 3 * h[0] + 1)
    elif len(h) == 2:
        a, b = h
        lb = max(lb, 2 * a + 2 * b + max(a, b) + 1)
    return lb


# Max |w| needed to realize ANY threshold function of f variables.
# f <= 6 are EXACT, computed by enumerating sorted-nonneg weight vectors
# (which cover every NPN class, since permuting/negating inputs is an NPN
# operation and preserves max|w|) and finding where the realizable set stops
# growing: 1, 1, 2, 3, 5, 9 for f = 1..6.
# f >= 7 fall back to Muroga-Toda-Takasu, (f+1)^((f+1)/2) / 2^f, which is a
# proven UPPER bound and therefore safe -- just not tight. Computing f=7
# exactly would need ~10^6 weight vectors x 10^2 thresholds; the literature
# value is smaller than 32, so this is where the remaining slack lives.
EXACT_MAXW = {0: 0, 1: 1, 2: 1, 3: 2, 4: 3, 5: 5, 6: 9}


def muroga_maxw(f):
    """Muroga-Toda-Takasu upper bound, ceil((f+1)^((f+1)/2) / 2^f). Computed,
    not tabulated -- a hand-written table drifted from the formula at f=9
    (197 vs 196) and the self-check caught it."""
    import math
    return math.ceil((f + 1) ** ((f + 1) / 2) / 2 ** f)


def max_fanin(cost, s, n, is_ltf=True):
    """Largest fan-in any gate can have in a TRIMMED circuit of cost <= cost
    on any architecture live for this class.

    Two limits, whichever is smaller: STRUCTURAL (a gate at layer L reads the
    earlier layers plus the n inputs) and BUDGETARY (all G gates are wired, so
    one gate can hold at most cost - G - (G-1) = cost - 2G + 1 of the wires).

    Not circular when `cost` is a certified value: a CHEAPER circuit has a
    tighter budget and hence a smaller fan-in, so this covers every circuit
    in the search space, not just the optimum.

    is_ltf=False sharpens the budget limit by one, and it is the difference
    between needing Muroga's f=8 bound (77) and its f=7 bound (32). TIGHTNESS
    ARGUMENT: if some gate's fan-in EQUALS cost - 2G + 1 then every other gate
    has exactly one in-wire, so each computes a constant or a literal of its
    single source — and by induction down the DAG, a literal or constant of a
    single RAW INPUT. Then either the wide gate is the output, in which case
    the output is an LTF of literals of inputs, i.e. an LTF of the inputs; or
    it is hidden, in which case the OUTPUT has one in-wire and the circuit is
    a literal of an LTF of literals, again an LTF. So a function that is
    proven NOT an LTF cannot attain the bound, and its fan-in is at most
    cost - 2G. (Same device that closed parity-4's [5] and [3,3] in the
    STRONG_L1 walkback, SKIP_log 2026-07-24.)"""
    best = 0
    for h in live_archs(cost + 1, s):
        G = sum(h) + 1
        layers = list(h) + [1]
        budget = cost - 2 * G + 1 - (0 if is_ltf else 1)
        for L in range(len(layers)):
            src = (sum(layers[:L]) if L else 0) + n
            best = max(best, min(src, budget))
    return best


def required_W(cost, s, n, is_ltf=True):
    """Weight bound at which the search PROVABLY loses nothing, so its answer
    is the true free-weight minimum rather than a W-relative one.

    Every gate is an LTF on at most max_fanin() of its sources. Replace each
    gate's weights by a minimal-max-weight realization of that gate's own
    function on its own support (zero elsewhere): the computed function is
    unchanged and wires never increase, since a weight can only be driven to
    zero, which lands on a narrower architecture that is enumerated
    separately. So a bound of maxw(max_fanin) excludes no circuit of cost
    <= cost.

    The point of the exact table: at n=5 the sweep already ran at W=9, which
    IS the exact requirement at fan-in 6 -- so every class whose gates cannot
    reach fan-in 7 is already free-weight certified, with no further compute."""
    f = max_fanin(cost, s, n, is_ltf=is_ltf)
    return EXACT_MAXW[f] if f in EXACT_MAXW else muroga_maxw(f)


def live_archs(incumbent, s, include_depth1=True):
    """EVERY architecture that could realise cost < incumbent, at any depth
    and any width -- the complete search family for a certified sweep.

    Finiteness: cost_lb >= s + 2G - 1, so G <= (incumbent - s) / 2 rounded,
    bounding the total gate count; the architectures are then the
    compositions of 1..Gmax-1 hidden gates into layers.

    Returned cheapest-bound-first, which is the order that lands a good
    incumbent early and lets the caller retire the rest without solving."""
    gmax = (incumbent + 1 - s) // 2          # cost_lb < incumbent => G <= gmax
    if gmax > MAX_SANITY_GATES:
        raise ValueError(f'family too large: incumbent={incumbent} s={s} '
                         f'implies up to {gmax} gates')
    out = []
    if include_depth1 and cost_lb([], s) < incumbent:
        out.append([])
    # compositions of t hidden gates into ordered layers, t = 1 .. gmax-1
    for t in range(1, max(gmax, 1)):
        stack = [[]]
        while stack:
            h = stack.pop()
            rem = t - sum(h)
            if rem == 0:
                if cost_lb(h, s) < incumbent:
                    out.append(list(h))
                continue
            for k in range(1, rem + 1):
                stack.append(h + [k])
    out.sort(key=lambda h: (cost_lb(h, s), len(h), tuple(h)))
    return out


if __name__ == '__main__':
    # self-check: the bounds must not contradict any CERTIFIED value we hold
    checks = [
        # (arch, support, known cost, source)
        ([], 4, 5, 'depth-1 LTF on 4 inputs: 4 wires + 1 gate -- bound tight'),
        ([1], 2, 7, 'XOR2@n4 skip, N4_SKIP_FULL'),
        ([1, 1], 4, 13, '0x17ac full-DAG, 2026-07-26'),
        ([1, 1], 4, 17, 'parity-4 one-step, N4_SKIP_FULL'),
        ([1, 1, 1], 4, 16, 'parity-4 full-DAG, 2026-07-26 (15 infeasible)'),
        ([2], 4, 9, 'N4_SKIP_FULL band'),
    ]
    for h, s, c, why in checks:
        lb = cost_lb(h, s)
        print(f'  cost_lb({h}, s={s}) = {lb:2d} <= {c}  [{why}]')
        assert lb <= c, (h, s, lb, c, why)

    print('\nfamily sizes (the point: a certified sweep is finite)')
    for inc, s, lbl in ((17, 4, 'n=4 parity-4 incumbent'),
                        (13, 4, 'n=4 typical'),
                        (16, 5, 'n=5 worst known'),
                        (12, 5, 'n=5 modal'),
                        (8, 5, 'n=5 cheapest')):
        fam = live_archs(inc, s)
        deep = sum(1 for h in fam if len(h) >= 3)
        print(f'  incumbent {inc:2d}, s={s}: {len(fam):3d} archs '
              f'(max depth {max(len(h) for h in fam) + 1}, {deep} of depth>=4) '
              f'[{lbl}]')
    # the family must contain the shapes that actually won
    assert [1, 1, 1] in live_archs(17, 4), 'family misses the parity-4 winner'
    assert [1, 1] in live_archs(14, 4), 'family misses the 0x17ac winner'
    print('\nself-check OK')
