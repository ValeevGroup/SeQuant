# Provenance tool for tests/unit/ga_reference.hpp: regenerating the reference
# numbers requires the GA_scratch Python prototype chain (proto_csv_*) this
# script imports; it is checked in as the authoritative record of where those
# numbers came from, not as a runnable part of the build.
# -*- coding: utf-8 -*-
"""Dump exact reference numbers from the prototype chain for the C++ port."""
import json, time
from proto_csv_core import Diagram, decode_tree, encode_tree, tree_ranges
from proto_csv_sched import Problem
import proto_csv_flat as FL
import proto_csv_run as R
import proto_csv_run2 as R2
from proto_csv_gonly import build_g_scope, per_term_opt_dp

out = {}


def measure(name, scope_fn):
    t0 = time.time()
    terms, scope = scope_fn()
    D = [Diagram(t, i) for i, t in enumerate(scope)]
    P = Problem(D)
    print('%s: problem built in %.0fs' % (name, time.time() - t0), flush=True)
    F = FL.Flat(P, method='greedy')
    per_term, g_opt = [], []
    for d in range(P.m):
        c, code = per_term_opt_dp(P, d)
        per_term.append(c)
        g_opt += code
    l2_flat = sum((len(blks) - 1) * P.tface[t] for t, blks in P.targets)
    h0 = []
    for t, blks in P.targets:
        h0 += [0] * len(tree_ranges(len(blks)))
    rec = {
        'terms': [d.name for d in D],
        'n_nodes': sum(d.n for d in D),
        'per_term_opt': per_term,
        'l2_flat': l2_flat,
        'base_a': sum(per_term) + l2_flat,
        'g_opt': g_opt, 'h0': h0,
        'seed_greedy': F(g_opt, h0),
        'n_keys': len(set(P.U.key.values())),
        'key_pairs': len(P.U.key),
    }
    out[name] = rec
    print('%s: base_a=%d seed_greedy=%d keys=%d/%d' % (
        name, rec['base_a'], rec['seed_greedy'], rec['n_keys'],
        rec['key_pairs']), flush=True)
    return P, D, F


# ---- f-only ----------------------------------------------------------------
P, D, Fg = measure('fonly', R.build_scope)
Fx = FL.Flat(P, method='exact')
out['fonly']['seed_exact'] = Fx(out['fonly']['g_opt'], out['fonly']['h0'])
g_ctrl = R2.control_genome(P, D)
l2fam = frozenset([frozenset([i]) for i in range(6)] +
                  [frozenset([0, 1]), frozenset([4, 5]),
                   frozenset([0, 1, 4, 5]), frozenset([2, 3]),
                   frozenset(range(6))])
h_ctrl = list(encode_tree(l2fam, list(range(6))))
out['fonly']['g_ctrl'] = list(g_ctrl)
out['fonly']['h_ctrl'] = h_ctrl
out['fonly']['control_exact'] = Fx(g_ctrl, h_ctrl)
out['fonly']['control_greedy'] = Fg(g_ctrl, h_ctrl)
print('fonly: seed_exact=%d control_exact=%d control_greedy=%d' % (
    out['fonly']['seed_exact'], out['fonly']['control_exact'],
    out['fonly']['control_greedy']), flush=True)

# ---- g-only ----------------------------------------------------------------
P, D, Fg = measure('gonly', build_g_scope)
W = json.load(open('gonly_winner.json'))
tot, l1, l2, vals, pick, dem, st = Fg(W['g'], W['h'], detail=True)
out['gonly']['winner_g'] = W['g']
out['gonly']['winner_h'] = W['h']
out['gonly']['winner_total'] = tot
out['gonly']['winner_l1'] = l1
out['gonly']['winner_l2'] = l2
out['gonly']['winner_n_arrays'] = len(set(pick.values()))
print('gonly: winner=%d l1=%d l2=%d arrays=%d' % (tot, l1, l2,
      out['gonly']['winner_n_arrays']), flush=True)

json.dump(out, open('ga_reference.json', 'w'), indent=1)
print('written ga_reference.json', flush=True)
