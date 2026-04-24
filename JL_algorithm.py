from sage.all import *
import itertools
import numpy as np


def algebra_info(a, b):
    A = QuaternionAlgebra(QQ, a, b)
    i, j, k = A.gens()
    return A, i, j, k


def maximal_info(a, b):
    A, i, j, k = algebra_info(a, b)
    O = A.maximal_order()
    return O


def get_units(a, b):
    units = []
    A, i, j, k = algebra_info(a, b)
    for x, y, z, w in itertools.product([-1, 0, 1], repeat=4):
        q = A(x) + A(y)*i + A(z)*j + A(w)*k
        if q.reduced_norm() == 1:
            units.append(q)
    return units


def number_of_units(a, b):
    return len(get_units(a, b))


def sols_of_norm(a, b, m):
    sols = []
    A, i, j, k = algebra_info(a, b)
    bound = ceil(sqrt(m)) + 1

    for x in range(-bound, bound):
        for y in range(-bound, bound):
            for z in range(-bound, bound):
                for w in range(-bound, bound):
                    q = x + y*i + z*j + w*k
                    if q.reduced_norm() == m:
                        sols.append((x, y, z, w, q))
    return sols


def max_order_solutions_of_norm(a, b, m, maximal_prop):
    all_sols = sols_of_norm(a, b, m)
    return [
        (x, y, z, w, q) for (x, y, z, w, q) in all_sols
        if maximal_prop(x, y, z, w)
    ]


def nontrivial_filter_norm(a, b, m, maximal_prop, nontrivial_prop):
    return [
        (x, y, z, w, q)
        for (x, y, z, w, q)
        in max_order_solutions_of_norm(a, b, m, maximal_prop)
        if nontrivial_prop(x, y, z, w)
    ]


def orbits_mod_maximal_order(a, b, m, maximal_prop):
    sols = max_order_solutions_of_norm(a, b, m, maximal_prop)
    num_units = number_of_units(a, b)
    return len(sols) / num_units


def nontrivial_orbits_mod_maximal_order(a, b, m, maximal_prop, nontrivial_prop):
    sols = nontrivial_filter_norm(a, b, m, maximal_prop, nontrivial_prop)
    num_units = number_of_units(a, b)
    return len(sols) / num_units


def hecke_operator_matrix(a, b, maximal_prop, nontrivial_prop, p):
    N = 4 * p
    Tpf11 = orbits_mod_maximal_order(a, b, N, maximal_prop)
    Tpfc1 = p + 1 - Tpf11
    Tpf1c = nontrivial_orbits_mod_maximal_order(a, b, 2 * N, maximal_prop,
                                                nontrivial_prop)
    Tpfcc = p + 1 - Tpf1c
    return [[Tpf11, Tpfc1], [Tpf1c, Tpfcc]]


def maximal_order_congruences(a, b):
    O = maximal_info(a, b)
    M = Matrix(QQ, [[c for c in e] for e in O.basis()])
    L = lcm([M[i, j].denominator() for i in range(4) for j in range(4)])
    LM = (L * M).change_ring(ZZ)
    coeffs = [[LM[i, j] for i in range(4)] for j in range(4)]
    residues = [[int(coeffs[j][i]) % L for i in range(4)] for j in range(4)]

    checks = []
    for pvec in itertools.product(range(L), repeat=4):
        if all(p == 0 for p in pvec):
            continue
        signed = [p - L if p % L > L // 2 else p % L for p in pvec]
        if next((p for p in signed if p != 0), None) != 1:
            continue
        if all(sum(pvec[j] * residues[j][i] for j in range(4)) % L == 0 for i in range(4)):
            checks.append((pvec, L))

    lams = [
        lambda x, p=pvec, m=mod: sum(p[j] * x[j] for j in range(4)) % m == 0
        for pvec, mod in checks
    ]
    return lambda A, B, C, D, lams=lams: all(lam((A, B, C, D)) for lam in lams)


def bad_primes(a, b, N):
    A, i, j, k = algebra_info(a, b)
    d = A.discriminant()
    finite_primes = []
    for p in Primes():
        if p > abs(d) and p > abs(N):
            break
        finite_primes.append(p)

    bad = []
    for p in finite_primes:
        if d % p == 0 or N % p == 0:
            bad.append(p)

    return bad


def hecke_eigenvalues_by_prime(a, b, maximal_prop, nontrivial_prop, N, prime_bound):
    bad = set(bad_primes(a, b, N))
    data = {}

    for p in Primes():
        if p > prime_bound:
            break
        if p in bad:
            continue

        mat = hecke_operator_matrix(a, b, maximal_prop, nontrivial_prop, p)
        M = Matrix(QQ, mat)
        eigs = M.eigenvalues()
        data[p] = eigs

    return data