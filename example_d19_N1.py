from sage.all import *
from JL_algorithm import (
    maximal_order_congruences,
    hecke_eigenvalues_by_prime,
)


def run_example():
    a = -1
    b = -19
    N = 1
    prime_bound = 9

    maximal_prop = maximal_order_congruences(a, b)

    nontrivial_prop = lambda A, B, C, D: (
        (A + C) % 4 == 0 and
        (B - D) % 4 == 0
    )

    data = hecke_eigenvalues_by_prime(
        a, b, maximal_prop, nontrivial_prop, N, prime_bound
    )

    S = CuspForms(-b, 2)
    f = S.newforms()[0]
    qf = f.q_expansion(prime_bound + 1)

    rows = []
    for p in sorted(data):
        rows.append((p, data[p], qf[p]))

    return rows


def print_table(rows):
    print("{:<6} {:<26} {:<10}".format("p", "Calculated eigenvalues", "Real a_p"))
    print("-" * 42)
    for p, eigs, ap in rows:
        print("{:<6} {:<26} {:<10}".format(str(p), str(eigs[1]), str(ap)))


if __name__ == "__main__":
    rows = run_example()
    print_table(rows)