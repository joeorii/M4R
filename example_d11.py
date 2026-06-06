"""
This file tests the functions implemented in file JL_algorithm.py for the
quaternion algebra of discriminant 11 and level 1. It computes a table of
quaternionic Hecke eigenvalues for good primes and compares them with the 
corresponding classical modular-form coefficients.
"""
from sage.all import *
import pandas as pd

from JL_algorithm import (
    maximal_order_congruences,
    hecke_eigenvalues_by_prime,
)


def compute_dataframe():
    """
    This function returns a table containing the good primes, the quaternionic
    Hecke eigenvalues, and the corresponding classical Fourier coefficients.
    """
    a = -1
    b = -11
    N = 1
    prime_bound = 20
    maximal_prop = maximal_order_congruences(a, b)
    nontrivial_prop = lambda A, B, C, D: (
        (A + C + 2 * D) % 4 == 0 and
        (B + 2 * C - D) % 4 == 0
    )
    data = hecke_eigenvalues_by_prime(
        a, b, maximal_prop, nontrivial_prop, N, prime_bound
    )
    S = CuspForms(-b, 2)
    f = S.newforms()[0]
    qf = f.q_expansion(prime_bound + 1)

    rows = []
    for p in sorted(data):
        eigs = data[p]

        eig2 = eigs[1] if len(eigs) > 1 else ""

        rows.append({
            "Prime": p,
            "EV": eig2,
            "a_p": qf[p],
        })
    return pd.DataFrame(rows)


def print_table(df):
    """
    df is the table produced by the function compute_dataframe. 
    This function prints this table.
    """
    print("\nComputed data:\n")
    print(df.to_string(index=False))

def main():
    """This function constructs and prints the table."""
    df = compute_dataframe()
    print_table(df)


if __name__ == "__main__":
    main()