from sage.all import *
import pandas as pd

from JL_algorithm import (
    maximal_order_congruences,
    hecke_eigenvalues_by_prime,
)


def compute_dataframe():
    # Quaternion algebra parameters
    a = -1
    b = -19

    # Level parameter for bad primes
    N = 1

    # Compute up to this prime
    prime_bound = 200

    # Maximal order congruence conditions
    maximal_prop = maximal_order_congruences(a, b)

    # Extra congruence condition defining the nontrivial class
    nontrivial_prop = lambda A, B, C, D: (
        (A + C + 2 * D) % 4 == 0 and
        (B + 2 * C - D) % 4 == 0
    )

    # Compute Hecke eigenvalues from your quaternionic setup
    data = hecke_eigenvalues_by_prime(
        a, b, maximal_prop, nontrivial_prop, N, prime_bound
    )

    # Compute the classical modular form coefficients for comparison
    S = CuspForms(-b, 2)
    f = S.newforms()[0]
    qf = f.q_expansion(prime_bound + 1)

    rows = []
    for p in sorted(data):
        eigs = data[p]

        eig1 = eigs[0] if len(eigs) > 0 else ""
        eig2 = eigs[1] if len(eigs) > 1 else ""

        rows.append({
            "Prime": p,
            "Trivial EV": eig1,
            "Nontrivial EV": eig2,
            "a_p": qf[p],
        })

    return pd.DataFrame(rows)


def print_table(df):
    print("\nComputed data:\n")
    print(df.to_string(index=False))


def save_latex_table(
    df,
    filename="JLdisc11.tex",
    caption="Hecke eigenvalue comparison",
    label="tab:hecke"
):
    df_latex = df.copy()

    # Put entries in math mode for nicer Overleaf output
    for col in df_latex.columns:
        df_latex[col] = df_latex[col].apply(lambda x: f"${x}$")

    tabular = df_latex.to_latex(index=False, escape=False)

    full_table = (
        "\\begin{table}[ht]\n"
        "\\centering\n"
        + tabular +
        f"\\caption{{{caption}}}\n"
        f"\\label{{{label}}}\n"
        "\\end{table}\n"
    )

    with open(filename, "w") as f:
        f.write(full_table)

    print(f"\nLaTeX table written to {filename}")


def main():
    df = compute_dataframe()
    print_table(df)
    save_latex_table(df)


if __name__ == "__main__":
    main()