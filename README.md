# goldbach-cantarini

Code, outputs and paper for **"A missing term in an explicit formula for Goldbach averages"**
(Giovane de Morais, Instituto Tecnológico de Aeronáutica, ITA).

## Result

Let Λ be the von Mangoldt function, r_G(n) = Σ_{m1+m2=n} Λ(m1)Λ(m2), and let ρ run over the
nontrivial zeros of ζ. M. Cantarini, *Explicit formula for the average of Goldbach numbers*
(Indian J. Math. 61 (2019), 253-279; preprint [arXiv:1801.08475](https://arxiv.org/abs/1801.08475)),
proves in Theorem 1 an explicit formula for the average Σ'_{n≤2N} r_G(n) with a remainder
F(N) = O(N).

The paper shows that the bookkeeping of the proof drops the term **+P(N)²**, where
P(N) = Σ_ρ N^ρ/ρ. With this term restored the formula holds with remainder O(N). Without it
the remainder satisfies F(N) = P(N)² + O(N), so **F(N) ≠ O(N)** (Littlewood's oscillation
theorem gives limsup F(N)/N = +∞). Consequently the estimate in Cantarini's Conjecture 2 fails.
The correction relies on Lemmas 5 and 6 of Cantarini's paper, which are not re-proved.

Numerically (two independent implementations, 1969 values of N up to 6·10⁶, no truncation of
the sums over zeros): F(N)/N = c + β·P(N)²/N with β = 0.999933 (implementation A) and
0.999875 (implementation B), and c = −7.3515, against −4 log(2π) = −7.351508.

**Note.** The analysis uses the **arXiv version** (v1) of Cantarini's paper. The published
version in Indian J. Math. 61 was not available to the author, so numbering and statements
may differ there.

The paper is a **preprint, to be submitted**. Status and DOI will be added here.

## Contents

| Path | What it is |
|---|---|
| `paper/` | `main.tex`, `main.pdf`, `refs.bib`, `main.bbl`, `figs/` (PDF figures and `make_figs.py`) |
| `certificate/data/zeros1` | First 10⁵ zeros of ζ (ordinates γ), from Odlyzko's tables |
| `certificate/implementation_A/` | Implementation A: `exact_F.py` (exact F(N), no truncation), `fit_F.py` (fit of F/N against P²/N), `beta_identity.py` (Beta-function identity, mpmath), `outer_sum.py` (limit of the outer sum over zeros), with outputs |
| `certificate/implementation_B/` | Second, independent implementation (shares no code with A): `core.py`, `exact_F.py`, `fit_F.py`, `beta_identity.py`, `zeros_check.py` (T(N) against the double series truncated at K zeros), with outputs |
| `certificate/SHA256SUMS` | SHA-256 of every file in `certificate/` |

Output files are named after the script that wrote them plus the parameters
(for example `exact_F_5e6.out`).

## Requirements

- Python ≥ 3.12 (developed with 3.14.4)
- `numpy`, `scipy` (`implementation_A/outer_sum.py`), `mpmath` (identity checks, `core.py`), `matplotlib` (figures only)
- A LaTeX installation with `pdflatex` and `bibtex` for the paper

```bash
pip install numpy scipy mpmath matplotlib
```

No GPU needed. Memory stays below 4 GB (`ulimit -v 4000000`) and each run uses at most four
threads. Times below are on a few CPU cores of the author's workstation.

## Reproducing

Run from `certificate/`. Redirect the output to a scratch directory to compare with the
shipped files without overwriting them.

```bash
cd certificate
# Beta-function identity, mpmath at 30 digits (seconds)
python3 -I implementation_A/beta_identity.py       > beta_A.out
python3 -B implementation_B/beta_identity.py       > beta_B.out
# Limit of the outer sum over zeros, 1e5 zeros (seconds)
python3 -I implementation_A/outer_sum.py data/zeros1 > outer_sum.out
# Implementation A: exact F(N) at 1560 values of N, then the fit (about 2 min)
python3 -I implementation_A/exact_F.py 5000000 1500 > exact_F_5e6.out
python3 -I implementation_A/fit_F.py exact_F_5e6.out
# Implementation B: exact F(N) at 409 values of N, then the fit (about 4 min)
python3 -B implementation_B/exact_F.py 2e5 6e6 280 20261017 > exact_F_2e5_6e6.tsv
python3 -B implementation_B/fit_F.py exact_F_2e5_6e6.tsv
# T(N) against the double series truncated at K zeros (slow: the first takes well over 30 min, the second about 13 min of CPU)
python3 -B implementation_B/zeros_check.py data/zeros1 25,58,74 1000,4000,16000,40000 > zeros_check_N25_58_74.out
python3 -B implementation_B/zeros_check.py data/zeros1 25,27,49,26,50 1000,4000,16000 > zeros_check_prime_powers.out
```

The last printed digit of a few rows of `exact_F_5e6.out` can differ between runs (summation order depends on the thread count); the fits are identical.

Expected key numbers: Beta identity relative error ≤ 4.07e−29 (implementation B) and
≤ 1.9e−29 (A); fit `F/N = -7.351498 + 0.999933*P^2/N` (A) and
`-7.351500 + 0.999875*P^2/N` (B).

Figures: `cd paper/figs && python3 -I make_figs.py` (reads the shipped outputs).
Paper: `cd paper && pdflatex main && bibtex main && pdflatex main && pdflatex main`.

## Verifying the checksums

```bash
cd certificate && sha256sum -c SHA256SUMS
```

## License

Code and text: MIT (see `LICENSE`). The zeros in `certificate/data/zeros1` are Odlyzko's
table, see <https://www-users.cse.umn.edu/~odlyzko/zeta_tables/>.
