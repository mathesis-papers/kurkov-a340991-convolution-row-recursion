# Kurkov's convolution-row conjectures over arbitrary semirings

Haobo Ma (ChronoAI Pte Ltd, Singapore) and Wenlin Zhang (National University of Singapore; ChronoAI Pte Ltd, Singapore).

## Result

For an arbitrary sequence in a unital semiring, set `A(z)=sum_{j>=1} a_j z^j` and `P_n(x)=sum_{r=0}^n [z^n]A(z)^r x^r`. Kurkov's array `R(1,k)=a_k`, `R(n,k)=R(n-1,k+1)+x R(n-1,1) a_k` satisfies `P_n(x)=x R(n,1)` for every positive `n`. The coefficient formula for every positive column is proved in `note.tex`.

The conjecture is the FORMULA line signed by Mikhail Kurkov, 5 April 2025, in OEIS A340991 revision 24. The written argument includes noncommutative coefficients with the indicated product order. The declaration `result` in the pinned Lean file proves the commutative-semiring case; `kurkov_coeff` is the corresponding coefficient invariant:

https://github.com/the-omega-institute/trureturing/blob/03150b25ffad5ced1bebdca2efc8c05b484c1eca/D5/S1/Recurrence/Algebraic/KurkovConvolutionRowRecursion.lean

The same signed conjecture in OEIS A357583 revision 10 concerns Bell numbers. Corollary 3 proves that specialization and recovers the established Bell row sums A129247 and central terms A357584. Corollary 4 deduces Bala's factorial generating function and continued fraction in A090238 from the classical INVERT identity and Euler's factorial fraction, with Bernstein--Sloane and Flajolet cited for those ingredients. Corollary 5 proves the binomial relation between shifted and unshifted triangles and recovers the factorial and Catalan coefficient recurrences. The note also gives both left and right finite-support row recurrences. Its coverage table identifies the statements proved, classical formulas recovered and topics outside the shifted-array identity.

The coefficient definition `RecursiveBoundaryRowSeries.tail` in the following pinned file is exactly the shifted tail `B_k` used in the proof:

https://github.com/the-omega-institute/trureturing/blob/03150b25ffad5ced1bebdca2efc8c05b484c1eca/D5/S1/Recurrence/Algebraic/RecursiveBoundaryRowSeries.lean

The noncommutative extension, Bell and factorial deductions, finite-support recurrences and nonzero-constant-term obstruction have written proofs in the note. The two Lean files supply the stated commutative theorem, its coefficient invariant and its tail definition; they are not a formalization of all the corollaries.

## Build the PDF

Use Tectonic 0.17.0:

```sh
tectonic --untrusted note.tex
```

Alternatively, a standard LaTeX installation with `amsmath`, `amssymb`, `amsthm`, `geometry`, `lmodern`, `booktabs`, `longtable` and `hyperref` suffices:

```sh
pdflatex -halt-on-error note.tex
pdflatex -halt-on-error note.tex
```

The bibliography is embedded in `note.tex`; no external bibliography tool or data download is needed.

## Reproduce the finite identities

Python 3.10 or later and its standard library suffice:

```sh
python3 verify.py --output finite-checks.recomputed.json
cmp finite-checks.json finite-checks.recomputed.json
```

The program compares independent convolution powers with the shifted array, tail coefficient formula and both ordered row recurrences. It covers the prime sequence, every length-six sequence over `{-1,0,1}`, Boolean and characteristic-two examples, noncommuting matrices, and free noncommuting symbols. It also compares six Bell/factorial/Catalan triangle variants to OEIS data, checks Bell row sums and central terms, and computes both factorial continued fractions and Bala's rational series with exact bivariate coefficients through degree 12. `finite-checks.json` gives exact ranges. These are finite checks; the all-index statements have the written proofs and exact classical references in the note. `fixtures.json` contains numerical data from the cited OEIS entries, with source revisions and offsets. The shifted Catalan triangle A059365 has a different row-zero convention, recorded explicitly.

## Files and license

`note.tex` contains the theorem, complete proof, consequences, coverage table and references. `metadata.json` gives the scientific metadata. `references.json` identifies the primary reference versions. `verify.py`, `fixtures.json` and `finite-checks.json` provide exact finite computations and their numerical source data.

The note and original accompanying code are licensed under Creative Commons Attribution 4.0 International (`CC BY 4.0`). Numerical OEIS fixtures are attributed to their respective entries; the OEIS source content is governed by the OEIS End-User License Agreement. See `LICENSE.txt`.
