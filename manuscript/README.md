# Prime-shift extensions of classical ménage identities

The note proves Mathar's six-step recurrence in OEIS A000159 and van Hoeij's
prime-shift congruence. Its main theorem covers every hit degree `1 <= h < p`
for every prime `p` and every index `n >= 0`; zero hits have range `n >= 1`.
The coverage table gives the exact OEIS contribution locators and the neighboring
statements that follow. The classical one-hit and two/three-hit identities are
credited to Kaplansky and Riordan, *The problème des ménages* (1946), p. 122,
equation (33), rather than presented as new.
The bivariate factorial-series generating function was recorded by Vladeta
Jovović in A058087, first FORMULA, 14 December 2009, in descending row order;
the note gives the change of variables to hit-degree order explicitly.

## Files and rebuilding

- `note.tex`: one main LaTeX source with its complete bibliography.
- `note.pdf`: the mathematical note.
- `metadata.json`: title, creators, affiliations, license and claim identity.
- `models.py`: the rational-affine factorial-series construction and complete
  source recurrence coefficient arrays.
- `verify_core.py`: the three-hit differential identities and finite checks.
- `verify.py`: sibling recurrence certificates, congruence comparisons and
  characteristic-two/three rational identities.
- `verification.txt`: expected concise checker output.

Build with Tectonic 0.17.0:

```sh
tectonic note.tex
```

Alternatively, run `pdflatex note.tex` twice. The bibliography is embedded,
so BibTeX is unnecessary. Standard AMS, geometry, booktabs, longtable,
lmodern and hyperref packages suffice.

Run the supporting checks with Python 3 and SymPy 1.14.0:

```sh
python3 verify.py
```

The expected last line is `ALL_SIBLING_CHECKS_OK`.

## Exact identities and finite ranges

The rational calculations establish the displayed Euler identities,
factorial-series reduction, differential identities and recurrence-operator
identities. `models.py` constructs the coefficients from the hit degree,
not by fitting sequence terms. For the six sibling recurrences it obtains
zero as the rational remainder, except for A000185's six-step recurrence,
whose remainder is `-3202597896*x**5`. Coefficient extraction gives all-index
recurrences on the domains in the note.

The finite comparisons independently reconstruct values from the defining
rook sum. Their ranges are:

- Defining polynomial rows 0 through 5; all eight hit degrees 0 through 7
  through index 300.
- Mathar's six-step recurrence at indices 6 through 300, including the
  zero-extension indices 6, 7 and 8.
- A000425 at 4 through 300; A000033 at 5 through 300; the third-order
  three-hit recurrence at 6 through 300; A000181 at 8 through 300; A000185's
  six-step and four-step recurrences at 11 and 9 through 300, respectively.
- The corrected one-hit identity at 1 through 300, and the two/three-hit
  relation at 3 through 299.
- 18,660 uniform-congruence comparisons: the 25 primes at most 97, hit
  degrees below `min(p,8)`, and `0 <= n <= 100`, except that zero hits
  begin at `n = 1`.
- 2,500 quotient-congruence comparisons: those 25 primes and `1 <= n <= 100`.
- Odd-hit parity through index 300; the two-hit characteristic-two
  congruence at 0 through 100; prime-square antiperiods whenever all
  required indices are at most 300; and the three-hit nine-shift identity
  at 0 through 291.
- The zero-hit four-prime-shift relation for every tested prime and every
  positive index with all terms at most 300.

The exact characteristic-two and characteristic-three rational identities
are also checked after clearing denominators in the appropriate finite field.
Finite comparisons are distinct from the written all-index proofs.

## Formal coverage

The pinned declaration
[`MatharA000159Recurrence.result`](https://github.com/the-omega-institute/trureturing/blob/8acba483bc80695fef00c6b8c491e2ae16d5145f/D5/S1/Recurrence/Invariants/MatharA000159Recurrence.lean#L656)
states the recurrence for the rational coefficient definition of `a(n)`
for every `n >= 9`. The same file defines `U`, `a` and `claim`.
The zero-extension cases, uniform prime shifts, antiperiods, quotient
congruence, corrected indexing and other recurrences have the written
proofs or exact rational certificates identified in the note.
The Bessel rounding bounds are outside its conclusions.

## Creators and license

Haobo Ma — ChronoAI Pte Ltd, Singapore.

Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore.

The note and original supporting material are licensed under Creative Commons
Attribution 4.0 International, CC BY 4.0 (`cc-by-4.0`).
