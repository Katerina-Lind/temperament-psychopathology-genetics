# GenomicSEM covariance validation

## Source

Seven-trait developmental joint LDSC analysis.

Traits:
EMO, ACT, SHY, CADHD, AINT, ASUD, ARISK.

## Matrix structure

- S: 7 x 7 genetic covariance matrix.
- V: 28 x 28 sampling covariance matrix.
- S_Stand: standardised genetic covariance matrix.
- V_Stand: standardised sampling covariance matrix.

The analysis uses column-major lower-triangular
covariance parameter ordering, including the diagonal.

## Validation results

The standardised matrix S_Stand equals cov2cor(S).

Selected standard errors from V_Stand:

| Pair | SE |
|---|---:|
| ACT–SHY | 0.075513 |
| ACT–CADHD | 0.088274 |
| AINT–ASUD | 0.030963 |
| SHY–ARISK | 0.039170 |

These reproduce the standard errors reported by
the joint GenomicSEM LDSC analysis.

## Statistical inference

The manuscript contrasts were calculated using
the unstandardised covariance matrix S and its
sampling covariance matrix V.

Genetic correlations were obtained by standardising S.

Contrast standard errors were estimated using
the multivariate delta method, propagating uncertainty
in genetic covariances and heritabilities.

These standard errors can differ from those obtained
directly from V_Stand.

## Remaining limitations

- Validation is based on the same LDSC estimates.
- The full contrast procedure has not been independently
  reproduced using a second software implementation.
- Multiple-testing corrections are exploratory.
