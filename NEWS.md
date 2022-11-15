# jackstraw 1.3.2 (2020-12-23)

* Fixed lots of minor things that caused `check` errors
  * Vignette: minor edits to fix recent build errors
  * Moved external package imports from R code to `DESCRIPTION`, adding missing packages, and fixed namespaces.
  * Corrected function arguments and usage discrepancies.
* Added this `NEWS.md` file to track changes to the package.

# jackstraw 1.3.3 (2021-03-09)

Overview:

- Added unit tests to almost all functions
  - Basic uses are covered (default optional parameters), and some but not all advanced uses
- Minor user-facing changes: 
  - Moved mandatory arguments to the front (no backward-incompatible problems since these arguments had to be called by name unless the earlier optional arguments were also present and were passed unnamed (rare))
  - Added more informative error messages for missing mandatory arguments, objects with wrong classes, and covariate dimension disagreements
  - Minor bug fixes (some `s = 1` edge cases)
  - Eliminated all non-error/warning verbosity when `verbose = FALSE`.

Additional details:

- Function `jackstraw_cluster`
  - Added option `pool` (default `TRUE`) to calculate p-values by pooling null statistics (to match `jackstraw_kmeans` option and default).  Previously the hardcoded behavior matched `pool = FALSE`.
- Functions `jackstraw_lfa`, `jackstraw_alstructure`, `jackstraw_pca`, `jackstraw_rpca`, `jackstraw_irlba`, `jackstraw_subspace`:
  - Argument `r` is now second argument.
    - `lfa`, `alstructure`, `subspace` versions: `r` does not have a default value and it is mandatory.
- Functions `jackstraw_kmeans`, `jackstraw_kmeanspp`, `jackstraw_MiniBatchKmeans`, `jackstraw_pam`, `jackstraw_cluster`:
  - Debugged `s = 1` edge case: null data used to be centered incorrectly.
    Bug only occurred in combination with `center = TRUE` (default).
- Internal (unexported) function `RSS` now returns actual residual sum of squares.  This change does not affect any exported functions that use it.  Previously `RSS` calculated a normalized version (equal to `1 - R^2`), but this normalization canceled out in `FSTAT` (its only downstream use), so the normalization had no user-facing effect.

Exclusive list of functions without unit tests (all are redundant with other packages, so they are candidates for removal in the near future):

- Exported: 
  - `lfa.corpcor` (redundant with `lfa::lfa`)
  - `pi0est_bootstrap` (redundant with `qvalue::pi0est`)
- Internal, no exported backward dependencies:
  - `getp` (redundant with `qvalue::empPvals`)
  - `devdiff_parallel` (redundant with `gcatest::gcat.stat`)

# jackstraw 1.3.4.9000 (2021-04-14)

- Functions `jackstraw_pca`, `jackstraw_rpca`, `jackstraw_irlba`:  Corrected documentation (parameter `r1` was incorrectly described as `PC` in parts of the documentation.  Thanks to Djordje Bajić (GitHub username `djbajic`) for reporting this error!
- Removed option `seed` from all functions that had it.  For the same behavior, call `set.seed(seed)` before calling the function.
- Functions `jackstraw_lfa` and `jackstraw_alstructure`: removed `devR` option.
- Removed redundant functions
  - `lfa.corpcor`: same as `lfa::lfa` with option `override = TRUE`
  - `pi0est_bootstrap`: redundant with `qvalue::pi0est` with option `pi0.method = 'bootstrap'`
  - `dev.R` (internal; functionality implemented in package `gcatest`)
  - `devdiff_parallel` (internal; redundant with `gcatest::gcat.stat`)
  - `getp` (internal; redundant with `qvalue::empPvals`)

# jackstraw 1.3.5.9000 (2021-05-14)

- Function `jackstraw_lfa` now accepts genotypes input as `BEDMatrix` objects.
  In this case, the function operates on a low-memory mode, keeping data on disk rather than memory as much as possible, and writes permuted data into temporary files as well.
  To enable this mode, the `BEDMatrix` and `genio` packages are now dependencies.
  Note only `jackstraw_lfa` supports `BEDMatrix` because `lfa` supports it too (most recent fork; see below).
- Removed function `devdiff`, which is redundant (and replaced internally) with `gcatest::delta_deviance_lf`, a function that supports more special cases, including genotypes accessed through a `BEDMatrix` object.
  The only internal dependencies were `jackstraw_lfa` and `jackstraw_alstructure`.
- Updated `README.md` to instruct users to install the most updated forks of `lfa` and `gcatest` on GitHub (under username `alexviiia`), rather than the Bioconductor versions that are lacking critical updates.

# jackstraw 1.3.6.9000 (2022-02-08)

- All `jackstraw_*` functions now return `NA` p-values for `NA` statistics.
  - Before `NA` statistics resulted in p-values of 1 instead, which is what `qvalue::empPvals` returns.  Now an internal wrapper function ensures the desired behavior.
- Removed two `jackstraw_pam` toy example unit tests that failed often due to colinearity.
- Reformatted this `NEWS.md` slightly to improve its automatic parsing.

# jackstraw 1.3.7 (2022-11-10)

- Temporary changes for CRAN resubmission
  - Removed `jackstraw_alstructure` because dependency `alstructure` is not on CRAN or Bioconductor.
  - Removed BEDMatrix functionality for function `jackstraw_lfa` because the latest versions of the dependencies `lfa` and `gcatest` on Bioconductor do not support BEDMatrix (the development versions that do support BEDMatrix are on GitHub only).
  - Removed `lfa` and `gcatest` devel version requirements.
  - Added internal functions `delta_deviance_lf`, `delta_deviance_snp`, `delta_deviance_snp_lf`, which are copies of the `gcatest` functions of the same name (available on development version only, hence this copying)
- Made testing a bit more lenient towards some NA cases
- Minor non-code edits
  - Spell checked documentation
  - Changed to single maintainer (Neo)
  - Removed `VignetteBuilder: knitr` since there's no vignette anymore
  - Fixed broken or outdated URLs

# jackstraw 1.3.8 (2022-11-15)

- Temporarily removed exported function `jackstraw_lfa` and package dependencies `lfa` and `gcatest`, since `lfa` is having build errors on bioc-devel.
  - Also removed internal functions `pseudo_Rsq`, `mcfadden_Rsq_snp`, `efron_Rsq`, `efron_Rsq_snp`, `delta_deviance_{lf,snp,snp_lf}.R` (all depended on `lfa`)
  - Also removed package dependencies `parallel` and `genio` which are temporarily not being used.
- Removed many package suggested dependencies (`knitr`, `rmarkdown`, `ggplot2`, `mutoss`, `Matrix`, `gridExtra`, `cowplot`, `scales`, `formatR`) that were only used in a vignette currently not being built.
- Reduced example dimensions for `jackstraw_irlba` and `jackstraw_rpca` by 5 to keep their runtime low.
