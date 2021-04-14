# 2020-12-23 - jackstraw 1.3.2

* Fixed lots of minor things that caused `check` errors
  * Vignette: minor edits to fix recent build errors
  * Moved external package imports from R code to `DESCRIPTION`, adding missing packages, and fixed namespaces.
  * Corrected function arguments and usage discrepancies.
* Added this `NEWS.md` file to track changes to the package.

# 2021-03-09 - jackstraw 1.3.3

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
    Bug only occured in combination with `center = TRUE` (default).
- Internal (unexported) function `RSS` now returns actual residual sum of squares.  This change does not affect any exported functions that use it.  Previously `RSS` calculated a normalized version (equal to `1 - R^2`), but this normalization canceled out in `FSTAT` (its only downstream use), so the normalization had no user-facing effect.

Exclusive list of functions without unit tests (all are redundant with other packages, so they are candidates for removal in the near future):

- Exported: 
  - `lfa.corpcor` (redundant with `lfa::lfa`)
  - `pi0est_bootstrap` (redundant with `qvalue::pi0est`)
- Internal, no exported backward dependencies:
  - `getp` (redundant with `qvalue::empPvals`)
  - `devdiff_parallel` (redundant with `gcatest::gcat.stat`)

# 2021-04-14 - jackstraw 1.3.4.9000

- Functions `jackstraw_pca`, `jackstraw_rpca`, `jackstraw_irlba`:  Corrected documentation (parameter `r1` was incorrectly described as `PC` in parts of the documentation.  Thanks to Djordje Bajić (GitHub username `djbajic`) for reporting this error!
- Removed option `seed` from all functions that had it.  For the same behavior, call `set.seed(seed)` before calling the function.
- Functions `jackstraw_lfa` and `jackstraw_alstructure`: removed `devR` option.
- Removed redundant functions
  - `lfa.corpcor`: same as `lfa::lfa` with option `override = TRUE`
  - `pi0est_bootstrap`: redundant with `qvalue::pi0est` with option `pi0.method = 'bootstrap'`
  - `dev.R` (internal; functionality implemented in package `gcatest`)
  - `devdiff_parallel` (internal; redundant with `gcatest::gcat.stat`)
  - `getp` (internal; redundant with `qvalue::empPvals`)
