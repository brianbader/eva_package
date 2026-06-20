# eva 0.2.7

## Bug fixes

* `gpdFit()`: fixed parameter initialization when starting values are supplied
  via the `start` argument; `start` is now correctly split into scale and shape
  components instead of being passed through unprocessed (#6).
* Sequential tests: fixed an erroneous double `rev()` in the ForwardStop /
  StrongStop computation of `gpdSeqTests()` and `gevrSeqTests()`. As a result,
  `pSeqStop()` now takes p-values in the order the hypotheses are tested and
  preserves that order (previously it expected pre-reversed input).
* `gevrFit()`: the variance-covariance matrix is now computed safely; a singular
  Hessian returns an `NA` covariance matrix instead of throwing an error.

## Other changes

* `gevrFit()`: added `extendInt = "yes"` to the internal `uniroot()` call to make
  the moment-based starting value for the shape parameter more robust.
* `gpdSeqTests()`: clarified documentation that candidate `thresholds` are tested
  in the order supplied and should be given from lowest to highest, with an
  example of applying the ForwardStop rule.
* Added formal package citations for the *Statistics and Computing* (2017) and
  *Annals of Applied Statistics* (2018) papers.
* Vignette now uses `SpatialExtremes` conditionally (guarded by
  `requireNamespace()`) so it builds even when that package is not installed.
