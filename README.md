# FBCModelTests.jl

| Build status | Documentation |
|:---:|:---:|
| ![CI status](https://github.com/LCSB-BioCore/FBCModelTests.jl/workflows/CI/badge.svg?branch=master) | [![stable documentation](https://img.shields.io/badge/docs-stable-blue)](https://lcsb-biocore.github.io/FBCModelTests.jl/) [![dev documentation](https://img.shields.io/badge/docs-dev-cyan)](https://lcsb-biocore.github.io/FBCModelTests.jl/dev) |

A collection of tests for constraint-based metabolic models.

There are currently 2 main test suites implemented:

- [FROG reproducibility and curation checks](https://www.ebi.ac.uk/biomodels/curation/fbc)
- [Memote](https://memote.readthedocs.io/)-style model consistency and annotation checks

## FROG

You can generate FROG reports using function `frog_generate_report` and compare
the "compatibility" of two reports using `frog_compare_reports`. See the inline
documentation for more details.

The supposed workflow with FROG reports is the following:
- If you are a model author, you generate a FROG report that is used as a
  reference solution of your metabolic model problem, and distribute it along a
  model.
- If you are a model curator, you demand the FROG report from the original
  author as a "certificate" of their model's intended functionality.
- If you use a model of someone other, you generate another FROG report with
  your analysis software, and then compare the report against the original
  model author's report to verify that your software is interpreting the model
  information correctly and that the solutions are compatible within a
  reasonable tolerance.

The implementation is based off the description of the [EBI's FBC curation
site](https://www.ebi.ac.uk/biomodels/curation/fbc), with some details
following the decisions in the [`fbc_curation` Python
tool](https://github.com/matthiaskoenig/fbc_curation) (working with COBRApy) by
Matthias König (@matthiaskoenig).

The implementation in FBCModelTests.jl is mostly authored by
Mirek Kratochvíl (@exaexa)
with parts contributed by
St Elmo Wilken (@stelmo).

## Memote-style tests

This part is currently under construction.

# Acknowledgements

FBCModelTests.jl package is developed at the
Luxembourg Centre for Systems Biomedicine of the University of Luxembourg
([uni.lu/lcsb](https://wwwen.uni.lu/lcsb))
and the
Institute for Quantitative and Theoretical Biology
at the
Heinrich Heine University in Düsseldorf
([qtb.hhu.de](https://www.qtb.hhu.de/)).
