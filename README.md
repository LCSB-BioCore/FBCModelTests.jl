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
tool](https://github.com/matthiaskoenig/fbc_curation) (working with [COBRApy](https://github.com/opencobra/cobrapy/)) by
Matthias König ([@matthiaskoenig](https://github.com/matthiaskoenig)).

The implementation in FBCModelTests.jl is mostly authored by
Mirek Kratochvíl ([@exaexa](https://github.com/exaexa))
with parts contributed by
St. Elmo Wilken ([@stelmo](https://github.com/stelmo)).

## Memote-style tests

This package exposes a number of tests aimed at quickly checking if certain
basic quality characteristics of a constraint-based metabolic model are
satisfied. Assuming your model will work with the default configuration
arguments of the test functions (see their docstrings), then you can test your
model with:

```
using FBCModelTests
using COBREXA, Tulip, Test

model = load_model("e_coli_core.json")

@testset "Test model" begin
  run_tests(model, Tulip.Optimizer; config = memote_config)
end

model_info = generate_memote_report(model, optimizer; config = memote_config)
```
You can set the configuration parameters through adjusting `memote_config`,
which is a type of `MemoteConfig`. Note, some of the original Memote tests
involve solving MILPs, these tests are not included in this package as the
issues they identify are often easier to identify by simpler means. For example, find the stoichiometrically inconsistent metabolites is simpler done by just checking the mass balance around each reaction.

The implementation in FBCModelTests.jl is mostly authored by
St. Elmo Wilken ([@stelmo](https://github.com/stelmo))
with parts contributed by
Vincent M. von Häfen ([@vm-vh](https://github.com/vm-vh))
and Flora Schlüter ([@Fl-Sch](https://github.com/Fl-Sch)).

# Acknowledgements

FBCModelTests.jl package is developed at the
Luxembourg Centre for Systems Biomedicine of the University of Luxembourg
([uni.lu/lcsb](https://wwwen.uni.lu/lcsb))
and the
Institute for Quantitative and Theoretical Biology
at the
Heinrich Heine University in Düsseldorf
([qtb.hhu.de](https://www.qtb.hhu.de/)).
