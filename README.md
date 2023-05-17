# FBCModelTests.jl

| Build status | Documentation |
|:---:|:---:|
| ![CI status](https://github.com/LCSB-BioCore/FBCModelTests.jl/workflows/CI/badge.svg?branch=master) [![codecov](https://codecov.io/gh/LCSB-BioCore/FBCModelTests.jl/branch/master/graph/badge.svg?token=NABDEWA380)](https://codecov.io/gh/LCSB-BioCore/FBCModelTests.jl) | [![stable documentation](https://img.shields.io/badge/docs-stable-blue)](https://lcsb-biocore.github.io/FBCModelTests.jl/) [![dev documentation](https://img.shields.io/badge/docs-dev-cyan)](https://lcsb-biocore.github.io/FBCModelTests.jl/dev) |

A collection of tests for constraint-based metabolic models.

There are currently 2 main test suites implemented:

- [FROG reproducibility and curation checks](https://www.ebi.ac.uk/biomodels/curation/fbc)
- [MEMOTE](https://memote.readthedocs.io/)-style model consistency and annotation checks

## FROG

You can generate FROG reports using function `FROG.generate_report` and compare
the "compatibility" of two reports using `FROG.compare_reports`. See the
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

### Running FROG from command line

You can use the supplied scripts to conveniently run FROG from a command line on
a system that has FBCModelTests.jl (and possibly a solver) already installed.

After copying the files from the `bin/` directory in this repository into your
`$PATH`, you can use them to generate and compare FROG reports as follows:
```sh
$ fbcmt-run-frog -s GLPK model.xml report_dir
$ fbcmt-compare-frog report_dir other_dir
```

### Running FROG with Docker

A pre-packaged dockerized version of the commands is available from GHCR. The
following commands run the dockerized versions of the above scripts:
```
$ docker run -ti --rm -v $PWD:/data -w /data ghcr.io/lcsb-biocore/docker/fbcmodeltests-run-frog -s GLPK model.xml report_dir
$ docker run -ti --rm -v $PWD:/data -w /data ghcr.io/lcsb-biocore/docker/fbcmodeltests-compare-frog report_dir other_dir
```

## MEMOTE-style tests

The primary entry point for the [MEMOTE](https://memote.readthedocs.io/) test
suit implemented here is the function `run_tests`. When building a model, it is
most convenient to incorporate it into the CI of the model. Another option is to
use the command line functionality, and save the output for later analysis.

To run the test suite on a toy model, use `run_tests`:
```julia
using FBCModelTests, GLPK, Distributed
addprocs(10)
FBCModelTests.Memote.run_tests("e_coli_core.json", GLPK.Optimizer; workers=workers())
```
Any optimizer supported by [JuMP](https://jump.dev/) can be used. The output of
`run_tests` is the standard Julia unit testing scheme. However, in the repl the
full output is usually truncated, and only a summary is shown. If you want more
details about where/why your model failed certain tests, it is best to capture
the output, and save it to a file. A convenient way to do this is with
[ansi2html](https://github.com/agnoster/ansi2html). Additionally, to make the
output more display friendly, we recommend `run_tests_toplevel` is used instead
of `run_tests`.

An example workflow entails using the scripts located in `bin/`:
```
fbcmt-memote-run --color=yes -s GLPK -w 6 e_coli_core.xml > e_coli_core.test.out
ansi2html < e_coli_core.test.out > e_coli_core.test.html
```
The resultant `html` can be inspected in any browser.

See the function documentation for additional test configuration information.
Note, the tests implemented here are significantly more conservative than in the
original Memote. In particular, no heuristics are used to guess reaction types,
e.g. biomass, atp maintenance, transporter, exchange, etc. Only [SBO
annotations](https://github.com/EBI-BioModels/SBO/blob/master/SBO_OBO.obo) are
used for this purpose, because only these are actually standardized.
Consequently, all tests that rely on properly annotated reactions will fail if
this is not incorporated into the model being tested.

The implementation in FBCModelTests.jl is mostly authored by
St. Elmo Wilken ([@stelmo](https://github.com/stelmo))
with parts contributed by
Mirek Kratochvíl ([@exaexa](https://github.com/exaexa)),
Vincent M. von Häfen ([@vm-vh](https://github.com/vm-vh)),
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
The development was supported by European Union's Horizon 2020 Programme under
PerMedCoE project ([permedcoe.eu](https://permedcoe.eu/)) agreement no. 951773.

<img src="https://lcsb-biocore.github.io/FBCModelTests.jl/dev/assets/unilu.svg" alt="Uni.lu logo" height="64px" style="height:64px; width:auto">   <img src="https://lcsb-biocore.github.io/FBCModelTests.jl/dev/assets/lcsb.svg" alt="LCSB logo" height="64px" style="height:64px; width:auto">   <img src="https://lcsb-biocore.github.io/FBCModelTests.jl/dev/assets/hhu.svg" alt="HHU logo" height="64px" style="height:64px; width:auto">   <img src="https://lcsb-biocore.github.io/FBCModelTests.jl/dev/assets/qtb.svg" alt="QTB logo" height="64px" style="height:64px; width:auto">   <img src="https://lcsb-biocore.github.io/FBCModelTests.jl/dev/assets/permedcoe.svg" alt="PerMedCoE logo" height="64px" style="height:64px; width:auto">
