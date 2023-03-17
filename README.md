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

You can use the supplied scripts to conveniently run FROG from a commandline on
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
$ docker run -u $UID -ti --rm -v $PWD:/data -w /data ghcr.io/lcsb-biocore/docker/fbcmodeltests-run-frog -s GLPK model.xml report_dir
$ docker run -ti --rm -v $PWD:/data -w /data ghcr.io/lcsb-biocore/docker/fbcmodeltests-compare-frog report_dir other_dir
```

## MEMOTE-style tests

You can use a number of tests that automatically check various basic quality
characteristics of a constraint-based metabolic model; the suite available in
FBCModelTests.jl is inspired by [MEMOTE](https://memote.readthedocs.io/).

To run the test suite on a toy model, use `run_tests`:
```julia
using FBCModelTests, Tulip
FBCModelTests.Memote.run_tests("e_coli_core.json", Tulip.Optimizer)
```

An overview of the model properties can be also gathered in a `Dict` for
mechanical inspection (or saving into JSON):
```julia
using FBCModelTests, Tulip
structured_report = FBCModelTests.Memote.generate_report("e_coli_core.json, Tulip.Optimizer)

using JSON
open("my_e_coli_core_report.json", "w") do io
    JSON.print(io, report)
end
```

See the function documentation for additional test configuration, and several
performance gotchas and notes about MEMOTE compability.

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
The development was supported by European Union's Horizon 2020 Programme under
PerMedCoE project ([permedcoe.eu](https://permedcoe.eu/)) agreement no. 951773.

<img src="https://lcsb-biocore.github.io/FBCModelTests.jl/dev/assets/unilu.svg" alt="Uni.lu logo" height="64px" style="height:64px; width:auto">   <img src="https://lcsb-biocore.github.io/FBCModelTests.jl/dev/assets/lcsb.svg" alt="LCSB logo" height="64px" style="height:64px; width:auto">   <img src="https://lcsb-biocore.github.io/FBCModelTests.jl/dev/assets/hhu.svg" alt="HHU logo" height="64px" style="height:64px; width:auto">   <img src="https://lcsb-biocore.github.io/FBCModelTests.jl/dev/assets/qtb.svg" alt="QTB logo" height="64px" style="height:64px; width:auto">   <img src="https://lcsb-biocore.github.io/FBCModelTests.jl/dev/assets/permedcoe.svg" alt="PerMedCoE logo" height="64px" style="height:64px; width:auto">
