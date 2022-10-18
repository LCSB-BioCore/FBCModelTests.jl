
# FBCModelTests.jl functions

```@autodocs
Modules = [FBCModelTests]
Pages = ["FBCModelTests.jl"]
```

# FROG

```@autodocs
Modules = [FBCModelTests, FBCModelTests.FROG]
Pages = ["frog.jl", "frog/structs.jl", "frog/frontend.jl"]
```

## Reading and writing reports
```@autodocs
Modules = [FBCModelTests.FROG.ReportIO]
Pages = ["frog/io.jl"]
```

## Generating and testing the reports
```@autodocs
Modules = [FBCModelTests.FROG.ReportGenerators, FBCModelTests.FROG.ReportTests]
Pages = ["frog/report.jl"]
```

# MEMOTE

```@autodocs
Modules = [FBCModelTests.Memote]
Pages = ["memote.jl", "memote/frontend.jl"]
```

## MEMOTE utilities

```@autodocs
Modules = [FBCModelTests.Memote.Utils, FBCModelTests.Memote.Config]
Pages = ["memote/utils.jl", "memote/config.jl"]
```

## MEMOTE checks

### Basic checks

```@autodocs
Modules = [FBCModelTests.Memote.Basic]
Pages = ["src/memote/checks/Basic.jl"]
```

### Consistency checks

```@autodocs
Modules = [FBCModelTests.Memote.Consistency]
Pages = ["src/memote/checks/Consistency.jl"]
```

### Metabolite checks

```@autodocs
Modules = [FBCModelTests.Memote.Metabolite]
Pages = ["src/memote/checks/Metabolite.jl"]
```

### Reaction checks

```@autodocs
Modules = [FBCModelTests.Memote.Reaction]
Pages = ["src/memote/checks/Reaction.jl"]
```

### GPRAssociation checks

```@autodocs
Modules = [FBCModelTests.Memote.GPRAssociation]
Pages = ["src/memote/checks/GPRAssociation.jl"]
```

### Biomass checks

```@autodocs
Modules = [FBCModelTests.Memote.Biomass]
Pages = ["src/memote/checks/Biomass.jl"]
```

### Network checks

```@autodocs
Modules = [FBCModelTests.Memote.Network]
Pages = ["src/memote/checks/Network.jl"]
```

### Annotation checks

```@autodocs
Modules = [FBCModelTests.Memote.Annotation]
Pages = ["src/memote/checks/Annotation.jl"]
```

### Energy checks

```@autodocs
Modules = [FBCModelTests.Memote.Energy]
Pages = ["src/memote/checks/Energy.jl"]
```

# Types and utilities

```@autodocs
Modules = [FBCModelTests]
Pages = ["common.jl"]
```
