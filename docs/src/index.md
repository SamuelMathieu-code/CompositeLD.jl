# CompositeLD.jl Documentation


## Overview 

CompositeLD is a package using [SnpArrays](https://github.com/OpenMendel/SnpArrays.jl) to perform some usefull calculations related LD. In particular :
- Calculate LD between two SNPs (given position or id in PLINK files)
- Get the LD matrix of a group of SNPs
- Get Strongly correlated SNPs to a given group of SNPs
- Clumping given SNPs in a prioritized order.

## Example

```julia
using SnpArrays
using CompositeLD

data = SnpData(SnpArrays.datadir(datapath))

ld = ld_r2("rs1234", "rs654321", data)

ld_mat = getLDmat(data, [(1, 459876), (1, 58735), (2, 97654)])

strong = getStrongLD(data, ["rs1234", "rs5678"])

kept_v_b = clump(data, [(1, 459876), (1, 58735), (2, 97654)])
```

## Contents

```@contents
```

## Functions

```@docs
formatSnpData!

ld_r2

getLDmat

getStrongLD

clump

tclump
```

## Index

```@index
```