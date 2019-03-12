Class11: Structural Bioinformatics 1
================

Section 1. The PDB
------------------

Q1: Determine the percentage of structures solved by X-Ray and Electron Microscopy. From the website what proportion of structures are protein? Aim to have a rendered GitHub document with working code that yields your answers.

``` r
stats <- read.csv("Data Export Summary.csv", row.names = 1)
stats
```

    ##                     Proteins Nucleic.Acids Protein.NA.Complex Other  Total
    ## X-Ray                 124770          1993               6451    10 133224
    ## NMR                    10988          1273                257     8  12526
    ## Electron Microscopy     2057            31                723     0   2811
    ## Other                    250             4                  6    13    273
    ## Multi Method             127             5                  2     1    135

``` r
percent.by.method <- stats$Total/sum(stats$Total)*100
names(percent.by.method) <- rownames(stats)
percent.by.method
```

    ##               X-Ray                 NMR Electron Microscopy 
    ##         89.43068692          8.40846082          1.88696977 
    ##               Other        Multi Method 
    ##          0.18325960          0.09062288

``` r
sum(stats$Proteins)/sum(stats$Total)*100
```

    ## [1] 92.76561

Section 3. Introduction to Bio3D in R
-------------------------------------

First we have to load the package we want to use

``` r
library(bio3d)
```

Now we can use the functions from this bio3d package

``` r
pdb <- read.pdb("1hsg")
```

    ##   Note: Accessing on-line PDB file

``` r
pdb
```

    ## 
    ##  Call:  read.pdb(file = "1hsg")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

``` r
#head(aa321(pdb$atom$resid))
```

``` r
ca.inds <- atom.select(pdb, "calpha")
ca.inds
```

    ## 
    ##  Call:  atom.select.pdb(pdb = pdb, string = "calpha")
    ## 
    ##    Atom Indices#: 198  ($atom)
    ##    XYZ  Indices#: 594  ($xyz)
    ## 
    ## + attr: atom, xyz, call

``` r
# Print details of the first few selected atoms
head( pdb$atom[ca.inds$atom, ] )
```

    ##    type eleno elety  alt resid chain resno insert      x      y     z o
    ## 2  ATOM     2    CA <NA>   PRO     A     1   <NA> 30.307 38.663 5.319 1
    ## 9  ATOM     9    CA <NA>   GLN     A     2   <NA> 30.158 36.492 2.199 1
    ## 18 ATOM    18    CA <NA>   ILE     A     3   <NA> 29.123 33.098 3.397 1
    ## 26 ATOM    26    CA <NA>   THR     A     4   <NA> 29.774 30.143 1.062 1
    ## 33 ATOM    33    CA <NA>   LEU     A     5   <NA> 27.644 27.003 1.144 1
    ## 41 ATOM    41    CA <NA>   TRP     A     6   <NA> 30.177 24.150 1.279 1
    ##        b segid elesy charge
    ## 2  40.62  <NA>     C   <NA>
    ## 9  41.30  <NA>     C   <NA>
    ## 18 34.13  <NA>     C   <NA>
    ## 26 30.14  <NA>     C   <NA>
    ## 33 30.12  <NA>     C   <NA>
    ## 41 30.82  <NA>     C   <NA>

``` r
# And selected xyz coordinates
head( pdb$xyz[, ca.inds$xyz] )
```

    ## [1] 30.307 38.663  5.319 30.158 36.492  2.199

``` r
## [1] 64.044 51.615 33.423 62.439 54.794 32.359
```

In addition to the common selection strings (such as ‘calpha’ ‘cbeta’ ‘backbone’ ‘protein’ ‘notprotein’ ‘ligand’ ‘water’ ‘notwater’ ‘h’ and ‘noh’) various individual atom properties can be used for selection.

``` r
# Select chain A
a.inds <- atom.select(pdb, chain="A")
# Select C-alphas of chain A
ca.inds <- atom.select(pdb, "calpha", chain="A")
# We can combine multiple selection criteria to return their intersection
cab.inds <- atom.select(pdb, elety=c("CA","CB"), chain="A", resno=10:20)
cab.inds
```

    ## 
    ##  Call:  atom.select.pdb(pdb = pdb, elety = c("CA", "CB"), chain = "A", 
    ##     resno = 10:20)
    ## 
    ##    Atom Indices#: 20  ($atom)
    ##    XYZ  Indices#: 60  ($xyz)
    ## 
    ## + attr: atom, xyz, call

We want to select out the protein and drug only parts of these molecular files

``` r
prot.inds <- atom.select(pdb, "protein")
prot.inds
```

    ## 
    ##  Call:  atom.select.pdb(pdb = pdb, string = "protein")
    ## 
    ##    Atom Indices#: 1514  ($atom)
    ##    XYZ  Indices#: 4542  ($xyz)
    ## 
    ## + attr: atom, xyz, call

``` r
prot.pdb <- trim.pdb(pdb, prot.inds)
write.pdb(prot.pdb, file = "protein.pdb")
```

Let's do the same things for our drug

``` r
lig.inds <- atom.select(pdb, "ligand")
lig.pdb <- trim.pdb(pdb, lig.inds)
write.pdb(lig.pdb, file = "ligand.pdb")
```
