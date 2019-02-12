
plotbio <- function(protein) {
  library(bio3d)
  s <- read.pdb(protein)
  s.chainA <- trim.pdb(s, chain = "A", elety = "CA")
  s.b <- s.chainA$atom$b
  plotb3(s.b, sse = s.chainA,type = 'l', ylab = "Bfactor")
}

