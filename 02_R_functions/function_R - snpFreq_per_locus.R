#function for checking snp freq per column/locus
Func_SnpFreq <- function(x) {
  cbind(freq = table(x), 
    percentage = prop.table(table(x))*100)
}