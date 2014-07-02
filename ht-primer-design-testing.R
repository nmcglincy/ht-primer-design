# Testing ht-primer design using callPrimer3 function
# 
# First I need a sequence
library("Biostrings")
?readDNAStringSet
ydr114c.seq = readDNAStringSet(filepath = "ydr117c.fasta")
ydr114c.seq
?DNAString
ydr114c.seqA = DNAString(ydr114c.seq[[1]])
ydr114c.seqA

# # Load the function
source("/Users/nmcglincy/Documents/computing/github/ingolia_lab/ht-primer-design/callPrimer3-NJM.R")
ls()

foo = p3.output(seq = ydr114c.seqA, name = "ydr117c")
f00 = .callP3NreadOrg(seq = ydr114c.seqA, name = "ydr117c")
ls()
f00

.callP3NreadOrg(seq = ydr114c.seqA, name = "ydr117c")
