# How do I skip the first line?
awk -F, '{print ">" $2 "_" $3 "_F\n" $4 "\n" ">" $2 "_" $3 "_R\n" $5}' tma-orf-primes.csv > foo