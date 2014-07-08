# How do I skip the first line?
awk -F, '{print ">" $1 "_" $2 "_F\n" $3 "\n" ">" $1 "_" $2 "_R\n" $4}' here-are-your-primers.csv > foo