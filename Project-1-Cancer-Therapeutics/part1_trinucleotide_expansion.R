library(Biostrings)
library(tidyverse)

##reading in input files
fasta_path <- "BNFOgenes.fasta"
csv_path <- "TRE_Cancer_Status.csv"

##reading in DNA sequence(fasta) into DNAStringSet
fasta_seq <- Biostrings::readDNAStringSet(fasta_path)

#reading in cancer csv
status_tbl <-readr::read_csv(csv_path, show_col_types = FALSE)

# headers for donor_id
donor_id <- sub("(Donor\\s+\\d+).*", "\\1", names(fasta_seq))

##inspecting in the console
length(fasta_seq)
head(donor_id)
fasta_seq[1:20]



##--------------------------------------------------
## counting 'CTG' trinucleotide occurence per donor
##--------------------------------------------------

##counting the exact 'CTG' matches in each sequence 
CTG_count <- Biostrings::vcountPattern(pattern = "CTG", subject = fasta_seq, fixed = TRUE)

##table to show results
counts_tbl <- data.frame(Donor = donor_id, CTG_count = as.integer(CTG_count))

##looking at it in the console
head(counts_tbl)
summary(counts_tbl$CTG_count)



##--------------------------------------------------
##computing the cutoff C1 at 0.01
##--------------------------------------------------
lambda0 <- 33
##upper cutoff C1 at alpha = 0.01
C1 <- qpois(p = 0.99, lambda = lambda0)
##cutoff value is 47!!!!



##calculating the p-value for each donor
##p-value: P(X >= observed | Poisson(33))
p_val <- ppois(q = counts_tbl$CTG_count - 1, lambda = lambda0, lower.tail = FALSE)
##adding p-values to my table and making the expanded 
results <- transform(counts_tbl, p_val = p_val, Expanded = p_val < 0.01)


##association with cancer

## 1. Joining expansion calls to cancer status by donor
analysis_tbl <- dplyr::inner_join(results, status_tbl, by = "Donor")
## 2. 2x2 table representing expanded vs cancer status
xtab <- table(Expanded = analysis_tbl$Expanded, Cancer = analysis_tbl$Cancer_Status )
## 3. using fisher's exact test
f_out <- fisher.test(xtab)

### Fisher's Exact Test for Count Data

## data:  xtab
## p-value = 0.06978
## alternative hypothesis: true odds ratio is not equal to 1
## 95 percent confidence interval: 0.007870555 1.133635839
## sample estimates: odds ratio 0.1226533 
