# concatanate multiple files

# The file where everything will be merged
output="merged_file.csv"

# Make sure the output file is empty
> "$output"

# Loop through all .txt files in the current directory
for file in *.PLS.ENCODE; do
  # Check if the file is not empty
  if [ -s "$file" ]; then
    # Add a new column with the filename to each line and append to the merged file
    awk -v fname="$file" 'BEGIN{OFS="\t"} {print $0, fname}' "$file" >> "$output"
  fi
done

R
library(tidyverse)
cols <- c("chr", "phenotype_id", "GeneSymbol", "tss_left", "tss_right", "num_var", "hits", "ref_cell")
df=read_tsv("merged_file.csv",F)

df <- df %>% mutate(X8 = gsub(".subtract.PLS.ENCODE", "", X8)) %>% 
separate(X8, into=c("phenotype_id", "ref_cell"), sep="_") %>% select(-phenotype_id)
colnames(df) <- cols

df %>% write_tsv("egenes_bulk.overlap.subtract.PLS.new.tsv.gz")
