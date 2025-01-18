# gather sldsc result
method="jaxqtl"

# The output file where everything will be merged
# output="${method}_sldsc_traits107_union.label.tsv"
# output="${method}_sldsc_traits107_union.label.newselection.tsv"
output="${method}_sldsc_traits107_union.label.addunion.tsv"

# Make sure the output file is empty
> "$output"

# Loop through all .txt files in the current directory
for file in *.results; do
  # Check if the file is not empty
  if [ -s "$file" ]; then
    # Add a new column with the filename to each line using a tab delimiter and append to the merged file
    awk -v fname="$file" 'BEGIN{OFS="\t"} {print $0, fname}' "$file" >> "$output"
  fi
done

gzip $output



