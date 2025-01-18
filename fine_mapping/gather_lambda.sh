# concatanate multiple files

# The file where everything will be merged
method="wald"
output="NB_lambda.${method}.tsv"

# Make sure the output file is empty
> "$output"

# Loop through all .txt files in the current directory
for file in *.lambda.tsv; do
  # Check if the file is not empty
  if [ -s "$file" ]; then
    # Add a new column with the filename to each line and append to the merged file
    tail -n +2 "$file" | awk -v fname="$file" 'BEGIN{OFS="\t"} {print $0, fname}' "$file" >> "$output"
  fi
done