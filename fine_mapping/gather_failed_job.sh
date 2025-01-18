# extract file name for unfinished LD calculation

celltype="CD4_SOX4"
touch ${celltype}.failed

# Iterate through all files in the current directory
for file in *.log; do
    # Check if the file contains the pattern. If not, process the filename.
    if ! grep -q "write out" "$file"; then
        # Extract the prefix of the filename. Adjust this line according to your definition of 'prefix'.
        prefix="${file%%.*}"
        
        # Write the prefix to the output file
        echo "$prefix" >> ${celltype}.failed
    fi
done
