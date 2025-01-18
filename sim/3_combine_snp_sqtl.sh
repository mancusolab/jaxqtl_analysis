# combine snp info together

# params_Va params_Vre params_beta0 params_type1 params_num_cells

params_file="params_sqtl_100cells"
cd /project/nmancuso_8/elezhang/projects/jaxqtl/code/sim/output/params_Va_Vre_cells_beta0

num_of_lines=$(< "../../params/${params_file}" wc -l)

# combine 
# for sim in {1..${num_of_lines}}; do
for IDX in `seq 1 ${num_of_lines}`; do
# run w/e you need to run
params=`sed "${IDX}q;d" ../../params/${params_file}`
echo "Running instance ${IDX} with params: ${params}"
set -- junk $params
shift

sim=$1
echo $sim
awk FNR!=2 sim1.pheno1_sqtl_causal_snp > header
awk FNR!=1 sim${sim}.pheno*_sqtl_causal_snp > sim${sim}_sqtl.tsv
cat sim${sim}_sqtl.tsv >> header
mv header sim${sim}_sqtl.tsv
done


# remove file with only one line
for filename in ./*_sqtl.tsv; do
    if [ "$( wc -l <"$filename" )" -eq 1 ]; then
        rm -f "$filename"
    fi
done
