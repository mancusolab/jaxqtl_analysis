# create collapse model of gtf

# 10x genomics single cell sequencing is strand specific
# counts sense strand only
# ref: https://kb.10xgenomics.com/hc/en-us/articles/360004396971-Are-10x-Single-Cell-gene-expression-libraries-strand-specific#:~:text=Answer%3A%20Yes%2C%2010x%20Single%20Cell,structure%20as%20the%20sense%20transcripts.

# here using release 84 from ensemble ftp:
# download from: https://ftp.ensembl.org/pub/grch37/release-84/gtf/homo_sapiens/

# note: collapse model gives same output for v82 and v87
python3 collapse_annotation.py \
 Homo_sapiens.GRCh37.82.gtf.gz \
 Homo_sapiens.GRCh37.82.genes.gtf \
 --collapse_only

# create 
python gtf_to_tss.py
