#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=10Gb
#SBATCH --array=1
#SBATCH --partition=main
#SBATCH --mail-type=all
#SBATCH --mail-user=zzhang39@usc.edu

#create annotation files for union
mkdir unioncells
R
for (chr in 1:22){
    unioncells = read.table(paste0("B_IN","/cs.",chr,".annot.gz"),h=T)[,1]
    for (CT in c("B_Mem", "CD4_NC", "CD8_ET", "CD8_NC", "Mono_C", "Mono_NC", "NK",
                  "CD4_ET", "CD4_SOX4", "CD8_S100B", "DC", "NK_R", "Plasma")){
        unioncells = unioncells + read.table(paste0(CT,"/cs.",chr,".annot.gz"),h=T)[,1]
    }
    unioncells = as.integer(unioncells>0)
    write.table(unioncells,file=paste0("unioncells/cs.",chr,".annot"),row.names=F,col.names=T,sep="\t",quote=F)
}
q()
n
gzip unioncells/cs.*.annot

# check should be equal 1012
ll */*annot.gz | wc -l
