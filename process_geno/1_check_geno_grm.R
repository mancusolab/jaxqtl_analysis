# check genotype

library(tidyverse)
library(genio)
grm <- read_grm("../result/check_geno/allchr")
ltriag_grm <- grm$kinship[lower.tri(grm$kinship, diag = FALSE)]

summary(ltriag_grm)
length(ltriag_grm)
length(ltriag_grm[ltriag_grm < max(ltriag_grm)])

tibble(grm_r=ltriag_grm) %>% ggplot(aes(x=ltriag_grm)) + geom_histogram()+theme_bw()
tibble(grm_r=ltriag_grm) %>% ggplot(aes(y=ltriag_grm, x="grm")) + geom_boxplot()+theme_bw()

# check extreme value for each row (person): 958_959 (row 867), 960_961 (row 869)
which(grm$kinship==max(ltriag_grm), arr.ind=TRUE)

which(grm$kinship>0.075, arr.ind=TRUE) %>% as_tibble() %>%  filter(row != col) %>% View

# 958_959
grm$kinship[869,][grm$kinship[869,] == max(ltriag_grm)]

summary(grm$kinship['958_959',][grm$kinship['958_959',] < 0.11])
summary(grm$kinship['960_961',][grm$kinship['960_961',] < 1.])

heatmap(grm$kinship, Rowv = NA, Colv = NA, scale="none")  
