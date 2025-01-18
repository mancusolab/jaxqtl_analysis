# colocalization
library(coloc)
library(susieR)

# hypothesis
# H0: no association with either trait in the region
# H1: association with trait 1 only
# H2: association with trait 2 only
# H3: both traits are associated, but have different single causal variants
# H4: both traits are associated and share the same single causal variant


### data structure ### 
# SNPs must * have summary data available in both studies 
# Please do include all SNPs you have in the region: do not trim by significance or MAF or other means. 

data(coloc_test_data)
attach(coloc_test_data)

str(D1)

minimum_data=D1[c("beta","varbeta","snp","position","type","sdY")]
str(minimum_data)

check_dataset(minimum_data,warn.minp=1e-10)

# -log10(p) check SNPs are dense in this region
plot_dataset(minimum_data)

# convert df to coloc data structure
as.list(df)

# create LD data structure (add dimension name)
str(D1$LD)

# check alignment (should be more positive values than negative values)
prop_pos <- check_alignment(D1) # proportion of pairs that are positive, eg. 0.9

# plot two datasets
str(D3)
str(D4)
par(mfrow=c(2,1))
plot_dataset(D3, main="Dataset D3")
plot_dataset(D4, main="Dataset D4")

# standard coloc
my.res <- coloc.abf(dataset1=D3, dataset2=D4)

# H4: both traits are associated and share the same single causal variant
sensitivity(my.res,"H4 > 0.9")

# can use same or different LD matrix betwen two traits
check_dataset(D3,req="LD")
check_dataset(D4,req="LD")
str(D3$LD)
str(D4$LD)

# Run coloc using SuSiE (allow multiple causal variants)
# First run susie on two data separately
str(D3)
D3$z <- D3$beta / sqrt(D3$varbeta)
D3$beta <- rep(1, length(D3$beta))
D3$varbeta <- rep(1, length(D3$varbeta))
# D3$sdY <- NULL # if quant outcome doesn't have sdY, or if standardized then set sdY = 1
# D3$type <- 'cc' # change to case-control study

check_dataset(D3)
S3 <- runsusie(D3, maxit = 100, L = 10, 
               estimate_residual_variance=FALSE, 
               refine=TRUE)

summary(S3)

D4$z <- D4$beta / sqrt(D4$varbeta)
D4$beta <- rep(1, length(D4$beta))
D4$varbeta <- rep(1, length(D4$varbeta))

S4 <- runsusie(D4, maxit = 100, L = 10, 
               estimate_residual_variance=FALSE, 
               refine=TRUE)
summary(S4)

susie.res <- coloc.susie(S3,S4)
susie.res$summary

sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=D3,dataset2=D4)
sensitivity(susie.res,"H4 > 0.9",row=2,dataset1=D3,dataset2=D4)
