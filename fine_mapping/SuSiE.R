## fine-mapping

library(susieR)
set.seed(1)

data(N3finemapping)
attach(N3finemapping)
n = nrow(X)
n

b <- true_coef[,1]
plot(b, pch=16, ylab='effect size')

which(b != 0)

sumstats <- univariate_regression(X, Y[,1])
z_scores <- sumstats$betahat / sumstats$sebetahat
susie_plot(z_scores, y = "z", b=b)

# number of individuals on rows
R <- cor(X)

# by default, assume L=10
fitted_rss2 = susie_rss(z = z_scores, R = R, n = n, L = 100,
                        estimate_residual_variance = TRUE)

# marginal PIP for each variant
fitted_rss2$pip
fitted_rss2$converged
dim(fitted_rss2$mu)
dim(fitted_rss2$alpha)
rowSums(fitted_rss2$alpha)

cs_df = summary(fitted_rss2)$cs
cs_df$cs_avg_r2
cs_df$cs_min_r2

colSums(fitted_rss2$alpha * fitted_rss2$mu)

susie_get_cs(fitted_rss2, coverage = 0.95)
susie_get_pip(fitted_rss2)

susie_plot(fitted_rss2, y="PIP", b=b, add_bar=TRUE)
susie_plot(fitted_rss2, y="PIP", main = "gene")

