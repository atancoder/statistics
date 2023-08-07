# install.packages("vcd")

# Q 2.1
vals <- rpois(100, lambda = .05)
barplot(table(vals))
library("vcd")
gf1 <- goodfit(vals, "poisson")
rootogram(gf1, xlab = "")

# Q 2.2
table(rpois(100, .6))

# Q 2.3
prod(dpois(c(0, 1, 2, 7), lambda = 0)^(c(58, 34, 7, 1)))
prod(dpois(c(0, 1, 2, 7), lambda = 1)^(c(58, 34, 7, 1)))
prod(dpois(c(0, 1, 2, 7), lambda = .4)^(c(58, 34, 7, 1)))

e100 <- c(rep(0, 58), rep(1, 34), rep(2, 7), rep(7, 1))
mean(e100)


# Task
lambda <- .55
# Determine p-val of getting a 7 with our results
pois_prob_over_n_experiments <- function(m, n = 100, lambda = .5) {
    # across n experiments, find probability that the maximum
    # of the n experiments is >= m
    # Pr(max >= m) = 1 - Pr(result of all n < m)
    prob_less_than_m <- ppois(m - 1, lambda)
    prob_all_less_than_m <- prob_less_than_m^n
    return(1 - prob_all_less_than_m)
}
pois_prob_over_n_experiments(m = 7, n = 100, lambda = .55)
# still a very low probability


library("Biostrings")
staph <- readDNAStringSet("bios221/data/staphsequence.ffn.txt", "fasta")

# Q 2.9
actual_freq <- letterFrequency(staph[[1]], letters = "ACGT", OR = 0)
# What's the probability that each nucleotide occurs with the same probability?
# H0: All nucleotides occur with the same probability (25%)
# HA: There's a different probability distribution for the nucleotides
expected_freq <- length(staph[[1]]) / 4


# We're going to use chi-squared test for goodness of fit

pvec <- rep(.25, 4)
expected_freq <- pvec * length(staph[[1]])
chi_score <- sum((actual_freq - expected_freq)^2 / expected_freq)
p_val <- pchisq(chi_score, 3, lower.tail = FALSE)
p_val # Very small


# To see if the first 10 genes have the same nucleotide distribution,
# use chi-squared test of homogeneity
letter_freq <- sapply(staph, letterFrequency, letters = "ACGT", OR = 0)
colnames(letter_freq) <- paste0("gene", seq(along = staph))
first_10_gene_letter_freq <- letter_freq[, 1:10]
letter_sum <- rowSums(first_10_gene_letter_freq)
letter_prob_distribution <- letter_sum / sum(letter_sum)

expected_freqs <- apply(first_10_gene_letter_freq, function(x) {
    outer(sum(x), letter_prob_distribution)
}, MARGIN = 2) # MARGIN=2 means apply vertically


chi_score <- sum(
    (first_10_gene_letter_freq - expected_freqs)^2 / expected_freqs
)
degrees_of_freedom <- (10 - 1) * (4 - 1)
p_val <- pchisq(chi_score, degrees_of_freedom, lower.tail = FALSE)
# very low p-val, which means unlikely NULL hypothesis
# is true. So each gene comes from a diff distribution

# 2.10
random_chisq <- rchisq(1000, 30)
hist(random_chisq, breaks = 50)
quantile(random_chisq)

# 2.11
# Median

# Quantile-Quantile plot to compare the simulation
# with the chi-sq distribution
B <- 1000
qqplot(qchisq(ppoints(B), df = 30), simulstat,
    main = "",
    xlab = expression(chi[nu == 30]^2), asp = 1, cex = 0.5, pch = 16
)
abline(a = 0, b = 1, col = "red")

# Chargaff's Rule
load("bios221/data/ChargaffTable.RData")

# Q 2.13
# The different organisms don't look like they
# come from the same distribution for letter frequency

# Q 2.15
# It's a table of 4X4 dimension

load("bios221/data/Deuteranopia.RData")
chisq.test(Deuteranopia)

# Q 2.18
# prior
rp <- rbeta(100000, 50, 350)
y <- rbinom(length(rp), rp, size = 300)
hist(y, breaks = 50, col = "orange", main = "", xlab = "")

# what probabilities led to getting y=40 (our data had y=40)?
pPostEmp <- rp[y == 40]
hist(pPostEmp,
    breaks = 40, col = "chartreuse4", main = "",
    probability = TRUE, xlab = "posterior p"
)

p_seq <- seq(0, 1, by = 0.001)
posterior_dist <- dbeta(p_seq, 50 + 40, 350 + 260)
lines(p_seq, posterior_dist, type = "l", lwd = 3)
# mean
mean(pPostEmp)


# MAP (max a posteriori)
p_seq[which.max(posterior_dist)]

# Q 2.20
# see how much prior affects final result
rp <- rbeta(100000, 1, 1)
y <- rbinom(length(rp), rp, size = 300)
hist(y, breaks = 50, col = "orange", main = "", xlab = "")

pPostEmp <- rp[y == 40]
hist(pPostEmp,
    breaks = 40, col = "chartreuse4", main = "",
    probability = TRUE, xlab = "posterior p"
)
p_seq <- seq(0, 1, by = 0.001)
posterior_dist <- dbeta(p_seq, 1 + 40, 1 + 260)
lines(p_seq, posterior_dist, type = "l", lwd = 3)

# mean
mean(pPostEmp)

# MAP
p_seq[which.max(posterior_dist)]

# looks like prior doesn't affect final result as much. Since the prior
# has small data (1,1), it was easily overshadowed with the new observed data
