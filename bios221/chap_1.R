# 1.3.3
distribution <- dpois(0:12, 5)
barplot(distribution, names.arg = 0:12, col = "blue")

# task to simulate mutations over 1000 replication cycles
num_mutations <- rbinom(1000, prob = 5 * 10^-4, size = 10^4)
mutation_freq <- table(num_mutations)
barplot(mutation_freq, col = "red")


# 1.4 task

generate_nucleotide <- function() {
    rv <- runif(1)
    if (rv < 1 / 8) {
        return("A")
    } else if (rv < 4 / 8) {
        return("C")
    } else if (rv < 7 / 8) {
        return("G")
    } else {
        return("T")
    }
}
generate_nucleotide()

# q 1.8
pvec <- rep(1 / 4, 4)
t(rmultinom(1, prob = pvec, size = 8))
# t stands for transpose. This overall logic distributes
# 8 balls into 4 containers
rmultinom(n = 8, prob = pvec, size = 1)
# the latter runs 8 separate simulations while the former
# does 1 simulation of throwing 8 balls

# 1.4.1 task
# Goal is to find n such that we can have >=80% power
# 1000 simulations should suffice
calc_stat <- function(observed, expected) {
    return(sum((observed - expected)^2 / expected))
}

calc_power_pct <- function(n, pvec) {
    null_h_pvec <- c(1/4, 1/4, 1/4, 1/4)
    observed_null <- rmultinom(1000, prob = null_h_pvec, size = n)
    stat_results_null <- apply(observed_null, MARGIN = 2, FUN=calc_stat, null_h_pvec*n)
    q95 = quantile(stat_results_null, probs = 0.95)
    print(q95)

    observed <- rmultinom(1000, prob = pvec, size = n)
    stat_results <- apply(observed, MARGIN = 2, FUN=calc_stat, null_h_pvec*n)
    num_greater <- sum(stat_results > q95)
    print(quantile(stat_results, probs=.95))
    power_pct <- 100 * num_greater/length(stat_results) 
    return(power_pct)
}

pvec <- c(3/8, 1/4, 1/4, 1/8)
calc_power_pct(n=90, pvec)

# ANS: n = 90 is sufficient to hit power of 80%


### Exercises
#1.1
pbinom, pnorm, ppois

#1.2
dbinom(2, 10, .3)
# ans: 23.3%

#1.3 & 1.4
pois_prob_over_n_experiments <- function(m, n=100, lambda=.5) {
    # across n experiments, find probability that the maximum
    # of the n experiments is >= m
    # Pr(max >= m) = 1 - Pr(result of all n < m)
    prob_less_than_m <- ppois(m-1, lambda)
    prob_all_less_than_m <- prob_less_than_m^n
    return(1-prob_all_less_than_m)
}
pois_prob_over_n_experiments(lambda = 0.5, n = 100, m = 7)

# 1.5
simulate <- function(num_simulations, target_val) {
    n <- 100
    results <- replicate(num_simulations, { max(rpois(lambda=.5, n=n))} )
    print(sum(results>=target_val))
    return(sum(results>=target_val) / length(results))
}
simulate(10^7, 9)
# 10,000,000 simulations 
# to prove P < 10^-7, we'll want 10^7 simulations

# 1.6
range <- 0:12
binom <- dbinom(range, prob=1/12, size=10)
barplot(binom, names.arg = range)
ois <- dpois(range, lambda=.5)
barplot(pois, names.arg = range)
chi_sq <- dchisq(0:100, 3*9)
barplot(chi_sq, names.arg = range)

# 5 non discrete distributions: look at distributions that don't require x to be an integer (e.g normal distribution)

# 1.7
instances <- rpois(100, 3)
mean(instances)
var(instances)

# 1.8
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "BSgenome.Celegans.UCSC.ce2"))

library(Biostrings)
library(BSgenome.Celegans.UCSC.ce2)
genome <- BSgenome.Celegans.UCSC.ce2
chrM <- genome[["chrM"]]
freq <- alphabetFrequency(chrM)[1:4]
# Let's find p val. This looks like testing goodness of fit

expected <- rep(.25*length(chrM), 4)
chi_sq_val <- sum((freq - expected)^2/expected)
df <- 4 - 1
pchisq(chi_sq_val, df, lower.tail=FALSE)
# p-val is 0, so we are very confident we can reject hypothesis
