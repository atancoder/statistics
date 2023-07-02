# Goal: Estimate the avg weight of chicks
# Population: ChickWeight dataset

library("dplyr")

data(ChickWeight)
weights <- ChickWeight$weight
sprintf("Number of weights is %d", length(weights))
population_mean_weight <- mean(weights)
sprintf("Population mean weight is %f grams", population_mean_weight)

# choose n to be large enough for C.L.T to apply
n <- 50
sampled_weights <- sample(weights, n, replace = TRUE)
estimated_mean_weight <- mean(sampled_weights)
sprintf("Sampled mean weight is %f grams", estimated_mean_weight)

# How confident are we about this prediction?
# There's a very low chance of predicting the exact
# mean, but we can provide a range with a
# confidence interval

# We first calculate SE for an estimate of how far we
# are from the actual mean weight
# SE = std dev / sqrt(n)
std_dev <- sd(sampled_weights)
std_error <- std_dev / (sqrt(n))

sprintf("Standard Error: %f", std_error)

# With SE, we can create a confidence interval
# For 95% confidence interval, we look 2 SEs away
lower_range <- estimated_mean_weight - 2 * std_error
upper_range <- estimated_mean_weight + 2 * std_error

sprintf(
    paste(
        "95%% confidence that the actual",
        "mean weight is in this interval [%f, %f]"
    ),
    lower_range,
    upper_range
)

is_in_interval <- between(population_mean_weight, lower_range, upper_range)
sprintf("Population mean in interval: %s", toString(is_in_interval))
