# Let's create a sample distribution plot
# Source
data(ChickWeight)
weights <- ChickWeight$weight

n <- 50
num_experiments <- 500

experiment_results <- replicate(
    num_experiments, sample(weights, n, replace = TRUE)
)

means <- apply(experiment_results, 2, mean)
sprintf("Means are [%s]", toString(means))

# Generate a sampling distribution histogram of the means
hist(
    means,
    main = "Sampling Distribution of the Mean",
    col = "lightblue",
    freq = FALSE,
    ylab = "Probability Density"
)
