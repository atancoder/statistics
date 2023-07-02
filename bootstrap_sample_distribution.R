# Let's create a sample distribution plot via bootstrapping
# Source
data(ChickWeight)
weights <- ChickWeight$weight

n <- 50
num_experiments <- 1000

single_experiment <- sample(weights, n, replace = TRUE)

# simulate via bootstrap
experiment_results <- replicate(
    num_experiments, sample(single_experiment, n, replace = TRUE)
)

means <- apply(experiment_results, 2, mean)

# Generate a sampling distribution histogram of the means
hist(
    means,
    main = "Sampling Distribution of the Mean",
    col = "lightblue",
    freq = FALSE,
    ylab = "Probability Density"
)

# When you run the experiment a lot of times
# you can create a normal sample distribution
# std_dev of the means will be similar to SE
std_dev_of_means <- sd(means)

one_experiment <- experiment_results[, 1]
std_dev_of_samples <- sd(one_experiment)
std_error <- std_dev_of_samples / sqrt(length(one_experiment))

sprintf(
    paste(
        "Std_dev of Means: %f", "Std_error: %f",
        sep = ". "
    ),
    std_dev_of_means,
    std_error
)
