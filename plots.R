# Box Plot?

data(ChickWeight)
weights <- ChickWeight$weight

pdf("chick_weight_boxplot.pdf")
boxplot(weights)
dev.off()
