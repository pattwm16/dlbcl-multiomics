# Adjacency Matrix Analysis Script

# Load in the .RData file
load("path_to_SmCCNetWeights.RData")

# Create list of threshold
remaining_edges <- list()
for (thresh in seq(0, 1, by=0.1)){
  remaining_edges <- c(remaining_edges, length(Abar@x[Abar@x > thresh]))
}

df <- as.data.frame(cbind("Threshold" = seq(0,1,by=0.1), 
                          "Remaining Edges" = remaining_edges))

library(ggplot2)
ggplot(df, aes(y=Remaining.Edges, x=Threshold)) + 
  geom_line() +
  geom_point(col="blue") +
  geom_label(aes(label=Remaining.Edges), nudge_y = 0.3) +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  ylab("Remaining Edges") +
  ggtitle("Edges remaining after trimming with threshold") +
  theme_classic()


