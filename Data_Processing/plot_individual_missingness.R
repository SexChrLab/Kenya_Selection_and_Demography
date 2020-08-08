library(ggplot2)

# high and low
data = data.frame(percentage = c(55.04087, 55.19126, 56.01093, 56.55738, 55.7377, 56.01093, 56.01093, 54.91803, 55.19126, 54.76839, 54.76839, 55.31335, 56.40327, 55.85831, 54.76839, 52.86104, 53.67847, 55.58583, 54.76839, 52.73224, 56.28415, 51.91257), chromosomes = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"))

data$chromosomes = as.character(data$chromosomes)
data$chromosomes = factor(data$chromosomes, levels = unique(data$chromosomes))

png("c://Users/tuyen/Documents/postdoc_asu/projects/Kenya_Selection_and_Demography/Data_Processing/plots/individual_missingness_percent.png", width = 10, height = 5.5, units = "in", res = 300)
ggplot(data, aes(x=chromosomes, y=percentage)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x="Chromosome", y="Percentage", title="Percentage of individuals with more than 10% missing data") +
  coord_cartesian(ylim=c(0, 100)) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# high
