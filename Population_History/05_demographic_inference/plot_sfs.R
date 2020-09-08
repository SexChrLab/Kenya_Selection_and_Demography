library(ggplot2)
library(ggpubr)

empirical = read.csv("c://Users/tuyen/Documents/postdoc_asu/projects/Kenya_Selection_and_Demography/Population_History/04_generate_sfs/sfs_nre/autosomes_sfs_nre.txt", header = FALSE, sep="\t")
model_1epoch = read.csv("c://Users/tuyen/Documents/postdoc_asu/projects/Kenya_Selection_and_Demography/Population_History/05_demographic_inference/1D_1Epoch_NRE/Turkana.dadi.inference.1D.1Epoch.runNum.1.expSFS.normalized.by.theta", header=FALSE, sep = "\t")
model_2epoch = read.csv("c://Users/tuyen/Documents/postdoc_asu/projects/Kenya_Selection_and_Demography/Population_History/05_demographic_inference/1D_2Epoch_NRE/Turkana.dadi.inference.1D.2Epoch.runNum.3.expSFS.normalized.by.theta", header = FALSE, sep="\t")
model_3epoch = read.csv("c://Users/tuyen/Documents/postdoc_asu/projects/Kenya_Selection_and_Demography/Population_History/05_demographic_inference/1D_1Bottleneck_NRE/Turkana.dadi.inference.1D.1Bottleneck.runNum.20.expSFS.normalized.by.theta", header = FALSE, sep="\t")

empirical_new = data.frame(bins=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10"), vals = c(empirical[,2][1:10], sum(empirical[,2][11:108])))
model_1epoch_new = data.frame(bins=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10"), vals = c(model_1epoch[,2][1:10], sum(model_1epoch[,2][11:108])))
model_2epoch_new = data.frame(bins=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10"), vals = c(model_2epoch[,2][1:10], sum(model_2epoch[,2][11:108])))
model_3epoch_new = data.frame(bins=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10"), vals = c(model_3epoch[,2][1:10], sum(model_3epoch[,2][11:108])))

empirical_labels = rep("Empirical", 11)
oneepoch_labels = rep("1Epoch", 11)
twoepoch_labels = rep("2Epoch", 11)
threeepoch_labels = rep("3Epoch", 11)

dat = data.frame(Bins = empirical_new[,1], Counts=c(empirical_new[,2], model_1epoch_new[,2], model_2epoch_new[,2], model_3epoch_new[,2]), category = c(empirical_labels, oneepoch_labels, twoepoch_labels, threeepoch_labels))

dat$category = as.character(dat$category)
dat$category = factor(dat$category, levels = unique(dat$category))

dat$Bins = as.character(dat$Bins)
dat$Bins = factor(dat$Bins, levels = unique(dat$Bins))

ggplot(data=dat, aes(x=Bins, y=Counts, fill=category)) + 
  theme_bw() +
  geom_bar(stat = "identity", position=position_dodge(), colour="black") + 
  scale_fill_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73")) + 
  theme(axis.text.x = element_text(size=14, colour="black"), axis.text.y = element_text(size=14, colour = "black"), axis.title.y = element_text(size=18), axis.line=element_line(colour = "black"),plot.title = element_text(hjust = 0.5, size=18), legend.title = element_blank(), panel.background = element_blank(), legend.text = element_text(size=14), legend.spacing.x = unit(0.15, 'cm')) +
  labs(y="Counts", x="")
ggsave("c://Users/tuyen/Documents/postdoc_asu/projects/Kenya_Selection_and_Demography/Population_History/05_demographic_inference/plots/compare_empirical_model_20200908.png", width = 9, height = 5, units = "in")

