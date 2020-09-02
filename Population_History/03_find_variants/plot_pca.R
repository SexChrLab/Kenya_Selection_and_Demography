# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SNPRelate")

# SNPrelate PCA
library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(ggpubr)

setwd("/scratch/tphung3/Kenya_Selection_and_Demography/Population_History/03_find_variants/filtered_variants/")

pop_code = c(rep('ASU', 83), rep('Princeton', 27))

vcf.fn = "autosomes.gatk.called.vqsr.select.variants.hwe.recode.pilot.masks.nre.neutral.vcf.gz"
snpgdsVCF2GDS(vcf.fn, "autosomes.gds", method="biallelic.only")
snpgdsSummary("autosomes.gds")
genofile <- snpgdsOpen("autosomes.gds")
pca = snpgdsPCA(genofile, autosome.only=TRUE)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

# Make a data frame
sample.id = pca$sample.id
autosomes_df <- data.frame(sample.id = pca$sample.id,
                      pop = factor(pop_code)[match(pca$sample.id, sample.id)],
                      EV1 = pca$eigenvect[,1],    # the first eigenvector
                      EV2 = pca$eigenvect[,2],    # the second eigenvector
                      EV3 = pca$eigenvect[,3],
                      EV4 = pca$eigenvect[,4],
                      EV5 = pca$eigenvect[,5],
                      stringsAsFactors = FALSE)


p1 = ggplot(autosomes_df, aes(x=EV1, y=EV2, fill=pop, shape=pop)) + 
  geom_point(size=3, colour="black") + 
  geom_text(aes(label=ifelse(abs(EV1)>0.6,as.character(sample.id),'')),hjust=0,vjust=0) +
  theme_bw()+
  labs(title="Autosomes", x="PC1", y="PC2") +
  theme(axis.text.x = element_text(size=14, colour="black"), axis.text.y = element_text(size=14, colour = "black"), axis.title.y = element_text(size=18, angle=0, vjust=0.5), axis.title.x = element_text(size=18), axis.line=element_line(colour = "black"),plot.title = element_text(hjust = 0.5, size=18), legend.title = element_blank(), panel.background = element_blank(), legend.text=element_text(size=14), legend.position = 'top') +
  scale_fill_manual(values=c("black", "red")) +
  scale_shape_manual(values=c(21, 21))
ggsave("/scratch/tphung3/Kenya_Selection_and_Demography/Population_History/03_find_variants/plots/pca.png", width = 9, height = 5, unit="in")

snpgdsClose(genofile)