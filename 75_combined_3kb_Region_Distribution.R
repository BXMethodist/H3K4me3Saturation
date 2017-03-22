library(ggplot2)
library(extrafont)

setwd('/Users/boxia/PycharmProjects/H3K4me3Saturation')
df = read.csv('75_combined_3kbstats.tsv', sep='\t')

df$width = log10(df$width)

ggplot(df, aes(x=factor(number.of.clusters), y=width))+ stat_boxplot(geom = "errorbar", width = 0.5) + geom_boxplot() + 
  labs(x='Number of Clusters', y = 'log10 of Width of Region (bp)') +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=10, face='bold'),
        axis.title=element_text(size=10, face='bold')) 
  
ggsave("Width_vs_NumCluster_75_3kb.eps", dpi = 600, width=3, height=3)


ggplot(df, aes(x=factor(number.of.clusters), y=number.of.samples))+ stat_boxplot(geom = "errorbar", width = 0.5) + geom_boxplot() + 
  labs(x='Number of Clusters', y = 'Number of Samples') +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=10, face='bold'),
        axis.title=element_text(size=10, face='bold')) 

ggsave("NumSample_vs_NumCluster_75_3kb.eps", dpi = 600, width=3, height=3)

