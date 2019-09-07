#Install packages
#Find more about mmgenome package in https://kasperskytte.github.io/mmgenome2/
install.packages("Rtsne") # Install Rtsne package from CRAN
if(!require(devtools)) install.packages("devtools") # If not already installed
devtools::install_github("jkrijthe/Rtsne") # more info in https://github.com/jkrijthe/Rtsne/tree/openmp

#check for remotes
if(!require(remotes))
  install.packages("remotes")

#install mmgenome2 using remotes
remotes::install_github("kasperskytte/mmgenome2")

#In case you havn't installed all the needed packages, they can be installed via e.g. `install_github("MadsAlbertsen/mmgenome/mmgenome")
#invoke libraries 
library (ggplot2)
library(mmgenome2)
library(Rtsne)
library(readxl)
options(scipen = 8)
library(knitr)
library(tidyverse)
library(magrittr)
library(gridExtra)
library(cowplot)
library(Biostrings)
library(kableExtra)


#load metagenomics data
#To load all the files neccesary for metagenomics analysis open the file "mmdata.Rdata" in R studio
assembly <- readDNAStringSet("data/assembly.fasta", format = "fasta")
ess    <- read.table("data/essential.txt", header = T)[,c(1,3)] #load essential genes
tax    <- read.delim("data/tax.txt", header = T) #taxonomy
#Load paired end connections
pe <- read.csv("data/network.txt",sep = "\t")# %>% 
#group_by(scaffold1,scaffold2)
#Load additional data (i.e., 16S)
#additional <- readDNAStringSet("data/16S.fa", format = "fasta")

#metadata
metadata<-readxl::read_excel(path = paste0("metadata.xlsx")) %>%
  mutate(Names = paste(Details, Reactor,sep = "_"))



#Merge data into a  single dataframe + Include coverage data
d <- mmload(
  assembly        = assembly, 
  list(Baseline1V_3= read.csv("data/18100FL-01-13_cov.csv", header = T)[,c("Name","Average.coverage")],
    Baseline1V_4= read.csv("data/18100FL-01-14_cov.csv", header = T)[,c("Name","Average.coverage")],
      Baseline1V_6= read.csv("data/18100FL-01-15_cov.csv", header = T)[,c("Name","Average.coverage")],
      Inoculum_0= read.csv("data/18100FL-01-16_cov.csv", header = T)[,c("Name","Average.coverage")],
      mth0.7V_3= read.csv("data/MQ180819-30_cov.csv", header = T)[,c("Name","Average.coverage")],
      mth0.7V_4=read.csv("data/MQ180819-31_cov.csv", header = T)[,c("Name","Average.coverage")],
      mth0.7V_6=read.csv("data/MQ180819-32_cov.csv", header = T)[,c("Name","Average.coverage")]
    ),
  taxonomy        = tax,
  essential_genes = ess,
  #additional = additional, 
  kmer_BH_tSNE    = T,
  kmer_pca        = T
  )


#general overview of the data frame
d

#Basic statistics of the metagenomics data
#Basic metagenome stats
#Basic stats of the meta genome libraries and assembly. 
#**Number of scaffolds (#)** is the number of scaffolds in the combined meta genome assembly 
#**Mean GC content (%)** is the mean GC content of all scaffolds in the assembly weighted by 
#scaffold length. **N50 (bp)** is a median statistic that indicates that 50% of the entire assembly 
#is contained in scaffolds equal to or larger than this value (bp).**Length Total (Mbp)** is the total
#combined length of the meta genome assembly in Mbp. **Length Maximum (bp)** is the length of the largest
#scaffold in the assembly. **Length mean** is the mean length of scaffolds in the assembly. 
mmstats(d)


#Mehtanobacteria bin extraction (bin 1)
#Differential coverage plot of the assembled meta genomic scaffolds (>5000 bp), 
#highlighting the initial extraction of scaffolds for **Bin 1**. 
#The size of the circles represent the length of the scaffolds. 
#Colours of the circles represent phylum level taxonomic classification, 
#based on the essential genes identified on the scaffolds, scaffolds with no colour 
#either contains no essential genes or could not be assigned a phylum level classification. 
#The x and y-axes show the sequencing coverage in the samples (log-scaled).

slct <- data.frame(
  cov_Baseline1V_3 = c(149.567, 149.567, 436.568, 953.968, 734.486, 222.012),
  cov_Baseline1V_4 = c(247.298, 130.046, 126.871, 425.993, 733.81, 577.864))

mmplot(mm = d,
       #locator = T,
       selection = slct,
       x = "cov_Baseline1V_3",
       y = "cov_Baseline1V_4",
       #min_length = 10000,
       color_by   = "class",
       x_scale    = "log10",
       y_scale    = "log10") +
  xlab("Coverage (Baseline -1V; Reactor 3)") +
  ylab("Coverage (Baseline -1V; Reactor 4)") +
  scale_x_log10(limits = c(1,1000),breaks = c(1, 10, 100,1000)) + 
  scale_y_log10(limits = c(1,1000),breaks = c(1, 10, 100,1000)) + 
  scale_size_area(max_size = 20, breaks = c(10000, 50000, 250000), name = "Scaffold Length (kbp)", label =  c(10, 50, 250))+
  scale_color_discrete(name = "Taxonomic Classification") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5, shape = 19), nrow = 3, title.position="top")) +
  guides(size = guide_legend(nrow = 3, title.position="top")) +
  theme(legend.position = "bottom", 
        legend.justification = "center",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10, face = "bold"))

#Save the differential coverage plot
ggsave(filename="Bin1_-1V_R3_R4.pdf", width = 8, height = 8, dpi=300) 

#Extract the bin
d_bin1 <- mmextract(d,slct)


#Expanding the bin by paired-end connections
#The initial bin extraction may not capture all scaffolds associated to the bin, 
#due to e.g. high-coverage repeats or a slightly different coverage profile than the rest of the genome. 
#To adjust for this, we use paired-end connections between scaffolds to draw in additional scaffolds.
#Figure shows all scaffolds in the expanded network, the size of the circles
#represent the length of the scaffolds and the width of the lines between scaffolds indicates the number 
#of paired-end connections between them. Points are coloured by their coverage, with several 
#large scaffolds of low coverage indicating microdiversity. In Figure the 
#expanded network are shown in context of the metagenome highlighted in dark red.

dNW <- mmexpand_network(mm = d, scaffolds = d_bin1, network = pe, min_connections = 20)

p1 <- mmnetwork(mm = dNW,
                network = pe,
                min_connections = 20,
                color_by = "cov_Baseline1V_3",
                seed = 80085) +
  scale_size_area(max_size = 20, 
                  breaks = c(5000,10000, 25000, 50000), 
                  name = "Scaffold \nLength (kbp)", 
                  label =  c(5, 10, 25, 50)) +
  scale_color_viridis_c(name = "Coverage \n(Fold)") +
  theme(legend.position = "right", 
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10, face = "bold"))


p2 <- mmplot(mm = dNW,             
             x = "cov_Baseline1V_3",
             y = "cov_Baseline1V_4",
             color_by = "class") +
  xlab("Coverage (Baseline -1V; Reactor 3)") +
  ylab("Coverage (Baseline -1V; Reactor 4)") +  
  xlim(0,700) +
  ylim(0,600) +
  scale_size_area(max_size = 20, 
                  breaks = c(5000,10000, 25000, 50000), 
                  name = "Scaffold \nLength (kbp)", 
                  label =  c(5, 10, 25, 50)) +
  scale_color_discrete(name = "Taxonomic\nClassification") +
  theme(legend.position = "right", 
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8))

plot_grid(p1,p2,nrow = 2)


# Export fasta file of the Bin
mmexport(dNW, assembly = assembly,file = "Bin1.fasta")



#Statistics of the extarcted bin:
#Basic stats of the extracted genome bin. **Number of scaffolds (#)** is the number of scaffolds
#in the extracted genome bin. **N50 (bp)** is a median statistic that indicates that 50% of the entire 
#assembly is contained in scaffolds equal to or larger than this value (bp).**Length Total (bp)** is the 
#total combined length of extracted genome bin. **Length max (bp)** is the length of the largest scaffold 
#in the extracted genome bin. **Length mean** is the mean length of scaffolds in the extracted genome bin.
#**Mean GC content (%)** is the mean GC content of all scaffolds in the bin weighted by scaffold length. 
#The **cov_** prefix shows the average coverage of the bin in all available samples. **Total essential genes
#(#)** is the total number of identified essential "single copy" genes in the genome bin, **Unique essential
#genes (#)** is the number of unique essential "single copy" genes. The list of essential "single copy" genes
#we look for contain 108 genes in total. If the total number of genes is much higher than the number of
#unique genes, it is likely that the bin contains scaffolds from more than one organism.

mmstats(dNW)




#Differential coverage plot Inoculum vs R3 (-1V)

slct <-  data.frame(cov_Baseline1V_3 = c(387.572, 992.919, 1635.732, 1111.437, 310.327, 119.198, 125.706),
                    cov_Inoculum_0 = c(2.716, 2.641, 1.641, 0.709, 0.362, 0.616, 1.552))

mmplot(mm = d,
       #locator = T,
       selection = slct,
       x = "cov_Baseline1V_3",
       y = "cov_Inoculum_0",
       #min_length = 10000,
       color_by   = "class",
       x_scale    = "log10",
       y_scale    = "log10") +
  xlab("Coverage (Baseline -1V; Reactor 3)") +
  ylab("Coverage (cov_Inoculum_0)") +
  scale_x_log10(limits = c(0.001,1700),breaks = c(0.01, 1, 100,1000)) + 
  scale_y_log10(limits = c(0.001,150),breaks = c(0.01, 1, 100)) + 
  scale_size_area(max_size = 20, breaks = c(10000, 50000, 250000), name = "Scaffold Length (kbp)", label =  c(10, 50, 250))+
  scale_color_discrete(name = "Taxonomic Classification") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5, shape = 19), nrow = 3, title.position="top")) +
  guides(size = guide_legend(nrow = 3, title.position="top")) +
  theme(legend.position = "bottom", 
        legend.justification = "center",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10, face = "bold"))

#Save the differential coverage plot
ggsave(filename="Inoculum_R3_1V.pdf", width = 8, height = 8, dpi=300) 




#Differential coverage plot  R3 (-1V) vs R6 (-1V)


slct <- data.frame(cov_Baseline1V_3 = c(517.589, 780.036, 873.309, 766.949, 316.09, 207.471, 189.873, 289.278),
                   cov_Baseline1V_6 = c(492.44, 473.45, 325.452, 192.569, 84.422, 111.131, 197.442, 333.688))
mmplot(mm = d,
       #locator = T,
       selection = slct,
       x = "cov_Baseline1V_3",
       y = "cov_Baseline1V_6",
       #min_length = 10000,
       color_by   = "class",
       x_scale    = "log10",
       y_scale    = "log10") +
  xlab("Coverage (Baseline -1V; Reactor 3)") +
  ylab("Coverage (Baseline -1V; Reactor 6)") +
  scale_x_log10(limits = c(1,1300),breaks = c(1, 10, 100,1000)) + 
  scale_y_log10(limits = c(1,1300),breaks = c(1, 10, 100,1000)) + 
  scale_size_area(max_size = 20, breaks = c(10000, 50000, 250000), name = "Scaffold Length (kbp)", label =  c(10, 50, 250))+
  scale_color_discrete(name = "Taxonomic Classification") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5, shape = 19), nrow = 3, title.position="top")) +
  guides(size = guide_legend(nrow = 3, title.position="top")) +
  theme(legend.position = "bottom", 
        legend.justification = "center",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10, face = "bold"))

#Save the differential coverage plot
ggsave(filename="R3_R6_1V.pdf", width = 8, height = 8, dpi=300) 



#Differential coverage plot Inoculum vs R3 (-0.7) vs R3 (-1V)


slct <- data.frame(cov_Baseline1V_3 = c(484.863, 844.619, 905.283, 692.573, 286.459, 184.23, 194.05, 293.995),
                   cov_mth0.7V_3 = c(1422.209, 1351.288, 661, 345.053, 180.124, 211.909, 462.304, 945.093))
mmplot(mm = d,
       #locator = T,
       selection = slct,
       x = "cov_Baseline1V_3",
       y = "cov_mth0.7V_3",
       #min_length = 10000,
       color_by   = "class",
       x_scale    = "log10",
       y_scale    = "log10") +
  xlab("Coverage (Baseline -1V; Reactor 3)") +
  ylab("Coverage (cov_mth0.7V_3)") +
  scale_x_log10(limits = c(1,1500),breaks = c(1, 10, 100,1000)) + 
  scale_y_log10(limits = c(1,1500),breaks = c(1, 10, 100,1000)) + 
  scale_size_area(max_size = 20, breaks = c(10000, 50000, 250000), name = "Scaffold Length (kbp)", label =  c(10, 50, 250))+
  scale_color_discrete(name = "Taxonomic Classification") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5, shape = 19), nrow = 3, title.position="top")) +
  guides(size = guide_legend(nrow = 3, title.position="top")) +
  theme(legend.position = "bottom", 
        legend.justification = "center",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10, face = "bold"))

#Save the differential coverage plot
ggsave(filename="R3_0month_1month.pdf", width = 8, height = 8, dpi=300) 



#Bin 1 abundance. 1V vs 1 month -0.7V

#Import data
bin1_abundance <- read_csv("~/Desktop/Dario Rangel Shaw/ANAMMOX research/Cathode set potential project/Metagenomics analysis cathode set potential/bin1_abundance.csv")

ggplot(bin1_abundance, aes(x=Coverage, y=Percentage, color= Coverage)) + 
  geom_boxplot()+
  scale_x_discrete(limits=c("Baseline1V", "1mth0.6V"))+ 
  theme_classic()+
  stat_boxplot() 

#ANOVA
mod1<-aov(Percentage~Coverage,data=bin1_abundance)
summary(mod1)
#Not siginificant difference between Bin 1 abundance at 1V and 1 month at -0.7V





