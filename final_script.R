#1 Load required libraries
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(VennDiagram)
library(dplyr)
library(tidyr)
library(vegan)
library(ggpubr)
library(FSA)
library(gridExtra)
library(pairwiseAdonis)

set.seed(0)

#2 Europe map visualization
europe <- ne_countries(scale = "medium", continent = "Europe", returnclass = "sf")
poland <- europe[europe$sovereignt == "Poland", ]
sampling_site <- data.frame(longitude = 21.0122, latitude = 52.2297)

ggplot(data = europe) +
  geom_sf(fill = "lightgray") +
  geom_sf(data = poland, fill = "lightblue", color = "black") +
  geom_point(data = sampling_site, aes(x = longitude, y = latitude), color = "#ffe599", size = 3) +
  labs(x = "", y = "") +
  coord_sf(xlim = c(-12, 50), ylim = c(35, 72), expand = FALSE) +
  theme_minimal()

#3 Microscopy
mikro <- read.csv('mikroskop.csv')
kruskal.test(nb_morph ~ nest, data = mikro)
dunnTest(nb_morph ~ nest, data = mikro, method = "bh")

mikro_morph <- ggplot(mikro, aes(x = nest, y = nb_morph)) +
  geom_boxplot() +
  labs(title = '', y = 'Number of Morphotypes', x = 'Nest ID')

#4 Culture based method
szalki <- read.csv('szalki.csv')
szalki$nest <- factor(szalki$nest, labels = c('KOR', 'FAL', 'CHO', 'ZAB', 'ZIE'), levels = c('CHO', 'FAL', 'KOR', 'ZAB', 'ZIE'))

#4.a. CFU per sample
kruskal.test(CFU ~ nest, data = szalki)

plot1 <- ggplot(szalki, aes(x = nest, y = CFU)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(title = '', y = 'Fungal CFUs per IBP', x = expression(paste(italic("Formica polyctena"), " colony"))) +
  theme_minimal(base_size = 16) +
  stat_compare_means()
ggsave("CFU_szalki_21.01.2025.png", plot1)

#4.b. Number of morphotypes
kruskal.test(nb_morph ~ nest, data = szalki)

szalki_morph <- ggplot(szalki, aes(x = Nest, y = nb_morph)) +
  geom_boxplot() +
  labs(title = '', y = 'Number of Morphotypes', x = 'Nest ID')

#5 Metabarcoding

#5.a. Venn diagram
venn <- read.csv('genus_poprawione_v3.csv')
set1 <- as.vector(venn$X18S)
set2 <- as.vector(venn$ITS2)
set3 <- as.vector(venn$ITS1)

venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("18S \n n = 121", "ITS2 \n  n = 464", "ITS1 n = 421"),
  main = 'Total',
  sub = 'n = 577',
  imagetype = 'svg',
  filename = 'venn.svg',
  output = TRUE
)

#5.b. Barplots
# Fungal phyla barplot
data <- read.csv('fungal phyla.csv')
long_data <- data %>%
  pivot_longer(cols = ITS2:X18S)

bar_colors_Type <- c('#f1a983', '#47d45a', '#992189', '#60cbf3', '#d8d8d8', '#7f7f7f')
long_data$name <- factor(long_data$name, levels = c('ITS2', 'ITS1', 'X18S'), labels = c('ITS2', 'ITS1', '18S'))

plot1 <- ggplot(long_data, aes(x = name, y = value, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  scale_fill_manual(values = bar_colors_Type) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "", x = "Barcode", y = "Average Relative Abundance", fill = "Phylum")

# Bacteria class barplot
data <- read.csv('bacteria_class.csv')

data$class <- factor(data$class, levels = c('unassigned', 'other Classes', 'Actinobacteria', 'Gammaproteobacteria', 'Alphaproteobacteria', 'Bacilli'))

plot3 <- ggplot(data, aes(x = rel_abu, y = '', fill = class)) +
  geom_bar(stat = "identity", position = "fill", width = 1) +
  scale_fill_manual(values = data$color) + # Custom colors
  labs(x = "Average Relative Abundance", y = "", fill = 'Class') +
  scale_x_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(), # Remove y-axis labels
    axis.ticks.y = element_blank(), # Remove y-axis ticks
    panel.grid = element_blank(), # Remove grid lines
    legend.position = "bottom" # Position legend on the right
  )

#5.c. OTU count
tax.tab <- read.csv('ITS1.csv') #ITS1.csv/ITS2.csv/18s.csv/16S.csv
otu.tab <- tax.tab[,-c(2,3,4)]
row.names(otu.tab) <- otu.tab[,1]
otu.tab <- otu.tab[,-1]
otu.tab.t <- t(otu.tab)

otu <- specnumber(otu.tab.t)
meta <- read.csv('meta_ITS1.csv') #meta_ITS.csv ...

data <- cbind(data, meta)

kruskal.test(otu ~ nest, data = data)

dunnTest(otu ~ nest, data = data, method = "holm")

ggplot(data, aes(x = nest, y = otu, fill = nest)) +
  geom_boxplot() +
  labs(title = '', y = 'Number of zOTUs', x = '', fill = 'Colony') +
  theme_minimal() +
  stat_compare_means()

#5.d. PERMANOVA
tax.tab <- read.csv('ITS1.csv') #ITS1.csv/ITS2.csv/18s.csv/16S.csv
otu.tab <- tax.tab[,-c(2,3,4,5)]
row.names(otu.tab) <- otu.tab[,1]
otu.tab <- otu.tab[,-1]
otu.tab.t <- t(otu.tab)

meta <- read.csv('meta_ITS.csv') #ITS ...

dist <- vegdist(otu.tab.t, method = "bray")

bd <- betadisper(dist, meta$nest)
anova(bd)

perma <- adonis2(dist ~ nest, data = meta)

# Pariwaise post hoc
pair.substrat <- pairwise.adonis(dist, factors = meta$nest, p.adjust.m = 'bonferroni')

#5.e. NMDS
tax.tab <- read.csv('ITS2.csv')
otu.tab <- tax.tab[,-c(2,3,4)]
row.names(otu.tab) <- otu.tab[,1]
otu.tab <- otu.tab[,-1]
otu.tab.t <- t(otu.tab)

grudki_NMDS <- metaMDS(otu.tab.t, distance = 'bray', autotransform = TRUE)

meta <- read.csv('meta_ITS.csv')
nmds <- cbind(scores(grudki_NMDS, display = "sites"), meta)

ggplot(nmds, aes(x = NMDS1, y = NMDS2, color = nest, shape = nest)) +
  geom_point(size = 1, alpha = 0.8) +
  stat_ellipse(linetype = 2, linewidth = 1) +
  labs(title = "", color = 'Colony', shape = 'Colony') +
  theme_minimal()

#6 ecological specificity of bacterial and fungal communities

eco=read.csv('eco.csv')

ecolong=eco%>%pivot_longer(cols=Ant.associated:No.identical.match)%>%group_by(Type, name)%>%
  summarize(count = sum(value>0)/n())

ecolong$Type = factor(ecolong$Type, labels = c('Bacteria', 'Fungi'))
ecolong$name = factor(ecolong$name, levels = c("Ant.associated", "Other.insect.associated", "Other.animals", "Plants.associated" ,"Wood.decay","Fungal.associated","Environmental","Food.industry","No.identical.match"),
                      labels = c("Ant", "Other insect", "Other animals", "Plants" ,"Dead wood","Fungi","Environmental","Food industry","No identical match"))

bar_colors=c('#27447a', '#dec264')

ggplot(ecolong, aes(x=name, y=count,  fill = Type))+
  geom_bar(stat = 'identity', position=position_dodge(), width = 0.7)+
  scale_y_continuous(limits = c(0,1), labels = scales::percent_format())+
  theme_minimal()+
  labs(y='Percentage of core ASVs', x = 'Isolation source', fill = 'Kingdom')+
  theme(axis.text.x=element_text(angle=35, hjust=1),
        legend.position = c(0.9,0.8))+
  scale_fill_manual(values= bar_colors)
