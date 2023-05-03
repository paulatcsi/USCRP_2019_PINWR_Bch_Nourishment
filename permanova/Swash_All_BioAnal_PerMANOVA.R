###################################################################################
# Script: PermANOVA Analysis of beach species abundance data                      #
# 		  created by Paul Paris                                    2020/02				  #
# 						                                        					  		  	  		  #
###################################################################################

### ------------ Load data and package -----------------------------------------
### Load vegan package
library(vegan)
library(dplyr)
library(ggplot2)
library(vegan3d)
library(xtable)

### Load the data
proj_path = '/Users/parisp15/OneDrive - East Carolina University/Projects/USCRP_Pea_Island_2020_2021/'
data_path = 'miscellaneous/Final Results Paper/Final_Data_Sets/'   #use_these_data_for_analysis/'
data_file_name = 'NCDOT_CSI_Bio_Counts_R.csv'      #'bio_data_2_R.csv'


### ------------ perMANOVAs (All Swash Species, No O quadrata) ----------------------
###  for inter-taxa comparisons:
df <- read.csv(file = paste(proj_path, data_path, data_file_name, sep=""))
df <- filter(df, E_talpoida != 0 | D_variabilis != 0 | Amphipods != 0 | S_squamata != 0)

taxa <- sqrt( df[c(10,11,12,13,14)])
taxa.env <- df[c(2:4,17)]

### distance matrix
taxa.dist <- vegdist(taxa, method='bray')   # all swash species/groups


# all swash species/groups
pmp <- adonis2(taxa.dist ~ epoch*area*season*year, data=taxa.env, permutations=999, method="bray")
pmp
pmp_tab <- xtable(pmp, caption='Permutational ANOVA Results - All Inverts')
print(pmp_tab,file="All_Species_PANOVA.tex",append=T,table.placement = "h", caption.placement="bottom", hline.after=seq(from=-1,to=nrow(pmp_tab),by=1))


### ------------ perMANOVAs (Individual Species, incld O quadrata) ----------------------
df <- read.csv(file = paste(proj_path, data_path, data_file_name, sep=""))

#df <- filter(df, year <= 2015)
taxa.et <- sqrt( df[c(11)])   # E talpoida
taxa.dv <- sqrt( df[c(12)])   # D variabilis
taxa.an <- sqrt( df[c(14)])  # S squamata
taxa.am <- sqrt( df[c(13)])  # Amphipods
taxa.oq <- sqrt( df[c(10)])   # O quadrata
taxa.env <- df[c(2:4,17)]    # environmernt vars.

### distance matrices
taxa.et.dist <- vegdist(taxa.et, method='euclidean')        # E talpoida only
taxa.dv.dist <- vegdist(taxa.dv, method='euclidean')        # Donax only
taxa.an.dist <- vegdist(taxa.an, method='euclidean')        # S squamata only
taxa.am.dist <- vegdist(taxa.am, method='euclidean')        # amphipods only
taxa.oq.dist <- vegdist(taxa.oq, method='euclidean')        # O quadrata only


# for E talpoida
et <- adonis2(taxa.et.dist ~ epoch*area*season*year, data=taxa.env, permutations=999, method="euclidean")
et
pmp_tab <- xtable(et, caption='Permutational ANOVA Results - E talpoida')
print(pmp_tab,file="E_talpoida_PANOVA.tex",append=T,table.placement = "h", caption.placement="bottom", hline.after=seq(from=-1,to=nrow(pmp_tab),by=1))

# for Donax variabilis
dv <- adonis2(taxa.dv.dist ~ epoch*area*season*year, data=taxa.env, permutations=999, method="euclidean")
dv
pmp_tab <- xtable(dv, caption='Permutational ANOVA Results - D variabilis')
print(pmp_tab,file="D variabilis_PANOVA.tex",append=T,table.placement = "h", caption.placement="bottom", hline.after=seq(from=-1,to=nrow(pmp_tab),by=1))

# for Solelepsis squamata (worms)
an <- adonis2(taxa.an.dist ~ epoch*area*season*year, data=taxa.env, permutations=999, method="euclidean")
an
pmp_tab <- xtable(an, caption='Permutational ANOVA Results - S squamata')
print(pmp_tab,file="S squamata_PANOVA.tex",append=T,table.placement = "h", caption.placement="bottom", hline.after=seq(from=-1,to=nrow(pmp_tab),by=1))

# for Amphipods
am <- adonis2(taxa.am.dist ~ epoch*area*season*year, data=taxa.env, permutations=999, method="euclidean")
am
pmp_tab <- xtable(am, caption='Permutational ANOVA Results - Amphipods undifferentiated')
print(pmp_tab,file="Amphipodss_PANOVA.tex",append=T,table.placement = "h", caption.placement="bottom", hline.after=seq(from=-1,to=nrow(pmp_tab),by=1))

#for O quadrata
oq <- adonis2(taxa.oq.dist ~ epoch*area*season*year, data=taxa.env, permutations=999, method="euclidean")
oq
pmp_tab <- xtable(oq, caption='Permutational ANOVA Results - O quadrata')
print(pmp_tab,file="O quadrata_PANOVA.tex",append=T,table.placement = "h", caption.placement="bottom", hline.after=seq(from=-1,to=nrow(pmp_tab),by=1))


### ------------ NMDS ----------------------------------------------------------
### Run NMDS (default is bray)
nmds <- metaMDS(taxa, k=2, distance='bray', autotransform = FALSE)
nmds$stress



### --------PLOT NMDS PRE vs POST-----------------------------------------------
data_scores <- as.data.frame(scores(nmds))

scores <- cbind(as.data.frame(data_scores), epoch = taxa.env$epoch)
centroids <- aggregate(cbind(NMDS1,NMDS2) ~ epoch, data=scores, FUN=mean)
seg <- merge(scores, setNames(centroids, c('epoch','NMDS1','NMDS2')),by='epoch', sort=FALSE)

ggplot(scores, aes(x = NMDS1, y = NMDS2, colour = epoch )) +
  geom_segment(aes(x=NMDS1.x, y=NMDS2.x, xend=NMDS1.y, yend=NMDS2.y), data=seg) +
  
  geom_point(data = centroids, size = 3) +                    # add centroids
  geom_point(size=1) +                                              
  coord_fixed()+                                              
  theme_bw()+ 
  theme(legend.position="right",legend.text=element_text(size=10),legend.direction='vertical')




### ------------ CHECK HOMOGENEITY OF VARIANCE ---------------------------------
### ------------ DISPERSION TEST FOR TIME (EPOCH) ------------------------------

## This is the only assumption that must be satisfied for PermANOVA...
bd <- betadisper(taxa.dist, taxa.env$epoch)    # computes PCoA result!
boxplot(bd)
anova(bd)   # F-Test
permutest(bd) # permutation test   (not significant means pre and post dispersion is alike)


#EOF