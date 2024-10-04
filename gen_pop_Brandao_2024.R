# Some packages may need to be downloaded from GitHub, so you may need to install devtools to do it

if (!require('devtools'))     install.packages('devtools');        library('devtools')


# Here are all the packages that will be required for this script


if (!require('vcfR'))          install.packages('vcfR');           library('vcfR')
if (!require('poppr'))         install.packages('poppr');          library('poppr')
if (!require('ape'))           install.packages('ape');            library('ape')
if (!require('RColorBrewer'))  install.packages('RColorBrewer');   library('RColorBrewer')
if (!require('igraph'))        install.packages('igraph');         library('igraph')
if (!require('ggplot2'))       install.packages('ggplot2');        library('ggplot2')
if (!require('tidyverse'))     install.packages('tidyverse');      library('tidyverse')
if (!require('pegas'))         install.packages('pegas');          library('pegas')
if (!require('hierfstat'))     install.packages('hierfstat');      library('hierfstat')
if (!require('LEA'))           install.packages('LEA');            library('LEA')
if (!require('adegenet'))      install.packages('adegenet');       library('adegenet')
if (!require('grid'))          install.packages('grid');           library('grid')
if (!require('ggpubr'))        install.packages('ggpubr');         library('ggpubr')
if (!require('dartR'))         install.packages('dartR');          library('dartR')
if (!require('ggtern'))        install.packages('ggtern');         library('ggtern')
if (!require("SNPfiltR"))      install.packages("SNPfiltR");       library("SNPfiltR")
if (!require("coda"))          install.packages("coda");           library("coda")
if (!require("networkD3"))     install.packages("networkD3");      library("networkD3")
if (!require("directlabels"))  install.packages("directlabels");   library("directlabels")
if (!require("pheatmap"))      install.packages("pheatmap");       library("pheatmap")
if (!require("gplots"))        install.packages("gplots");         library("gplots")

#==============================================================================#
####'.'.' PART 1 - Filtering and preparing dataset for other analysis '.'.'#####
#==============================================================================#
rm(list=ls())

#----------------------- Setting the working directory ------------------------#


setwd("D:/gen_pop_anemonia/")


#------------------------------- Loading VCF ----------------------------------#

# This VCF has already been pre-filtered using VCFTOOLS on Hydra, 
# for info on filters applied check last section of this script


# Loading the vcf file

vcf_all <- read.vcfR("filtered_vcf_brazil_65.recode.vcf")

# Checking general information about the file

vcf_all



#----------------------- Adding pop and geo info for VCF ----------------------#


# This bit was adapted from the script of Dalapicolla 2023 available in https://github.com/jdalapicolla/PopGenPipe/tree/main/Filtering

# Load information from a csv file

geo_data = read.csv("./snp_filtering/infopop.csv", sep = ";") 
#head(geo_data)
#tail(geo_data)

# Select rows in the same order than vcf

genomic_order = as.data.frame(colnames(vcf_all@gt)[-1]) #always 1st line will be FORMAT. You need to remove #1? linha sempre ? esse FORMAT e n?o a primeira amostra
names(genomic_order) = "Sample_ID"
#genomic_order

geo_data_reorder = semi_join(geo_data, genomic_order, join_by(x$id == y$Sample_ID))
#geo_data_reorder
geo_data_reorder = inner_join(genomic_order, geo_data_reorder, join_by(x$Sample_ID == y$id))
#head(geo_data_reorder)

# Check line incongruencies:

identical(as.character(colnames(vcf_all@gt)[-1]), as.character(geo_data_reorder$Sample_ID))
setdiff(as.character(colnames(vcf_all@gt)[-1]), as.character(geo_data_reorder$Sample_ID))
setdiff(as.character(geo_data_reorder$Sample_ID), as.character(colnames(vcf_all@gt)[-1]))

popmap = geo_data_reorder[,1:2]#[,c(2,1)] #col "pop" as population/location.
#popmap
names(popmap) = c("id", "pop")
#head(popmap)


####------------------------ Removing clones from VCF ----------------------####

# Removing samples if necessary (in case there are clones) - Details on last section

colnames(vcf_all@gt)
remove_ind = c("Moufpe_880_8" , "Moufpe_883_3" , "Moufpe_883_7")
ind_to_keep = colnames(vcf_all@gt)[!colnames(vcf_all@gt) %in% remove_ind]

vcf_noclone = vcf_all[samples=ind_to_keep]
#vcf_noclone
vcf_noclone

####---------------------------- Filtering MAC -----------------------------####


MAF_5 = (length(colnames(vcf_noclone@gt)[-1])*5)/100 #5% das amostras

vcf_nc_maf = min_mac(vcf_noclone, min.mac = MAF_5)
vcf_nc_maf

####-------------------- Removing alleles out of balance --------------------####

vcf_nc_maf_ab = filter_allele_balance(vcf_nc_maf, min.ratio = 0.2, max.ratio = 0.8)

vcf_nc_maf_ab 

#-------------------------------- Converting VCF ------------------------------#


# Converting vcf to genlight

noclone_gl <- vcfR2genlight(vcf_nc_maf_ab) #used for DAPC, PCA

noclone_gl


#---------------------------------- Adding pop info ---------------------------#
       

pop.data <- read.table("snp_filtering/infopop.csv", sep = ";", header = TRUE)
pop.data <- pop.data[pop.data$id %in% noclone_gl@ind.names,]  

pop.data

# population
pop(noclone_gl) <- pop.data$pop

noclone_gl@pop


####------------------------- Removing loci out of HWE --------------------#####

# This bit was also adapted from the script of Dalapicolla 2023 available in https://github.com/jdalapicolla/PopGenPipe/tree/main/Filtering


gl_vcf_fil = gl.filter.hwe(
  noclone_gl,
  subset = 'each',
  n.pop.threshold = 2, # most of pops 
  method_sig = "Exact", # Wigginton et al. (2005)
  multi_comp = FALSE,
  multi_comp_method = "BH", # Benjamini & Hochberg, 1995
  alpha_val = 0.05,
  pvalue_type = "midp", # Graffelman & Moreno, 2013
  cc_val = 0.5,
  min_sample_size = 5, #minimum sample per population size 
  verbose = NULL
)

head(gl_vcf_fil@loc.names)
head(noclone_gl@loc.names)

noclone_hwe_gl = gl_vcf_fil

####------------------------- Removing monomorphs -------------------------#####


monomorph_out <- gl.filter.monomorphs(noclone_hwe_gl)

nc_hwe_nmono_gl = monomorph_out
nc_hwe_nmono_gl # No clone, Hardy-Weinberg equilibrium, no monomorphs genlight

nc_hwe_nmono_gl_nona <- gl.filter.allna(nc_hwe_nmono_gl)
nc_hwe_nmono_gl_nona

# Keep in vcf the same snps as in the filtered genlight.

filtered_snps = which(nc_hwe_nmono_gl_nona@loc.names %in% gl_vcf_fil@loc.names)
nc_hwe_nmono_vcf= vcf_noclone[filtered_snps, ]
nc_hwe_nmono_vcf

missing_by_snp(nc_hwe_nmono_vcf) #verify statistics  
nc_hwe_nmono_vcf1 = missing_by_snp(vcfR=nc_hwe_nmono_vcf, cutoff = 0.70)
nc_hwe_nmono_vcf1

popmap = popmap[popmap$id %in% colnames(nc_hwe_nmono_vcf1@gt),]


###-------------------------------- Filtering NAs ---------------------------###

missing_by_sample(vcfR=nc_hwe_nmono_vcf1, popmap = popmap)  #verify statistics #verificar as estat?sticas antes de escolher um limiar. No exemplo, 60%
vcf_noclone_clean1 = missing_by_sample(vcfR=nc_hwe_nmono_vcf1, cutoff = .2)
vcf_noclone_clean1

nc_hwe_nmono_gl<-vcfR2genlight(vcf_noclone_clean1)

pop(nc_hwe_nmono_gl) <- pop.data$pop
nc_hwe_nmono_gl@pop

####----------------------------- Saving files ----------------------------#####

# Converting to BayeScan format and saving on current working directory

#nc_hwe_nmono_gl <- gl.recalc.metrics(nc_hwe_nmono_gl) #I was having some issues with the gl.complience.check() funtion so I had to recalculate the metrics here for it to work

nc_hwe_nmono_gl <- gl.compliance.check(nc_hwe_nmono_gl)

# This file will be used on BayeScan

gl2bayescan(nc_hwe_nmono_gl, outpath = '.', outfile = '1267snps_36ind_asarg_bayescan.txt')

# Saving the filtered gl as VCF to be used when ploting BayScan graphs

vcfR::write.vcf(vcf_noclone_clean1, file = '1267snps_36ind_asarg.vcf') # saving this way causes issues

gl2vcf(nc_hwe_nmono_gl, outfile = '1267snps_36ind_asarg1', plink_path = 'D:/plink', outpath = './') #this one I saved to be used for creating the id fille used to exclude outliers

# Next step is to run BayeScan and obtain the outliers



####--------------- Removing the outliers (BayeScan results) ---------------####

# This function bellow comes along with the installation of BayeScan and it is used to 
# generate the graphs

# Run it first

plot_bayescan<-function(res,FDR=0.05,size=1,pos=0.35,highlight=NULL,
                        name_highlighted=F,add_text=T)
{
  if (is.character(res))
    res=read.table(res)
  
  colfstat=5
  colq=colfstat-2
  
  highlight_rows=which(is.element(as.numeric(row.names(res)),highlight))
  non_highlight_rows=setdiff(1:nrow(res),highlight_rows)
  
  outliers=as.integer(row.names(res[res[,colq]<=FDR,]))
  
  ok_outliers=TRUE
  if (sum(res[,colq]<=FDR)==0)
    ok_outliers=FALSE;
  
  res[res[,colq]<=0.0001,colq]=0.0001
  
  # plot
  plot(log10(res[,colq]),res[,colfstat],xlim=rev(range(log10(res[,colq]))),
       xlab="log10(q value)",ylab=names(res[colfstat]),type="n")
  points(log10(res[non_highlight_rows,colq]),res[non_highlight_rows,colfstat],
         pch=19,cex=size)
  
  if (name_highlighted) {
    if (length(highlight_rows)>0) {
      text(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],
           row.names(res[highlight_rows,]),col="red",cex=size*1.2,font=2)
    }
  }
  else {
    points(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],
           col="red",pch=19,cex=size)
    # add names of loci over p and vertical line
    if (ok_outliers & add_text) {
      text(log10(res[res[,colq]<=FDR,][,colq])+pos*(round(runif(nrow(res[res[,colq]<=FDR,]),1,2))*2-3),
           res[res[,colq]<=FDR,][,colfstat],row.names(res[res[,colq]<=FDR,]),cex=size)
    }
  }
  lines(c(log10(FDR),log10(FDR)),c(-1,1),lwd=2)
  
  return(list("outliers"=outliers,"nb_outliers"=length(outliers)))
}

#----------------------- Setting the working directory ------------------------#


#Should be set where BayeScan results are saved

setwd('D:/BayeScan2.1/')


#--------------------------- Check chain convergence --------------------------#

# Here I've run BayeScan two times

chain <- read.table("run1_1267_od100/bayescanLongRunOD100.sel",header=TRUE)
chain <- chain[-c(1)]
chain <- mcmc(chain,thin=10)


chain2 <- read.table("run2_1267_od100/bayescanLongRunOD100.sel",header=TRUE)
chain2 <- chain2[-c(1)]
chain2 <- mcmc(chain2,thin=10)


combined = mcmc.list(chain,chain2)
plot(combined)
gelman.diag(combined)
gelman.plot(combined)


# for explanation on convergence results, check end of document

plot_bayescan('run1_1267_od100/bayescanLongRunOD100_fst.txt', FDR = 0.05)

#Naming the row position of the outliers

outliers <- plot_bayescan('run1_1267_od100/bayescanLongRunOD100_fst.txt', FDR = 0.05)$outliers
outliers

# saving images

tiff(file="chain_convergence_bayescan.tiff",
     width=20, height=16, units="in", res=300)
gelman.plot(combined)


tiff(file="bayescan_qvalue_plot.tiff",
     width=20, height=16, units="in", res=300)
plot_bayescan('run1_1267_od100/bayescanLongRunOD100_fst.txt', FDR = 0.05)

dev.off()



#------------------ Preparing table to identify outliers ----------------------# 

setwd('D:/gen_pop_anemonia')

# Creating a table with the fst.txt file

bayescan <- read.table('D:/BayeScan2.1/run1_1267_od100/bayescanLongRunOD100_fst.txt')

# vcf used to create file below was the one saved from the nc_hwe_nmono_vcf

# bash command used was: grep -v "#" YOUR_VCF.vcf | cut -f 3 > OUTPUT_NAME.txt


SNPids=read.table('id-1267snps.txt',header=FALSE)

# Generating and saving the table

bayescan=cbind(SNPids, bayescan) 
colnames(bayescan)=c("SNP","PROB","LOG_PO","Q_VALUE","ALPHA","FST") 

#ploting the names of the outliers

outlier <- bayescan[outliers,]$SNP
outlier

#removing outlier(s)

null_set_vcf <- read.vcfR('1267snps_36ind_asarg.vcf')

noout_gl <- vcfR2genlight(null_set_vcf)

null_set_gl <- gl.drop.loc(noout_gl, outlier)
null_set_gl

# Removing outliers from VCF (important to run BayeAss later on)

null_snps = which(null_set_gl@loc.names %in% noout_gl@loc.names)
null_set_vcf = null_set_vcf[null_snps, ]

####-------------------------- Last NA filtering ---------------------------####

#resetting popmap

geo_data = read.csv("./snp_filtering/infopop.csv", sep = ";") 
genomic_order = as.data.frame(colnames(null_set_vcf@gt)[-1]) 
genomic_order

names(genomic_order) = "Sample_ID"


geo_data_reorder = semi_join(geo_data, genomic_order, join_by(x$id == y$Sample_ID))
geo_data_reorder = inner_join(genomic_order, geo_data_reorder, join_by(x$Sample_ID == y$id))


popmap = geo_data_reorder[,1:2]#[,c(2,1)] #col "pop" as population/location.
names(popmap) = c("id", "pop")


popmap = popmap[popmap$id %in% colnames(null_set_vcf@gt),]

identical(as.character(colnames(null_set_vcf@gt)[-1]), as.character(geo_data_reorder$Sample_ID))
setdiff(as.character(colnames(null_set_vcf@gt)[-1]), as.character(geo_data_reorder$Sample_ID))
setdiff(as.character(geo_data_reorder$Sample_ID), as.character(colnames(null_set_vcf@gt)[-1]))

# Filtering missing data

missing_by_sample(vcfR=null_set_vcf, popmap = popmap)  #verify statistics to take decision
null_set_final_vcf = missing_by_sample(vcfR=null_set_vcf, cutoff = .05)


missing_by_snp(vcfR=null_set_final_vcf)  #verify statistics to take decision
null_set_final_vcf = missing_by_snp(vcfR=null_set_final_vcf, cutoff = .75)

null_set_final_vcf

null_final <- vcfR2genlight(null_set_final_vcf)


####----------------------------- Saving NULL set --------------------------####

null_final <- gl.compliance.check(null_final)

gl2bayescan(null_final, outfile = 'asarg_final_convertion.txt')

vcfR::write.vcf(null_set_final_vcf, file = 'null_set_final.vcf')


# Now that we have a putative null set of allels we can proceed with other analysis



#........................... END OF FILTERING PROCESS .........................#

# ......................................

# ...................

# .........

# .....

# ..


#==============================================================================#
######.'.'.'.'.'.'.'. PART 2 - Analysis using NULL data set '.'.'.'.'.'.'.######
#==============================================================================#




#------------------------- Setting the working directory ----------------------#


setwd("D:/gen_pop_anemonia/")


#---------------------------- Setting color scheme ----------------------------#

# Here I've generated a color pallet using https://coolors.co/
# I've chosen to generate a pallete for people with protanomaly (one of the most commmom types of colorblindnes)

colorspop <- c("#5BC0EB","#E6DF0F","#44C10A","#AB6C3F","#FE5F55")


#--------------------------------- Loading VCF --------------------------------#

# Loading the NULL vcf file

asarg_brazil <- read.vcfR("null_set_final.vcf")

# Checking general information about the file

asarg_brazil



#---------------------- Conversions (software version 7.0) --------------------#

# Converting null set vcf to genind

asarg_brazil_gind <- vcfR2genind(asarg_brazil, return.alleles = TRUE)

asarg_brazil_gind


# Converting null set vcf to genlight

asarg_brazil_gl <- vcfR2genlight(asarg_brazil) 

asarg_brazil_gl

#-------------------------- Adding pop, geo, strata ---------------------------#

pop.data <- read.table("snp_filtering/infopop.csv", sep = ";", header = TRUE)
pop.data <- pop.data[pop.data$id %in% asarg_brazil_gl@ind.names,]  

pop.data

# pop

pop(asarg_brazil_gl) <- pop.data$pop
pop(asarg_brazil_gind) <- pop.data$pop

# lat,long

geovec <- pop.data[,c(3,4)]
other(asarg_brazil_gl)$latlong <- pop.data[,c(pop.data$lat,pop.data$lon)]
other(asarg_brazil_gind)$latlong <- pop.data[,c(pop.data$lat,pop.data$lon)]

# stratus

pop.data[1,] #to check headers

strata(asarg_brazil_gl) <- pop.data[ ,c(5,6,7)]
strata(asarg_brazil_gind) <- pop.data[ ,c(5,6,7)]


#-------------- Saving str format to convert to BayesAss format ---------------#



gl2structure(asarg_brazil_gl, outpath = '.', outfile = 'asarg_null_1178_final')

write_delim(popmap, file = 'popmap.txt', delim = 'tab')


####--------------------------- BASIC METRICS -----------------------------#####

asarg_brazil_gl <- gl.compliance.check(asarg_brazil_gl)

######------------------------ Private alleles ---------------------------######


gl.report.pa(asarg_brazil_gl, method = 'one2rest')

gl.report.pa(asarg_brazil_gl)
######--------------------- Loci that are under HWE ----------------------######


hw <- hw.test(asarg_brazil_gind)
sum(hw[,3] < 0.01)


######--------------------------- Basic Stats ----------------------------######

table <- locus_table(asarg_brazil_gind)

table


######------------------------- Allele richness --------------------------######

richness <- hierfstat::allelic.richness(asarg_brazil_gind)

mean(richness$Ar$PE, na.rm = T) # Pernambuco
mean(richness$Ar$BA, na.rm = T) # Bahia
mean(richness$Ar$ES, na.rm = T) # Espirito Santo
mean(richness$Ar$RJ, na.rm = T) # Rio de Janeiro
mean(richness$Ar$SC, na.rm = T) # Santa Catarina



######-------------------------- Ho,He, Inbreeding -----------------------######


# Inbreeding varies from -1 to 1: (-1) excess of heterozigozity, (1) excess of homozigozity (inbreeding)

div_pop = gl.report.heterozygosity (asarg_brazil_gl, method = "pop")
div_pop$pop
div_pop$FIS

#percent of polymorphic loci

(div_pop$polyLoc[1]*100)/1178  #Bahia
(div_pop$polyLoc[2]*100)/1178  #Espirito Santo
(div_pop$polyLoc[3]*100)/1178  #Pernambuco
(div_pop$polyLoc[4]*100)/1178  #Rio de Janeiro
(div_pop$polyLoc[5]*100)/1178  #Santa Catarina


######-------------------------- Fst & Nei' dist -------------------------######


gl.dist.pop(asarg_brazil_gl, method = "nei")

fst <- gl.fst.pop(asarg_brazil_gl, nboots = 1000, percent = 95, nclusters = 8)
a <- as.matrix(as.dist(fst$Fsts))


heatmap.2(a, dendrogram = 'none',col = brewer.pal(7, 'YlOrRd'),
          cellnote = round(a,3), notecex =1.3, notecol = "black",
          trace = 'none', density.info = 'none', key.title = "Fst")


######-------------------------------- IBD -------------------------------######

# Nei's

gl.ibd(asarg_brazil_gl,distance = 'D', coordinates = geovec, permutations = 1000, 
       Dgeo_trans = 'Dgeo/1000', plot_theme = theme_light()) 
  
# Fst

gl.ibd(asarg_brazil_gl, coordinates = geovec, permutations = 1000, 
       Dgeo_trans = 'Dgeo/1000', plot_theme = theme_light()) 


######------------------------------ Summary -----------------------------######

sumario <- summary(asarg_brazil_gind)
sumario$pop.n.all

# Alelles per population




####                            MULTIFACTOR ANALISIS                        ####



######--------------------------------- DAPC -----------------------------######


#Choosing the number of groups

k_clust = find.clusters(asarg_brazil_gl,n.pca = 35, max.n.clust = 20, n.iter=1e6, n.start=1000) # 36 PCs and 10 pops maximum
k_clust$Kstat

#K was already determined by find.clusters function before, so I'm going to use it here

pnw.dapc <- dapc(asarg_brazil_gl, n.pca = 35, k_clust$grp, n.da = 3) # k_clust$grp is where k is stored


## TESTING THE NUMEBR OF INFORMATIVE PCs


# 1 - Alpha score

alpha<-optim.a.score(pnw.dapc)

# 2 - Cross validation

# Converting genlight to genclone

mygenclone <- as.genclone(asarg_brazil_gind)


#this is to check quickly, from that pick an interval for a more serious run

set.seed(999)
asargx <- xvalDapc(tab(mygenclone, NA.method = "mean"), pop(mygenclone)) 
asargx[-1]$`Number of PCs Achieving Lowest MSE`

# now a more robust test
set.seed(999)
system.time(asargx <- xvalDapc(tab(mygenclone, NA.method = "mean"), pop(mygenclone),
                               n.pca = 1:17, n.rep = 1000, #set n.pca a reasonable number between 2 and a number above the result of last test
                               parallel = "snow", ncpus = 8L))
asargx[-1]$`Number of PCs Achieving Lowest MSE`

# Since alpha-score and crossentropy gave back similar results (6 and 5), 
# I chose to use the higher (one or two extra pcs will not create that much noise)

## Running DAPC again using k and alfa score

pnw.dapc <- dapc(asarg_brazil_gl, k_clust$grp, n.pca = 6, n.da = 4)

eig_percent <- 100 * pnw.dapc$eig / sum(pnw.dapc$eig)

eig_percent

colorsdapc <- c( "#E6DF0F", "#AB6C3F", "#44C10A","#FE5F55", "#5BC0EB") #Had to rearrange color to keep patterns between PCA and DAPC

scatter(pnw.dapc, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomright", scree.pca = TRUE,
        posi.da = "bottomleft", posi.pca = "topleft",pch = 18, label.inds = TRUE, mstree = TRUE, col = colorsdapc)

table.value(table(pop(asarg_brazil_gl), k_clust$grp), 
            col.lab=paste("Gen.Clust", 1:length(k_clust$size)))


# creating a dataframe fo include as strata on gl

dapc_grp <- as.data.frame(pnw.dapc$assign)
names(dapc_grp) = 'dapc'
dapc_grp

addStrata(asarg_brazil_gl) <- dapc_grp

asarg_brazil_gl@strata

######-------------------------------- PCA -------------------------------######

# PCA eigenvalues plot

asarg.pca <- glPca(asarg_brazil_gl, nf = 6, returnDotProd = T) 
eigenvlues_plot <- barplot(100*asarg.pca$eig/sum(asarg.pca$eig), col = heat.colors(39), main="PCA Eigenvalues") +
  title(ylab="Percent of variance\nexplained", line = 2) +
  title(xlab="Eigenvalues", line = 1)

# PCA pc plot

eigenval <- asarg.pca$eig

# variance of each PC

for (i in eigenval) {
print(i/sum(eigenval)*100)
}

12.72+10.83+10.63

# coloring by pop
asarg.pca.scores_pop <- as.data.frame(asarg.pca$scores)


asarg.pca.scores_pop$pop <- asarg_brazil_gl@pop

pca_plot_12 <- ggplot(asarg.pca.scores_pop, aes(x=PC1, y=PC2, colour=pop)) + 
  geom_point(size=2.3) + 
  stat_ellipse(level = 0.95, linewidth = 0.6, alpha = 0.4) + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  color_palette(palette = colorspop) +
  theme_bw()

pca_plot_13 <- ggplot(asarg.pca.scores_pop, aes(x=PC1, y=PC3, colour=pop), color_palette(colorspop)) + 
  geom_point(size=2.3) + 
  stat_ellipse(level = 0.95, linewidth = 0.6, alpha = 0.4) + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  color_palette(palette = colorspop) +
  theme_bw()

pca_plot_23 <- ggplot(asarg.pca.scores_pop, aes(x=PC2, y=PC3, colour=pop), color_palette(pallet = colorspop)) + 
  geom_point(size=2.3) + 
  stat_ellipse(level = 0.95, linewidth = 0.6, alpha = 0.4) + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  color_palette(palette = colorspop) +
  theme_bw()


pca_bypop <- ggarrange(pca_plot_12,pca_plot_13,pca_plot_23, ncol = 2, nrow = 2)

pca_bypop


# coloring by strata

# bin

bin <- as.data.frame(asarg.pca$scores)
bin$pop <- asarg_brazil_gl@strata$bin

pca_bin_plot <- ggplot(bin, aes(x=PC1, y=PC2, colour=pop)) + 
  geom_point(size=2.3) + 
  stat_ellipse(level = 0.95, linewidth = 0.6, alpha = 0.4) + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  theme_bw()

#dapc

dapc <- as.data.frame(asarg.pca$scores)
dapc$pop <- asarg_brazil_gl@strata$dapc

pca_dapc_plot <- ggplot(dapc, aes(x=PC1, y=PC2, colour=pop)) + 
  geom_point(size=2.3) + 
  stat_ellipse(level = 0.95, linewidth = 0.6, alpha = 0.2) + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  theme_bw()

dapc_bin_pop <- ggarrange(pca_plot_12,pca_dapc_plot,pca_bin_plot, ncol = 2, nrow = 2)

dapc_bin_pop

#                     [ [ [ Saving plots as tiff images ] ] ]                  #

setwd("D:/gen_pop_anemonia/images")

tiff(file="IBD.tiff",width=5, height=4, units="in", res=300)
gl.ibd(asarg_brazil_gl, coordinates = geovec, permutations = 1000, 
       Dgeo_trans = 'Dgeo/1000', plot_theme = theme_light())

dev.off()

tiff(file="IBD_Nei.tiff",width=5, height=4, units="in", res=300)
gl.ibd(asarg_brazil_gl,distance = 'D', coordinates = geovec, permutations = 1000, 
       Dgeo_trans = 'Dgeo/1000', plot_theme = theme_light()) 
dev.off()


tiff(file="Fst_heatmap.tiff",width=10, height=8, units="in", res=300)
heatmap.2(a, dendrogram = 'none',col = brewer.pal(7, 'YlOrRd'),
          cellnote = round(a,3), notecex =1.3, notecol = "black",
          trace = 'none', density.info = 'none', key.title = "Fst")
dev.off()


tiff(file="PCA.tiff",width=10, height=8, units="in", res=300)
pca_bypop
dev.off()

tiff(file="PCA_by strata.tiff",width=10, height=8, units="in", res=300)
dapc_bin_pop
dev.off()

tiff(file="Eigenvalues_PCA.tiff",width=10, height=8, units="in", res=300)
eigenvlues_plot <- barplot(100*asarg.pca$eig/sum(asarg.pca$eig), 
                           col = heat.colors(39), main="PCA Eigenvalues") +
  title(ylab="Percent of variance\nexplained", line = 2) +
  title(xlab="Eigenvalues", line = 1)
dev.off()


tiff(file="BIC.tiff",width=10, height=8, units="in", res=300)
k_clust = find.clusters(asarg_brazil_gl,n.pca = 35, max.n.clust = 20, n.iter=1e6, n.start=1000) 
dev.off()

tiff(file="Alpha_scores.tiff",width=10, height=8, units="in", res=300)
optim.a.score(pnw.dapc)
dev.off()

tiff(file="DAPC_crossentropy.tiff",width=10, height=8, units="in", res=300)
asargx$
dev.off()

tiff(file="DAPC.tiff",width=10, height=8, units="in", res=300)
scatter(pnw.dapc, cex = 2, ldw =3, legend = TRUE, clabel = F, posi.leg = "bottomright", scree.pca = TRUE,
        posi.da = "bottomleft", posi.pca = "topleft",pch = 18, label.inds = TRUE, mstree = TRUE, col = colorsdapc)
dev.off()


tiff(file="DAPC_distributions.tiff",width=10, height=8, units="in", res=300)
table.value(table(pop(asarg_brazil_gl), k_clust$grp), 
            col.lab=paste("Gen Clust", 1:length(k_clust$size)))


dev.off()



#.......................... END OF INITIAL ANALYSIS ...........................#

# ......................................

# ...................

# .........

# .....

# ..

#==============================================================================#
####.'.'.'.'.'.'.'.'.'. EXTRA INFO AND OPTIONAL FILTERINGs .'.'.'.'.'.'.'.'.####
#==============================================================================#

# -------------------------------- PLINK 2.0 --------------------------------- #

# clone-checking was done using PLINK2 on "raw" VCF (previous to first filtering)
# plink2 --make-king-table --king-table-filter 0.345 
#check https://www.cog-genomics.org/plink/2.0/distance#make_king for more details


# ---------------------------------- VCFTOOLS -------------------------------- #

# First VCF filtering was done using VCFTOOLs, 

#vcftools --vcf /scratch/genomics/brandaor/projetos/anemonia/anemonia_moufpe_nat_gskimming/results/gen_pop/out/genotyp$         --keep keep_snp_filtering_brazil.txt \

#--min-alleles 2 \
#--max-alleles 2 \
#--thin 1000 \
#--remove-indels \
#--min-meanDP 10 \
#--max-meanDP 150 \
#--max-non-ref-af 0.99 \
#--recode \


# -------------------------------- BayeScan ---------------------------------- #


# Chain convergence was verified using this tutorial as reference: 

# https://theoreticalecology.wordpress.com/2011/12/09/mcmc-chain-analysis-and-convergence-diagnostics-with-coda-in-r/

# Prior odd value was chosen based on manual and this tutorial:

# https://evomics.org/wp-content/uploads/2016/01/BayeScan_BayeScEnv_exercises.pdf



# ---------------------------- Filtering reference ----------------------------#

# Some of the filters used here were adapted from the author bellow

#  Dalapicolla 2023 available in https://github.com/jdalapicolla/PopGenPipe/tree/main/Filtering

# the PCAdpt code was also taken from there. It is an alternative to BayeScan



#----------------------------------- PCAdapt ----------------------------------#


## Removing SNPs possibly under selection using PCAdapt (BayeScan is another option for that)


# Convert VCF to Genotypes #Converter para Genotypes

asarg_brazil
genotypes= t(extract.gt(asarg_brazil, element = "GT", mask = FALSE, 
                        as.numeric = TRUE, return.alleles = FALSE, 
                        IDtoRowNames = TRUE, extract = TRUE, convertNA = TRUE))

#check:

genotypes[1:5,1:5]
row.names(genotypes)

# Convert Genotypes to LFMM #Converter para LFMM

dim(genotypes)
((sum(is.na(genotypes)))/(dim(genotypes)[1]*dim(genotypes)[2]))*100  #Amount of missing data 

genotypes[is.na(genotypes)] = 9 #The missing genotypes have to be encoded with the value 9 or -9 #Os gen?tipos faltantes devem ser codificados com o valor 9 ou -9
genotypes[1:10,1:10]
write.lfmm(genotypes,"asarg_brazil.lfmm")

# Load LFMM

lfmm_input = read.lfmm("asarg_brazil.lfmm")
class(lfmm_input)

# Convert LFMM to a PCadapt matrix and run the analyses.

x <- read.pcadapt(lfmm_input, type="lfmm")
y <- pcadapt(x, K=20)

# With names

poplist.names <- c(rep("PE", 8), rep("ES", 7), rep("RJ", 8), rep("BA", 7), rep("SC", 6))
print(poplist.names)

#"Cattell's Rule" for interpreting the scree plot (PC to the left of the flat line)
#"Regra de Cattell" para interpretar o gr?fico de scree (PC ? esquerda da linha plana)

plot(y, option="screeplot") #2 PCs
plot(y, option="scores", pop = poplist.names) #PC1 and PC2 #2 clusters
plot(y, option="scores", pop = poplist.names, i = 1, j= 3) #PC1 and PC3 #2 clusters
plot(y, option="scores", pop = poplist.names, i = 2, j= 3) #PC2 and PC3 #3 clusters

# With PC2 we see some clustering pattern 

# K = 3. Run PCadapt again. Test K and see the p-values graphs #Rodamos novamente fixando o melhor valor de PCs que ? 3

gen.pcadapt = read.pcadapt(lfmm_input, type = c("lfmm"))
class(gen.pcadapt)
pcadapt.test = pcadapt(gen.pcadapt, K=3, min.maf=0.01, method="mahalanobis")
summary(pcadapt.test)

# Use Mahalanobis distance to compute the (rescaled) test statistics (z-scores in this case).
# The robust Mahalanobis distance is a metric that identifies outliers in multidimensional space. "Robust" means the estimation is not sensitive to outliers in the covariance matrix of the z-scores.


# Graphical tools. #gr?ficos para visualiza??o. Tem muito outliers

plot(pcadapt.test, option = "manhattan")
plot(pcadapt.test, option = "qqplot")
hist(pcadapt.test$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(pcadapt.test, option = "stat.distribution")

# Choosing a cutoff for outlier detection

pcadapt.test$gif #   1.470604

# The GIF indicates how well the test is "calibrated".
# It corrects for inflation of the test score at each locus, which can occur when population
# structure or other confounding factors are not appropriately accounted for in the model.
#
# GIF of 1=well calibrated, >1 =liberal (too many small p-values), <1=conservative (too few small p-values) # Note: GIFs > 2 indicate poor test calibration.
#
# For a given alpha (real valued number between 0 and 1), SNPs with q-values less than alpha will be considered as outliers with an expected false discovery rate bounded by alpha. The false discovery rate is defined as the percentage of false discoveries among the list of candidate SNPs. Here is an example of how to provide a list of candidate SNPs, for an expected false discovery rate lower than 10%.


# q-values

qval = qvalue(pcadapt.test$pvalues)$qvalues
alpha = 0.1
outliers1 = which(qval < alpha)
length(outliers1) #It will be eliminated 3769 SNPs 

#Benjamini-Hochberg Procedure

padj = p.adjust(pcadapt.test$pvalues, method="BH")
alpha = 0.1
outliers2 = which(padj < alpha)
length(outliers2) #Number of SNPs that will be eliminated

#Bonferroni correction

padj2 = p.adjust(pcadapt.test$pvalues, method="bonferroni")
alpha = 0.1
outliers3 = which(padj2 < alpha)
length(outliers3) #Number of SNPs that will be eliminated


#Choose one approach to eliminate outlier SNPs. In the case Benjamini-Hochberg Procedure the "outliers2"

#Which approach to chose???

#Remove SNPs with adaptive signals


snps_to_keep = which(!c(1:dim(asarg_brazil@gt)[1]) %in% outliers2)
asarg_brazil = asarg_brazil[snps_to_keep, ]
asarg_brazil



# Remove by missing data for 20%

missing_by_sample(vcfR=asarg_brazil, popmap = popmap)  #verify statistics #verificar as estat?sticas
asarg_brazil = missing_by_sample(vcfR=asarg_brazil, cutoff = .2)


missing_by_snp(vcfR=asarg_brazil)  #verify statistics #verificar as estat?sticas
asarg_brazil = missing_by_snp(vcfR=asarg_brazil, cutoff = .8)
asarg_brazil



# ------------------------------- DAPC --------------------------------------- #

# Some info for DAPC I took from this tutorial

# https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html



# -------------------- Migration rate matrix (BAYESASS) ---------------------- #


setwd('D:/BayesAss/BayesAss3-SNPs-master')

nb.cols <- 80
mycolors <- colorRampPalette(brewer.pal(8, "Blues"))(nb.cols)

migration <-  as.matrix(read.delim("migration_rate_matrix1.txt" , row.names= 1, sep=",")) 

setwd('D:/gen_pop_anemonia/Images')

tiff(file="Heatmap_migration.tiff",width=5, height=4, units="in", res=300)
heatmap.2(migration, dendrogram = 'none', notecex =1, notecol = "black", 
          Colv = NA, Rowv = NA, cellnote = migration, trace = 'none', col = mycolors, 
          density.info = 'none', key.title = "m", scale = 'none', breaks=seq(0.01,0.09,0.001))
dev.off()
