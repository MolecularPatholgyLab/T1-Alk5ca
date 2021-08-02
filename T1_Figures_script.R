####################################################################################################################################################
##################################################    Script of figures for   ######################################################################
# Epithelial TGFβ activation engages growth-factor signalling to circumvent cell death and drive intestinal tumorigenesis with aggressive features #
###########################---------------------------------------------------------------------------------------##################################
##############################################      Authour: Raheleh Amirkhah       ################################################################
##############################################    Email: r.amirkhah@qub.ac.uk       ################################################################
####################################################################################################################################################
###############################################################################################
#                                                                                             #
# Figure 1A : Lines 88 to 128                                                                 #
# Figure 1B : Lines 135 to 157                                                                #
# Figure 1C : Lines 161 to 215                                                                #
# Figure 1D, E, G  & Extended Data Figure 1C : Lines 219 to 250                               #
# Figure 1F : Lines 251 to 281                                                                #
# Extended Data Figure 1A : Lines 288 to 315                                                  #
# Extended Data Figure 1B : Lines 317 to 482                                                  #
# Extended Data Figure 1D : Lines 485 to 561                                                  #
# Figure 6A : Lines 568 to 594                                                                #
# Figure 6B : Lines 598 to 624                                                                #
# Figure 6D : Lines 628 to 683                                                                #
# Figure 6E : Lines 687 to 717                                                                #
# Figure 6F : Lines 723 to 748                                                                #
# Figure 6G and Extended Data Figure 7d : Lines 752 to 826                                    #
# Figure 6H : Lines 828 to 904                                                                #
# Figure 6I : Lines 912 to 936                                                                #
# Figure 7F : Lines 939 to 978                                                                #
# Extended Data Figure 7A : Lines 985 to 1024                                                 #
# Extended Data Figure 7B : Lines 1029 to 1065                                                #
# Extended Data Figure 7C : Lines 1068 to 1102                                                #


###############################################################################################

## Load libraries
library(MCPcounter)
library(ggsignif)
library(ggplot2)
library(survival)
library(survminer)
library(survivalAnalysis)
library(magrittr)
library(fgsea)
library(GSVA)
library(msigdbr)
library(tibble)
library(plyr)
library(dplyr)
library(pheatmap)
library(grid)
library(biomaRt)
library(reshape)
library(RColorBrewer)
library(gplots)
options(stringsAsFactors = F)


##################################################################################################################################################
##################################################### Import T1 expression data ##################################################################
## read T1 data
t1 <- read.csv("T_cohort/T1_collapsed.csv")
head(t1[1:10])
t1 <- column_to_rownames(t1, var = "group")[-1]

## read meta file
meta <- read.csv("T1_cohort/T1_meta.csv")

## check if samples in T1 and meta file are in the same order
all(colnames(t1)==meta$Sample_ID)


###################################################################################################################################################
################################### Estimate Microenvironment Cell Population using  MCP-counter  #################################################
##################################################--------------------------------------###########################################################

MCP <- MCPcounter.estimate(t1,featuresType="HUGO_symbols",
                           probesets=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt"),sep="\t",stringsAsFactors=FALSE,colClasses="character"),
                           genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE))

###  Transpose the result
MCP <-as.data.frame(t(MCP))
MCP$group <- meta$group

###################################################################################################################################################
########################################################### Figure 1 ############################################################################
####################################################################################################################################################

################################### Figure 1A #############################################

## wilcox.test to test significant differencce of fibroblast population between relapse and non-relapse group 
wilcox.test(MCP[MCP$group=="Relapse","Fibroblasts"], MCP[MCP$group=="Non-Relapse","Fibroblasts"])

## plot
tiff("fibroblast_new.tiff", units="mm", width=54, height=54, res=300)
set.seed(121)
p_f <- ggplot(MCP, aes(x = group, y = Fibroblasts, fill = group))+geom_boxplot(size = 0.3)+geom_jitter(size=0.65, width=.1)+
  labs(y = "MCP-counter scores")+theme(legend.position = "none", axis.title.x=element_blank(), 
                                       axis.title.y = element_text(size = 12, family = "Arial"), axis.text.x = element_text(size = 10, family = "Arial"),
                                       plot.title = element_text(hjust = 0.5, size = 12, family = "Arial"),axis.text.y = element_text(size = 8),
                                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), 
                                       axis.line = element_line())+ggtitle("Fibroblasts")+scale_fill_manual(values = c('Relapse' = '#ff4242', 'Non-Relapse' = '#1f1fd8'))+
  ylim(c(8,10.5))

#### Add p-value (more beautiful)
p_f+geom_signif(comparisons = list(c("Relapse", "Non-Relapse")), test="wilcox.test",
                map_signif_level=TRUE, color = "black", annotations = 'P=0.17', fontface = "italic")
dev.off()

#############--------------------------------------------------

## wilcox.test to test significant differencce of Cytotoxic lymphocytes population between relapse and non-relapse group 
wilcox.test(MCP[MCP$group=="Relapse", "Cytotoxic lymphocytes"], MCP[MCP$group=="Non-Relapse", "Cytotoxic lymphocytes"])

## plot
tiff("CTL_new.tiff", units="mm", width=54, height=54, res=300)
set.seed(121)
p_f <- ggplot(MCP, aes(x = group, y = `Cytotoxic lymphocytes`, fill = group))+geom_boxplot(size = 0.3)+geom_jitter(size=0.65, width=.1)+
  labs(y = "MCP-counter scores")+theme(legend.position = "none", axis.title.x=element_blank(), 
                                       axis.title.y = element_text(size = 12, family = "Arial"), axis.text.x = element_text(size = 10, family = "Arial"),
                                       plot.title = element_text(hjust = 0.5, size = 12, family = "Arial"),axis.text.y = element_text(size = 8),
                                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), 
                                       axis.line = element_line())+ggtitle("Cytotoxic lymphocytes")+scale_fill_manual(values = c('Relapse' = '#ff4242', 'Non-Relapse' = '#1f1fd8'))+
  ylim(c(3.5,4.3))

#### Add p-value (more beautiful)
p_f+geom_signif(comparisons = list(c("Relapse", "Non-Relapse")), 
                map_signif_level=TRUE, color = "black", annotations = 'P=0.57', fontface = "italic")
dev.off()



################################### Figure 1B #############################################

## read QuPath score
QuPath <- read.delim("T1_cohort/scort id paired with stroma and tumour epithelium scores.txt")

### join fibroblast score from MCP with digi-patological score
mcp_qpath <- cbind(MCP[,c(10,11)], QuPath_stromal = QuPath$Stroma.fibroblast.score..pixel.)

### measure correlation between transcriptome-based stromal scores and the digital pathology-based stromal scores
cor.test(mcp_qpath$Fibroblasts, mcp_qpath$QuPath_stromal, method = "spearman")


tiff("cor_mcp_qpath_2.tiff", units="mm", width=75, height=55, res=300)

g <- ggplot(mcp_qpath, aes(Fibroblasts, QuPath_stromal))+geom_point(size = 0.65, aes(color=group))+geom_smooth(method = "lm", size=0.2)+
  labs(x = "MCP-counter scores", y = "QuPath-scores")+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
                   panel.border = element_blank(),axis.line = element_line(), legend.position= c(0.87, 0.2),legend.text = element_text(size = 7),
                   legend.title = element_blank(),legend.key.width = unit(0.002, "mm"),legend.key.size = unit(0.3, "cm"),
                   axis.title = element_text(size = 12, family = "Arial"),plot.title = element_text(hjust = 0.5, size = 12, 
                                                                                                    family = "Arial"))+ggtitle("Fibroblasts")

g+annotate("text", x=8.5, y=99, label="rho = 0.69", size = 3.5)+annotate("text", x=8.5, y=90, label="P = 0.00009", size = 3.5, fontface = "italic")+
  scale_color_manual(values = c('Relapse' = '#ff4242', 'Non-Relapse' = '#1f1fd8'))

dev.off()



################################### Figure 1C #############################################

## read T1 mutation profile
mutationData_t1 <- read.delim("T1_cohort/T1_mutation_for.custom.csv", stringsAsFactors = F, sep = ",")

selected_mutationData_t1 <- mutationData_t1[which(mutationData_t1$gene %in% c("APC", "TP53", "KRAS", "SMAD4")),]
colnames(selected_mutationData_t1)[2] <- "Mutation"


### keep unique
selected_mutationData_t1_2 <- selected_mutationData_t1[,-c(4,5)]
selected_mutationData_t1_3 <- unique(selected_mutationData_t1_2)
df <- selected_mutationData_t1_3 %>%
  group_by(group, Mutation) %>%
  summarise(n())

## caculate proportion
df <- as.data.frame(df)
df$Proportion <- ifelse(df$group=="non-Relapse", (df$n*100)/13, (df$n*100)/8)


## plot
tiff("Figure 1C.tiff", units="mm", width=40, height=55, res=300)
p <- ggplot(data=df, aes(x=Mutation, y=Proportion, fill=group)) +
  geom_bar(stat="identity", position=position_dodge())+xlab(" ") + ylab("Mutation Proportions")+ 
  scale_y_continuous(expand = c(0, 0))+
  ## Plot background theme
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 7, family = "Helvetica", angle = 45, margin = margin(t = 7, r = 20, b = 0, l = 0)),
        axis.title.y = element_text(size = 12, family = "Helvetica"),
        axis.text.y = element_text(size = 6, colour = "black", angle = 90, hjust = 0.5),
        axis.ticks.y = element_line(size = 0.5, colour = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.ticks.x = element_line(size = 0.5, colour = "black"),
        panel.background = element_rect(fill = 'white'),
        
        ##legend
        legend.background = element_blank(),
        #legend.text = element_text(size = 18),
        legend.title = element_blank(),
        legend.title.align = 0,
        legend.position = "top",
        legend.key.size = unit(0.1, "cm"),
        legend.margin = margin(c(5, 5, 5, 0)),
        legend.box.margin=margin(-7,-10,-10,-10),
        legend.text = element_text(margin = margin(r = 0.15, unit = "cm"), size=6),
        legend.spacing.x = unit(0, 'cm'))+
  scale_fill_manual(values = c('Relapse' = '#ff4242', 'non-Relapse' = '#1f1fd8'))
p
dev.off()



################################# Figure 1D, E, G  & Extended Data Figure 1C  ######################

## read differentially expressed genes (relapse VS non-relapse)
DEG <-read.delim("T1_cohort/DEGs_Relapse-Nonrelapse.txt")
head(DEG)
## provide a ranked gene list based on log FC
gene_list <- as.vector(DEG$logFC)
names(gene_list) <- row.names(DEG)
gene_list <- sort(gene_list, decreasing = T)

################ load all genesets
load("T1_cohort/all_genesets.RData")  #all_genesets


############# pre-ranked GSEA 
set.seed(127)
fgsea_re <- fgsea(all_genesets, stats = gene_list , 1000, minSize = 1, maxSize = Inf, nproc = 0, gseaParam = 1, BPPARAM = NULL)

## plot gsea result for MTORC1 geneset
tiff("Fig.1D_MTORC1.tiff", units="mm", width=55, height=40, res=300)
plotEnrichment(all_genesets[["HALLMARK_MTORC1_SIGNALING"]], gene_list, gseaParam = 1, ticksSize = 0.1)+labs(title="HALLMARK_MTORC1")+
  labs(y = "Enrichment scores")+
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", size = 12), panel.background = element_blank(), 
        axis.title.y = element_text(size = 12, family = "Arial"),axis.text.y = element_text(size = 8, family = "Arial"),
        panel.border = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
grid.text(label = "NES = 2.26\npadj = 0.002", gp = gpar(fontsize = 8), x = unit(0.8, "npc"), y = unit(0.68, "npc"))
grid.text(label = "Relapse", gp = gpar(fontsize = 8), x = unit(0.3, "npc"), y = unit(0.04, "npc"))
grid.text(label = "non-Relapse", gp = gpar(fontsize = 8), x = unit(0.85, "npc"), y = unit(0.04, "npc"))
dev.off()



########################################### Figure 1F ###########################################

### Get Hallmark genesete from msigdb
H.df <- msigdbr(species = "Homo sapiens", category = "H")
H_list <- H.df %>% split(x = .$gene_symbol, f = .$gs_name)

### ssGSEA --- GSVA
set.seed(127)
ssgsea_hallmark <- gsva(as.matrix(t1), H_list, 
                        min.sz=10, max.sz=Inf, verbose = T, method = 'ssgsea')


### Combine TGFb ssGSEA score with MCP data
MCP$HALLMARK_TGF_BETA_SIGNALING <- ssgsea_hallmark["HALLMARK_TGF_BETA_SIGNALING",]

### measure correlation between transcriptome-based stromal scores and TGFb ssGSEA score
cor.test(MCP$Fibroblasts, MCP$HALLMARK_TGF_BETA_SIGNALING, method = "spearman")

## plot
tiff("Fig.1F.tiff", units="mm", width=75, height=55, res=300)
p <- ggplot(MCP, aes(Fibroblasts, HALLMARK_TGF_BETA_SIGNALING))+geom_point(size = 0.65, aes(color=group))+geom_smooth(method = "lm", size=0.2)+
  labs(x = "MCP-score_Fibroblasts", y = "ssGSEA-score_TGFb")+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
                   panel.border = element_blank(),axis.line = element_line(), axis.text = element_text(size = 8, family = "Arial"),
                   legend.position= c(0.87, 0.15),legend.text = element_text(size = 6),legend.title = element_blank(),legend.key.width = unit(0.002, "mm"),
                   legend.key.size = unit(0.3, "cm"),axis.title = element_text(size = 12, family = "Arial"))

p+annotate("text", x=8.5, y=0.63, label="rho = -0.07", size = 3.5)+annotate("text", x=8.5, y=0.615, label="P = 0.7", size = 3.5, fontface = "italic")+
  scale_color_manual(values = c('Relapse' = '#ff4242', 'Non-Relapse' = '#1f1fd8'))

dev.off()


#######################################################################################################################################################
#########################################################  Extended Data Figure 1  ####################################################################
#######################################################################################################################################################

####################################### Extended Data Figure 1A ##########################

## read the T1 normalized data 
exp <- read.delim("T1_cohort/T1_normalized.txt")
exp <- column_to_rownames(exp, var = "ID")

## Row-wise centring expression data 
exp.scale <- t(scale(t(exp), scale = F))

####PCA analysis
pc <- prcomp(exp.scale, scale. = F, center = F)
pc$rotation
pcr <- data.frame(pc$rotation[,1:2], Group = meta$group)
pcr$label <- meta$label

## plot
tiff("Extended Data Figure 1A.tiff", units="mm", width=60, height=50, res=300)
p=ggplot(pcr, aes(PC1, PC2, color = Group)) + geom_point(size=0.65)+labs(x='PC1=31.82%', y='PC2=8.91%')+
  theme_bw(base_size = 7)+theme(axis.text.x = element_text(size = 7.5),legend.key.width = unit(2, "mm"),
                                axis.text.y = element_text(size = 7.5),axis.title = element_text(size = 9.5, family = "Helvetica"),
                                legend.text = element_text(margin = margin(r = 0.15, unit = "cm"), size=7),
                                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                legend.margin=margin(5,0,0,0),legend.box.margin=margin(-10,-10,-7,-10),
                                panel.background = element_blank(), legend.title = element_blank(), legend.position = "top",
                                legend.spacing.x = unit(0, 'cm'))

p+scale_color_manual(values = c('Relapse' = '#ff4242', 'Non-Relapse' = '#1f1fd8'))
dev.off()

####################################### Extended Data Figure 1B ##########################

## read CRIS classification result
cris <- read.delim("T1_cohort/CRIS_predictions.txt")
cris$label <- meta$group
cris <- cris[,c(2, 8)]

cris[,3] <- ifelse(cris$predict.label==1, "CRIS-A",
                   ifelse(cris$predict.label==2, "CRIS-B",
                          ifelse(cris$predict.label==3, "CRIS-C",
                                 ifelse(cris$predict.label==4, "CRIS-D",
                                        ifelse(cris$predict.label==5, "CRIS-E", NA)))))

table(cris$V3)

######### CRIS Relapse
cris.relapse <- data.frame(Group = c("CRIS-A", "CRIS-B", "CRIS-C", "CRIS-D", "CRIS-E"), Relapse = c(30,10,10,20,30))
cris.relapse <- melt(cris.relapse)

## save plot
tiff("Extended Data Figure 1B_CRIS_relapse.tiff", units="mm", width=40, height=40, res=300)

pie = ggplot(cris.relapse, aes(x="", y=value, fill=Group)) + geom_bar(stat="identity", width=1)
# Convert to pie (polar coordinates) and add labels
pie = pie + coord_polar("y", start=0) + geom_text(aes(label = paste0(value, "%")), 
                                                  position = position_stack(vjust = 0.5), size = 2.5)
# Add color scale (hex colors)
pie = pie + scale_fill_manual(values=c("#ef4c2a", "#d42127", "#1c275c", "#018647", "#00ad9b")) 
# Remove labels and add title
pie = pie + labs(x = NULL, y = NULL, fill = NULL, title = "Relapse")
# Tidy up the theme
pie + theme_classic() + theme(axis.line = element_blank(),
                              axis.text = element_blank(),
                              axis.ticks = element_blank(),
                              plot.title = element_text(hjust = 0.5, color = "black", size = 12),
                              legend.text = element_text(margin = margin(r = 0.15, unit = "cm"), size=6),
                              legend.spacing.x = unit(0, 'cm'),
                              legend.key.height = unit(3, "mm"),legend.key.width = unit(2, "mm"),
                              legend.box.margin=margin(-15,5,-5,15),
                              axis.ticks.length = unit(0, "mm"),
                              #axis.title=element_text(size=14,face="bold"),
                              #panel.grid.major = element_blank(),
                              panel.spacing.x=unit(0, "lines"),
                              panel.spacing.y=unit(0, "lines"),
                              legend.position="right",
                              legend.box.spacing = unit(-0.8, 'cm'),
                              legend.margin = margin(0, 0, 0, 0, "cm"),
                              legend.title = element_blank(),
                              plot.caption = element_text(hjust = 0,margin = unit(c(-5,0,0,0), "mm")),
                              plot.margin = margin(0, 0, -0.7, -0.3, "cm"))

dev.off()


########### CRIS  Non-relapse
cris.nonrelapse <- data.frame(Group = c("CRIS-A", "CRIS-B", "CRIS-C", "CRIS-D", "CRIS-E"), Relapse = c(3,2,5,3,4))
cris.nonrelapse <- melt(cris.nonrelapse)

#### Save plot
tiff("Extended Data Figure 1B_CRIS_nonrelapse.tiff", units="mm", width=40, height=40, res=300)

pie = ggplot(cris.nonrelapse, aes(x="", y=value, fill=Group)) + geom_bar(stat="identity", width=1)
# Convert to pie (polar coordinates) and add labels
pie = pie + coord_polar("y", start=0) + geom_text(aes(label = paste0(round((value*100)/17), "%")), position = position_stack(vjust = 0.5), size =2.5)
# Add color scale (hex colors)
pie = pie + scale_fill_manual(values=c("#ef4c2a", "#d42127", "#1c275c", "#018647", "#00ad9b")) 
# Remove labels and add title
pie = pie + labs(x = NULL, y = NULL, fill = NULL, title = "non-Relapse")
# Tidy up the theme
pie + theme_classic() + theme(axis.line = element_blank(),
                              axis.text = element_blank(),
                              axis.ticks = element_blank(),
                              plot.title = element_text(hjust = 0.5, color = "black", size = 12),
                              legend.text = element_text(margin = margin(r = 0.15, unit = "cm"), size=6),
                              legend.spacing.x = unit(0, 'cm'),
                              legend.key.height = unit(3, "mm"),legend.key.width = unit(2, "mm"),
                              legend.box.margin=margin(-15,5,-5,15),
                              axis.ticks.length = unit(0, "mm"),
                              panel.spacing.x=unit(0, "lines"),
                              panel.spacing.y=unit(0, "lines"),
                              legend.position="right",
                              legend.box.spacing = unit(-0.8, 'cm'),
                              legend.margin = margin(0, 0, 0, 0, "cm"),
                              legend.title = element_blank(),
                              plot.caption = element_text(hjust = 0,margin = unit(c(-5,0,0,0), "mm")),
                              plot.margin = margin(0, 0, -0.7, -0.3, "cm"))
dev.off()


##########--------------------------------------------------------------------------
## read CMS classification result using RF
CMS_RF <- read.delim("T1_cohort/CMS_predictions.txt")
rownames(CMS_RF) <- CMS_RF$Name
table(CMS_RF$RF)

#### CMS ---- non-Relapse
df.nonRe <- data.frame(Group = c("CMS1", "CMS2", "CMS3", "CMS4", "UC"), nonRelapse = c(1,6,2,2,6))
df.nonRe <- melt(df.nonRe)

## save plot
tiff("Extended Data Figure 1B_CMS_nonrelapse.tiff", units="mm", width=40, height=40, res=300)

pie = ggplot(df.nonRe, aes(x="", y=value, fill=Group)) + geom_bar(stat="identity", width=1)
# Convert to pie (polar coordinates) and add labels
pie = pie + coord_polar("y", start=0) + geom_text(aes(label = paste0(round((value*100)/17), "%")), position = position_stack(vjust = 0.5), size = 2.5)
# Add color scale (hex colors)
pie = pie + scale_fill_manual(values=c("#e69c00", "#0070af", "#c877a5", "#089b73", "#999999")) 
# Remove labels and add title
pie = pie + labs(x = NULL, y = NULL, fill = NULL, title = "non-Relapse")
# Tidy up the theme
pie + theme_classic() + theme(axis.line = element_blank(),
                              axis.text = element_blank(),
                              axis.ticks = element_blank(),
                              plot.title = element_text(hjust = 0.5, color = "black", size = 12),
                              legend.text = element_text(margin = margin(r = 0.15, unit = "cm"), size=6),
                              legend.spacing.x = unit(0, 'cm'),
                              legend.key.height = unit(3, "mm"),legend.key.width = unit(2, "mm"),
                              legend.box.margin=margin(-15,5,-5,15),
                              axis.ticks.length = unit(0, "mm"),
                              panel.spacing.x=unit(0, "lines"),
                              panel.spacing.y=unit(0, "lines"),
                              legend.position="right",
                              legend.box.spacing = unit(-0.8, 'cm'),
                              legend.margin = margin(0, 0, 0, 0, "cm"),
                              legend.title = element_blank(),
                              plot.caption = element_text(hjust = 0,margin = unit(c(-5,0,0,0), "mm")),
                              plot.margin = margin(0, 0, -0.7, -0.3, "cm"))
dev.off()



##### CMS ---- Relapse
df.Re <- data.frame(Group = c("CMS1", "CMS2", "CMS3", "CMS4", "UC"), Relapse = c(0,30,10,20,40))
df.Re <- melt(df.Re)

## save plot
tiff("Extended Data Figure 1B_CMS_relapse.tiff", units="mm", width=40, height=40, res=300)

pie = ggplot(df.Re, aes(x="", y=value, fill=Group)) + geom_bar(stat="identity", width=1)
# Convert to pie (polar coordinates) and add labels
pie = pie + coord_polar("y", start=0) + geom_text(aes(label = paste0(value, "%")), position = position_stack(vjust = 0.5), size = 2.5)
# Add color scale (hex colors)
pie = pie + scale_fill_manual(values=c("#e69c00", "#0070af", "#c877a5", "#089b73", "#999999")) 
# Remove labels and add title
pie = pie + labs(x = NULL, y = NULL, fill = NULL, title = "Relapse")
# Tidy up the theme
pie + theme_classic() + theme(axis.line = element_blank(),
                              axis.text = element_blank(),
                              axis.ticks = element_blank(),
                              plot.title = element_text(hjust = 0.5, color = "black", size = 12),
                              legend.text = element_text(margin = margin(r = 0.15, unit = "cm"), size=6),
                              legend.spacing.x = unit(0, 'cm'),
                              legend.key.height = unit(3, "mm"),legend.key.width = unit(2, "mm"),
                              legend.box.margin=margin(-15,5,-5,15),
                              axis.ticks.length = unit(0, "mm"),
                              panel.spacing.x=unit(0, "lines"),
                              panel.spacing.y=unit(0, "lines"),
                              legend.position="right",
                              legend.box.spacing = unit(-0.8, 'cm'),
                              legend.margin = margin(0, 0, 0, 0, "cm"),
                              legend.title = element_blank(),
                              plot.caption = element_text(hjust = 0,margin = unit(c(-5,0,0,0), "mm")),
                              plot.margin = margin(0, 0, -0.7, -0.3, "cm"))

dev.off()


####################################### Extended Data Figure 1D ##########################

## read FOCUS data
exp <- read.csv("FOCUS_cohort/scort_data_ws2_expression_mean_per_gene_v1.3.csv")
sum(duplicated(exp$Gene))
#exp[exp$Gene=="Mar-01",]
exp <- exp[!duplicated(exp$Gene),]
rownames(exp) <- NULL
##
exp.2 <- column_to_rownames(exp, var = "Gene")
head(exp.2[1:5])

## read meta data
meta <- read.csv("FOCUS_cohort/clin.csv")
meta.2 <- column_to_rownames(meta, var = "Sample_ID")

###### only get part of meta data which have have expression value for that
meta.3 <- meta.2[colnames(exp.2),]
dim(meta.3)

#### check if the samples in meta and expression file are in the same order
all(rownames(meta.3)==colnames(exp.2))

### ssGSEA --- GSVA
set.seed(124)
ssgsea_focus <- gsva(as.matrix(exp.2), H_list, ssgsea.norm=T,
                     min.sz=1, max.sz=Inf, verbose = T, method = 'ssgsea')

t_ssgsea_focus <- as.data.frame(t(ssgsea_focus))
t_ssgsea_focus <- rownames_to_column(t_ssgsea_focus, var = "sample_ID")

### MCP
MCP_FOCUS <- MCPcounter.estimate(exp.2,featuresType="HUGO_symbols",
                                 probesets=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt"),sep="\t",stringsAsFactors=FALSE,colClasses="character"),
                                 genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE))

MCP_FOCUS[1:5]

##
t_MCP_FOCUS <-as.data.frame(t(MCP_FOCUS))
t_MCP_FOCUS <- rownames_to_column(t_MCP_FOCUS, var = "sample_ID")
dim(t_MCP_FOCUS)

################################# JOIN fibroblast, Cytotoxic lymphocytes and TGFb scores in one file
df <- join(t_ssgsea_focus[,c(1,45)], t_MCP_FOCUS[, c(1,4,11)], by="sample_ID", type="inner")
df <- column_to_rownames(df, var="sample_ID")
##
all(rownames(meta.3)==rownames(df))

### calculate correlation between fibroblast and TGFb score
cor.test(df$Fibroblasts, df$HALLMARK_TGF_BETA_SIGNALING)

## plot
tiff("Extended Data Figure 1D_(ssGSEA-TGFb)_(MCP-Fibroblast).tiff", units="mm", width=70, height=60, res=300)
p <- ggplot(df, aes(Fibroblasts, HALLMARK_TGF_BETA_SIGNALING))+geom_point(size = 0.65)+geom_smooth(method = "lm", size=0.2)+
  labs(x = "MCP score_Fibroblast", y = "ssGSEA score_TGFB")+theme_bw()+theme(panel.grid.major = element_blank(), 
                                                                             panel.grid.minor = element_blank(),panel.background = element_blank(),
                                                                             legend.position= c(0.2, 0.85),legend.text = element_text(size = 6),legend.title = element_text(size = 6, face = "bold"),
                                                                             legend.key.size = unit(0.3, "cm"),panel.border = element_blank(),axis.line = element_line(), 
                                                                             axis.text = element_text(size = 6, family = "Arial"),axis.title = element_text(size = 12, family = "Arial"))
p+annotate("text", x=7.3, y=0.4, label="r = 0.55", size = 3)+annotate("text", x=7.3, y=0.38, label="P < 2.2e-16", size = 3, fontface = "italic" )
dev.off()

### calculate correlation between Cytotoxic lymphocytes and TGFb score
cor.test(df$`Cytotoxic lymphocytes`, df$HALLMARK_TGF_BETA_SIGNALING)

## plot
tiff("Extended Data Figure 1D_(ssGSEA-TGFb)_(MCP-Cytotoxic).tiff", units="mm", width=72, height=60, res=300)
p <- ggplot(df, aes(`Cytotoxic lymphocytes`, HALLMARK_TGF_BETA_SIGNALING))+geom_point(size = 0.65)+geom_smooth(method = "lm", size=0.2)+
  labs(x = "MCP score_Cytotoxic lymphocytes", y = "ssGSEA score_TGFB")+theme_bw()+theme(panel.grid.major = element_blank(), 
                                                                                        panel.grid.minor = element_blank(),panel.background = element_blank(),legend.position= c(0.2, 0.85),
                                                                                        legend.text = element_text(size = 6),legend.title = element_text(size = 6, face = "bold"),legend.key.size = unit(0.3, "cm"),
                                                                                        panel.border = element_blank(),axis.line = element_line(), axis.text = element_text(size = 6, family = "Arial"),
                                                                                        axis.title = element_text(size = 12, family = "Arial"))
p+annotate("text", x=2.25, y=0.4, label="r = -0.3", size = 3)+annotate("text", x=2.25, y=0.38, label="P < 1.09e-09", size = 3, fontface = 'italic')

dev.off()


######################################################################################################################################################
#######################################################################  Figure 6 ####################################################################
######################################################################################################################################################

####################################### Figure 6A ##########################
# read the gene list obtained from PAMR
pamr_output <- read.csv("mouse_models/PAMR_result/Alk5_GeneSignature_pamr_thres.4.98.csv")
# read the expression profile of Alk5ca_WT mouse models
Alk5_collapsed <- read.delim("mouse_models/PAMR_result/Alk5ca_WT_mouse_collapsed.txt")
colnames(Alk5_collapsed)[1] <- "symbol"
### join both data to get expression value for pamr gene list
pamr_exp <- join(Alk5_collapsed, pamr_output, by = "symbol", type = "inner")

## heatmap plot
pamr_exp.2 <- column_to_rownames(pamr_exp, var = "symbol")[-c(1,8:10)]
## row-wise centring
pamr_exp.3 <- pamr_exp.2 - rowMeans(pamr_exp.2)


#### provide meta file for heatmap
meta <- data.frame(Genotype=c(rep("WT",3), rep("Alk5ca", 3)))
rownames(meta) <- colnames(pamr_exp.3)

ann_colors = list(Genotype = c(WT = "#E69F00", Alk5ca = "darkgreen"))

tiff("Figure 6A.tiff", units="mm", width=90, height=70, res=300)
my_palette <- colorRampPalette(c("#1f1fd8", "#fffcfc", "#fb4141"))(n = 1000)
pheatmap(pamr_exp.3, color=my_palette, fontsize_row =6, border_color = NA, family = "Arial",
         treeheight_row = 8, cluster_cols = F, legend = T, show_colnames=F,
         annotation_col = meta, annotation_colors=ann_colors)
dev.off()



####################################### Figure 6B ##########################

classifier_re <- matrix(c(1, 9, 13, 4),
                        nrow = 2, dimnames = list( c("WT", "Alk5"),c("Relapse", "nonRelapse")))
test_res <- fisher.test(classifier_re)

## convert the number to percent
df <- data.frame(label = c("WT-like", "Alk5-like", "WT-like", "Alk5-like"), type = c("Relapse", "Relapse", "non-Relapse", "non-Relapse"), value = c(10, 90, 76.5, 23.5))

### plot
tiff("Figure 6B.tiff", units="mm", width=41, height=42, res=300)
p <- ggplot(df, aes(x = type, y = value))+geom_bar(stat = "identity", aes(fill = label))+labs(y = "Percentage of samples")+
  scale_y_continuous(expand = c(0, 0))+ # start plot from 0
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.title = element_blank(), legend.position = "top",
        legend.text = element_text(margin = margin(r = 0.15, unit = "cm"), size=6),
        legend.key.size = unit(0.1, "cm"),
        legend.margin=margin(2,0,0,0),
        legend.box.margin=margin(-7,-10,-10,-10),
        legend.spacing.x = unit(0, 'cm'),
        panel.border = element_blank(), axis.line = element_line(), panel.background = element_blank(), 
        axis.text.x = element_text(size = 7, family = "Arial"), axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 11, family = "Arial"),axis.title.x = element_text(size = 7, family = "Arial")
  )+scale_fill_manual(values = c('Alk5-like' = '#ff4242', 'WT-like' = '#1f1fd8'))+labs(x = "Fisher's exact test = 0.001")
p
#+annotate("text", x=2, y=6, label="Relapse 3", size = 2, color = "white", fontface = "bold")
dev.off()



####################################### Figure 6D ##########################

stage2CRC <- read.delim("stageII_CRC_Sanghee/patient.4clusters_2.csv", sep = ",")
stage2CRC <- stage2CRC[,-c(2:9,11)]
head(stage2CRC)

recsurv.2<-with(stage2CRC, Surv(RFS, Rec.Status))
fit_KM.2 <- survfit(recsurv.2~1,type="kaplan-meier", conf.type="log-log")
summary(fit_KM.2)
cutoff <- 50
stage2CRC$perc.BM_Clust.4.4 <- ifelse(stage2CRC$perc.BM_Clust.4.4 > cutoff, 'High','Low')
stage2CRC$perc.BM_Clust.4.4 <- factor(stage2CRC$perc.BM_Clust.4.4, levels = c('High', 'Low'))


## save plot
tiff("Figure 6D.tiff", units="mm", width=70, height=50, res=300)
k <- ggsurvplot(survfit(recsurv.2~perc.BM_Clust.4.4, data = stage2CRC),size = 0.4, censor.shape="|", censor.size = 2,
                data = stage2CRC, ylim = c(0, 1),font.legend = list(size = 8),legend.labs = c("High", "Low"),
                palette = c("#fb4141", "#1f1fd8"), # custom color palette 
                risk.table = TRUE,
                pval = F,
                conf.int = F,xlab="Time (years)",legend = c(0.5,0.42),
                ggtheme = theme(axis.title = element_text(size = 12, family = "Arial"),
                                axis.text = element_text(size = 6),
                                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
                                panel.border = element_blank(), axis.line = element_line(), legend.key = element_rect(fill = NA), 
                                legend.text = element_text(size = 6.1),legend.key.width = unit(0.5, "mm"), legend.margin=margin(2,0,0,0),
                                legend.box.margin=margin(-7,-10,-10,-10), 
                                plot.margin = margin(0, 0, 0, 0.1, "cm")),
                risk.table.y.text.col = T,
                risk.table.y.text = F,
                risk.table.height = 0.1,
                legend.title="",  #remove Strata
                fontsize = 2)+   # adjust size for risk table
  guides(color = guide_legend(nrow = 1))


k$table <- k$table+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(),legend.key.width = unit(0.3, "mm"),
                         axis.line = element_blank(), axis.ticks = element_blank())
k
############### Adding values to the plot manually
grid.text(label = "T1_signature", gp = gpar(fontsize = 8), x = unit(0.325, "npc"), y = unit(0.515, "npc"))  
grid.text(label = "HR = 3.36(CI 1.5-7.54)", gp = gpar(fontsize = 8), x = unit(0.419, "npc"), y = unit(0.425, "npc"))  
grid.text(label = "P = 0.003", gp = gpar(fontsize = 8, fontface = "italic"), x = unit(0.3, "npc"), y = unit(0.35, "npc"))  

dev.off()


##### calculate HR using "survivalAnalysis" package
#### convert years to days (because the survivalAnalysis package doesn't work with years)
time2 <- (stage2CRC$RFS)*30.4
stage2CRC$Time <- round(time2)

stage2CRC %>%
  mutate(perc.BM_Clust.4.4=recode_factor(perc.BM_Clust.4.4, `Low`="1", `High`="2")) %>% 
  analyse_survival(vars(Time, Rec.Status), by=perc.BM_Clust.4.4) -> resultCRC #store the result object for later use



####################################### Figure 6E ##########################

#### read DEGs of AKA vs WT mouse models
dif <- read.delim("mouse_models/shrunkenLogFC_AKA_WT.csv", sep = ",")
head(dif)
## make a ranked list of ENSEMBL_IDs based on Shrunken log2 foldchanges (LFC)
gene_list <- as.vector(dif$log2FoldChange)
names(gene_list) <- dif$ENSEMBL_ID
ensemble_list <- sort(gene_list, decreasing = T)
head(ensemble_list)

##### Fessler mouse signature with Ensemble IDs
fessler <- read.delim("mouse_models/FESSLER_2016_Ensembl_ID_TGFB_UP_2.gmt", sep = "\t", header = T)
fessler <- t(fessler)
## convert the fessler signature to 'list' file
fessler <- list(TGFB_UP = as.character(row.names(fessler)))
## apply fgsea to do pre-ranked GSEA
set.seed(127)
fgsea(fessler, stats = ensemble_list , 1000, minSize = 1, maxSize = Inf, nproc = 0, gseaParam = 1, BPPARAM = NULL)

# fessler
tiff("Fig.6E_fessler.tiff", units="mm", width=55, height=40, res=300)
plotEnrichment(fessler[["TGFB_UP"]], gene_list, gseaParam = 1, ticksSize = 0.1)+labs(title="TGFB_UP")+labs(y = "Enrichment scores")+
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", size = 12), panel.background = element_blank(), 
        axis.title.y = element_text(size = 12, family = "Arial"),axis.text.y = element_text(size = 8, family = "Arial"),
        panel.border = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
grid.text(label = "NES = 1.98 \npadj = 0.01", gp = gpar(fontsize = 8), x = unit(0.8, "npc"), y = unit(0.68, "npc"))
grid.text(label = "Fessler et al. 2016", gp = gpar(fontsize = 8, fontfamily = "Arial"), x = unit(0.6, "npc"), y = unit(0.82, "npc"))
grid.text(label = "AKACa", gp = gpar(fontsize = 8), x = unit(0.3, "npc"), y = unit(0.04, "npc"))
grid.text(label = "WT", gp = gpar(fontsize = 8), x = unit(0.85, "npc"), y = unit(0.04, "npc"))
dev.off()

####################################### Figure 6F ##########################

## make a ranked list of gene symbols based on Shrunken log2 foldchanges (LFC)
gene_list <- as.vector(dif$log2FoldChange)
names(gene_list) <- dif$Gene_Symbol
gene_list <- sort(gene_list, decreasing = T)
head(gene_list)

### Get Hallmark genesets from msigdb
H.mus <- msigdbr(species = "Mus musculus", category = "H")
H_list_mus <- H.mus %>% split(x = .$gene_symbol, f = .$gs_name)  ## 7526

## apply fgsea to do pre-ranked GSEA
set.seed(127)
fgsea(H_list_mus, stats = gene_list , 1000, minSize = 1, maxSize = Inf, nproc = 0, gseaParam = 1, BPPARAM = NULL)

# tgfb
tiff("Fig.6F_tgfb_hallmark.tiff", units="mm", width=55, height=40, res=300)
plotEnrichment(H_list_mus[["HALLMARK_TGF_BETA_SIGNALING"]], gene_list, gseaParam = 1, ticksSize = 0.1)+labs(title="HALLMARK_TGFB")+labs(y = "Enrichment scores")+
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", size = 12), panel.background = element_blank(), 
        axis.title.y = element_text(size = 12, family = "Arial"),axis.text.y = element_text(size = 8, family = "Arial"),
        panel.border = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
grid.text(label = "NES = 2.28 \npadj = 0.01", gp = gpar(fontsize = 8), x = unit(0.8, "npc"), y = unit(0.68, "npc"))
grid.text(label = "AKACa", gp = gpar(fontsize = 8), x = unit(0.3, "npc"), y = unit(0.04, "npc"))
grid.text(label = "WT", gp = gpar(fontsize = 8), x = unit(0.85, "npc"), y = unit(0.04, "npc"))
dev.off()



####################################### Figure 6G and Extended Data Figure 7d ##########################

## Read mouse CRIS template
load("mouse_models/Mouse_Ensembl_ID_CRIS_NTP_template.Rdata") #template_CRIS
## convert template_CRIS to 'list' file
cris <- list(CRIS_A = template_CRIS[template_CRIS$class=="CRIS-A",2], 
             CRIS_B = template_CRIS[template_CRIS$class=="CRIS-B",2], 
             CRIS_C = template_CRIS[template_CRIS$class=="CRIS-C",2], 
             CRIS_D = template_CRIS[template_CRIS$class=="CRIS-D",2], 
             CRIS_E = template_CRIS[template_CRIS$class=="CRIS-E",2])

## apply fgsea to do pre-ranked GSEA
set.seed(127)
fgsea(cris, stats = ensemble_list , 1000, minSize = 1, maxSize = Inf, nproc = 0, gseaParam = 1, BPPARAM = NULL)

## plot crisA
tiff("Extended Data Figure 7d_crisA.tiff", units="mm", width=55, height=40, res=300)
plotEnrichment(pathways[["CRIS_A"]], gene_list, gseaParam = 1, ticksSize = 0.1)+labs(title="CRIS-A signature")+labs(y = "Enrichment scores")+
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", size = 12), panel.background = element_blank(), 
        axis.title.y = element_text(size = 12, family = "Arial"),axis.text.y = element_text(size = 8, family = "Arial"),
        panel.border = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
grid.text(label = "NES = 1.43 \npadj = 0.01", gp = gpar(fontsize = 8), x = unit(0.8, "npc"), y = unit(0.68, "npc"))
grid.text(label = "AKACa", gp = gpar(fontsize = 8), x = unit(0.3, "npc"), y = unit(0.04, "npc"))
grid.text(label = "WT", gp = gpar(fontsize = 8), x = unit(0.85, "npc"), y = unit(0.04, "npc"))
dev.off()

## plot crisB
tiff("Fig.6G_crisB.tiff", units="mm", width=55, height=40, res=300)
plotEnrichment(pathways[["CRIS_B"]], gene_list, gseaParam = 1, ticksSize = 0.1)+labs(title="CRIS-B signature")+labs(y = "Enrichment scores")+
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", size = 12), panel.background = element_blank(), 
        axis.title.y = element_text(size = 12, family = "Arial"),axis.text.y = element_text(size = 8, family = "Arial"),
        panel.border = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
grid.text(label = "NES = 1.87 \npadj = 0.01", gp = gpar(fontsize = 8), x = unit(0.8, "npc"), y = unit(0.68, "npc"))
grid.text(label = "AKACa", gp = gpar(fontsize = 8), x = unit(0.3, "npc"), y = unit(0.04, "npc"))
grid.text(label = "WT", gp = gpar(fontsize = 8), x = unit(0.85, "npc"), y = unit(0.04, "npc"))
dev.off()


########## CMS
## Read mouse CMS template
load("mouse_models/Mouse_Ensembl_ID_CMS_NTP_template.Rdata") # template_CMS

## convert template_CMS to 'list' file
cms <- list(CMS1 = template_CMS[template_CMS$class=="CMS1",2], 
            CMS2 = template_CMS[template_CMS$class=="CMS2",2], 
            CMS3 = template_CMS[template_CMS$class=="CMS3",2], 
            CMS4 = template_CMS[template_CMS$class=="CMS4",2])

## apply fgsea to do pre-ranked GSEA
set.seed(124)
fgsea(cms, stats = ensemble_list , 1000, minSize = 1, maxSize = Inf, nproc = 0, gseaParam = 1, BPPARAM = NULL)


## plot CMS1
tiff("Extended Data Figure 7d_CMS1.tiff", units="mm", width=55, height=40, res=300)
plotEnrichment(cms[["CMS1"]], gene_list, gseaParam = 1, ticksSize = 0.1)+labs(title="CMS1 signature")+labs(y = "Enrichment scores")+
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", size = 12), panel.background = element_blank(), 
        axis.title.y = element_text(size = 12, family = "Arial"),axis.text.y = element_text(size = 8, family = "Arial"),
        panel.border = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
grid.text(label = "NES = 1.68 \npadj = 0.01", gp = gpar(fontsize = 8), x = unit(0.8, "npc"), y = unit(0.68, "npc"))
grid.text(label = "AKACa", gp = gpar(fontsize = 8), x = unit(0.3, "npc"), y = unit(0.04, "npc"))
grid.text(label = "WT", gp = gpar(fontsize = 8), x = unit(0.85, "npc"), y = unit(0.04, "npc"))
dev.off()


## plot CMS2
tiff("Extended Data Figure 7d_CMS2.tiff", units="mm", width=55, height=40, res=300)
plotEnrichment(cms[["CMS2"]], gene_list, gseaParam = 1, ticksSize = 0.1)+labs(title="CMS2 signature")+labs(y = "Enrichment scores")+
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", size = 12), panel.background = element_blank(), 
        axis.title.y = element_text(size = 12, family = "Arial"),axis.text.y = element_text(size = 8, family = "Arial"),
        panel.border = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
grid.text(label = "NES = -1.25 \npadj = 0.28", gp = gpar(fontsize = 8), x = unit(0.42, "npc"), y = unit(0.25, "npc"))
grid.text(label = "AKACa", gp = gpar(fontsize = 8), x = unit(0.3, "npc"), y = unit(0.04, "npc"))
grid.text(label = "WT", gp = gpar(fontsize = 8), x = unit(0.85, "npc"), y = unit(0.04, "npc"))
dev.off()

####################################### Figure 6H #############################

## read mouse normalized data
df <- read.csv("mouse_models/normalized_all_models.csv")
df  <- column_to_rownames(df , var = "Gene_Symbol")[-1]

#### read meta data
meta <- read.delim("mouse_models/meta_all_models.csv")
## reorder names
meta.3 <- meta[c(45:47, 1:6, 16:18, 38:44, 7:12, 19:25, 35:37, 32:34, 26:31),]
rownames(meta.3) <- NULL
meta.3 <- column_to_rownames(meta.3, var = "Sample_ID")

## rearrange colnames in expression file to be in the same order as rownames in meta file
df.2 <- df[, rownames(meta.3)]
## check the order
all(colnames(df.2)== rownames(meta.3))

## read CRIS-B gene signature 
cris <- read.delim("mouse_models/template_CRIS_symbol.txt")
crisB <- list(crisB=cris[cris$class=="CRIS-B",3])


### ssGSEA --- GSVA
set.seed(124)
ssgsea_crisB <- gsva(as.matrix(df.2), crisB, 
                     min.sz=10, max.sz=Inf, verbose = T, method = 'ssgsea')

t_ssgsea_crisB <- as.data.frame(t(ssgsea_crisB))
t_ssgsea_crisB <- rownames_to_column(t_ssgsea_crisB, var = "label")
t_ssgsea_crisB$Type <- factor(meta.3$Abbreviation, levels = unique(meta.3$Abbreviation))


#################################---------------------------------------------------------------------------------------

##### only keep AK, AKA, AKAlk5ca
t_ssgsea_crisB.2 <- t_ssgsea_crisB[which(t_ssgsea_crisB$Type %in% c("AK", "AKA", "AKAlk5KO")),]

#t_ssgsea_crisB.2$group <- c(rep("AK", 3), rep("AKA", 7), rep("AKAlk5caKO", 3))
t_ssgsea_crisB.2$Type <- factor(t_ssgsea_crisB.2$Type, levels = unique(t_ssgsea_crisB.2$Type))


################## plot to start from zero
## calculate median for AK
AKA_median = median(t_ssgsea_crisB.2[t_ssgsea_crisB.2$Type=="AKA",2])
## calculate relative score to force plot to start from zero
t_ssgsea_crisB.2[,4] <- (t_ssgsea_crisB.2$crisB)-(AKA_median)


tiff("Figure 6H.tiff", units="mm", width=67, height=80, res=300)
p_crisB <- ggplot(t_ssgsea_crisB.2, aes(x = Type, y =  V4))+
  geom_boxplot(fill="black")+geom_quasirandom(width=.1, aes(color=Type), size=3)+
  labs(y = "CRIS.B / ssGSEA score")+geom_hline(aes(yintercept=0), linetype="dotted")+
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 3))+
  theme_bw()+theme(
    legend.title = element_blank(),
    legend.position = c(0.7835,0.80),
    legend.text = element_text(margin = margin(r = 0.15, unit = "cm"), size=7.5),
    legend.key.size = unit(0.13, "cm"),
    legend.spacing.x = unit(0, 'cm'),
    #legend.position = "none",
    axis.title.x=element_blank(), 
    #plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.y = element_text(size = 14, family = "Helvetica"),
    #axis.text.x = element_text(size = 11, family = "Helvetica", angle = 45, margin = margin(t = 25, r = 20, b = 0, l = 0)),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
    axis.line = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank())+
  scale_color_manual(values = c('AK'='#5fa171', 'AKAlk5KO'='#a15da6', 'AKA' = '#bac1c2'))

#### Add p-value (more beautiful)
p_crisB+geom_signif(comparisons = list(c("AKA", "AK"), c("AKA", "AKAlk5KO")),
                    step_increase = 0.1, tip_length = 0, vjust=0.2,
                    map_signif_level=TRUE, test = "t.test", color = "black")

dev.off()


####################################### Figure 6I #############################

selected_genes <- df.2[c("Snai3", "Cbl", "Mycl", "Itgb6", "Vegfa", "Egfr", "Itgav", "Tgfbr1", "Mycn", "Fos", "Smurf1", 
                         "Cblb", "Skil", "Map3k13", "Smad7", "Pdgfb", "Mycbp2", "Ret"),]

################## scale
selected_genes_scaled <- t(scale(t(selected_genes)))

### meta for heatmap
Genotype <- factor(meta.3$Abbreviation, levels = unique(meta.3$Abbreviation))
Genotype <- as.data.frame(Genotype)
rownames(Genotype) <- rownames(meta.3)


## plot heatmap
ann_colors = list(Genotype = c(WT="#967518", A = "#7570B3", AK = "#5fa171", Alk5ca = "#47d1c6",  AA="#2802fa", AKA="#bac1c2", 
                               AKAS="#993404", AKAlk5KO="#a15da6", AKA_AZD6244="#ff4242", AKA_EGFRi="#1f1fd8"))

tiff("Fig.6I.tiff", units="mm", width=190, height=90, res=300)
my_palette <- colorRampPalette(c("#1f1fd8", "#fffcfc", "#fb4141"))(n = 1000)
pheatmap(selected_genes_scaled, angle_col = "45", color=my_palette, fontsize_row =12, border_color = NA, family = "Arial",
         cluster_cols = F, cluster_rows = F, legend = T,  gaps_col=c(3,9,12,19,25,32, 35,38,41,44), show_colnames=F,
         annotation_col = Genotype, annotation_colors=ann_colors)
dev.off()


####################################### Figure 7F #############################

##### only keep AK, AKA, AKAlk5ca from 't_ssgsea_crisB' file
t_ssgsea_crisB.3 <- t_ssgsea_crisB[which(t_ssgsea_crisB$Type %in% c("AKA_AZD6244", "AKA", "AKA_EGFRi")),]
t_ssgsea_crisB.3$Type <- factor(t_ssgsea_crisB.3$Type, levels = unique(t_ssgsea_crisB.3$Type))

################## plot to start from zero
## calculate median for AK
AKA_median = median(t_ssgsea_crisB.3[t_ssgsea_crisB.3$Type=="AKA",2])
## calculate relative score to force plot to start from zero
t_ssgsea_crisB.3[,4] <- (t_ssgsea_crisB.3$crisB)-(AKA_median)


tiff("Figure 7F.tiff", units="mm", width=65, height=80, res=300)

p_crisB <- ggplot(t_ssgsea_crisB.3, aes(x = Type, y =  V4))+
  geom_boxplot(fill="black")+geom_quasirandom(width=.1, aes(color=Type), size=3)+
  labs(y = "CRIS.B / ssGSEA score")+geom_hline(aes(yintercept=0), linetype="dotted")+
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 3))+
  theme_bw()+theme(
    legend.title = element_blank(),
    legend.position = c(0.71,0.70),
    legend.text = element_text(margin = margin(r = 0.15, unit = "cm"), size=7.5),
    legend.key.size = unit(0.13, "cm"),
    legend.spacing.x = unit(0, 'cm'),
    axis.title.x=element_blank(), 
    axis.title.y = element_text(size = 14, family = "Helvetica"),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
    axis.line = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank())+
  scale_color_manual(values = c('AKA_EGFRi'='#5fa171', 'AKA_AZD6244'='#a15da6', 'AKA' = '#bac1c2'))
#geom_text_repel()+

#### Add p-value (more beautiful)
p_crisB+geom_signif(comparisons = list(c("AKA", "AKA_EGFRi"), c("AKA", "AKA_AZD6244")),
                    step_increase = 0.1, tip_length = 0, vjust=0.2,
                    map_signif_level=TRUE, test = "t.test", color = "black")

dev.off()


#####################################################################################################################################################
##################################################### Extended Data Figure 7 ########################################################################
#####################################################################################################################################################

################################### Extended Data Figure 7A #########################

stage2.cell <- read.csv("stageII_CRC_Sanghee/CellData_Standardized_stage2.csv")

### select only columns with "pS6", "4EBP1", "TGFB1", "pmTOR", "Clus" value
df <- stage2.cell[,c(6,10,14,17,21)]
## renames colnames
colnames(df) <- c("pS6", "4EBP1", "TGFB1", "pmTOR", "Clus")

###################### calculate median for each marker per cluster
df_cl1.median <- apply(df[df$Clus == "1",-5], 2, median)
df_cl2.median <- apply(df[df$Clus == "2",-5], 2, median)
df_cl3.median <- apply(df[df$Clus == "3",-5], 2, median)
df_cl4.median <- apply(df[df$Clus == "4",-5], 2, median)

df.median <- data.frame(Cluster1 = df_cl1.median, Cluster2 = df_cl2.median, Cluster3 = df_cl3.median, Cluster4 = df_cl4.median) 


tiff("Extended Data Figure 7A.tiff", units="mm", width=80, height=85, res=300)
h.col <- colorRampPalette(c('blue', 'white', 'red'))(100)
heatmap.2(as.matrix(df.median), 
          col = h.col, 
          scale="none",
          margins=c(8,0), # ("margin.Y", "margin.X")
          trace='none',
          srtCol=45,
          symkey=FALSE, 
          symbreaks=FALSE, 
          dendrogram='none',
          density.info='none', 
          denscol="black",
          key.title = NA,
          key.xlab="Median_value",
          keysize=1, 
          #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
          key.par=list(mar=c(4.5,2,2,2)),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2.5, 5), lwid=c(3, 10, 7))

dev.off()


################################### Extended Data Figure 7B #########################

stage2CRC <- read.delim("stageII_CRC_Sanghee/patient.4clusters_2.csv", sep = ",", stringsAsFactors = F)
dim(stage2CRC)
head(stage2CRC)
## plot heatmap
stage2.heat <- stage2CRC
rownames(stage2.heat) <- stage2CRC[,1]
stage2.heat <- stage2.heat[,-c(1:6,11:13)]
colnames(stage2.heat) <- paste0("Cluster", 1:4)
head(stage2.heat)

## plot
tiff("Extended Data Figure 7B.tiff", units="mm", width=90, height=95, res=300)
#par(oma=c(1,5,1,1))

h.col <- colorRampPalette(c('blue', 'white', 'red'))(100)
heatmap.2(as.matrix(stage2.heat), 
          col = h.col,
          #ylab = NULL,
          labRow = FALSE,
          ylab = "Patients",
          scale="none",
          margins=c(8,2), # ("margin.Y", "margin.X")
          trace='none',
          srtCol=45,
          symkey=FALSE, 
          symbreaks=FALSE, 
          dendrogram='none',
          density.info='none', 
          denscol="black",
          key.title = NA,
          key.xlab="Cluster proportion",
          keysize=1, 
          #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
          key.par=list(mar=c(4.5,2,2,3.5)),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2.5, 5), lwid=c(3, 10, 7))
dev.off()


################################### Extended Data Figure 7C #########################

stage2.cell_1 <- stage2.cell[,-c(1:5)]
head(stage2.cell_1)

## select only values for cluster 4
stage2.cell_cl4 <- stage2.cell_1[stage2.cell_1$Clus == "4",]
head(stage2.cell_cl4)
## renames colnames
name_df <- colnames(stage2.cell_cl4)
colnames(stage2.cell_cl4) <- substr(colnames(stage2.cell_cl4), 11, nchar(name_df))
stage2.cell_cl4 <- stage2.cell_cl4[,-16]

## calculate mean per column
m <- colMeans(stage2.cell_cl4)
##### calculate quantile instead of min & max
Quan <- apply(stage2.cell_cl4, 2, quantile)
lower <- Quan["25%",]
upper <- Quan["75%",]

df_Quan <- data.frame(label=names(stage2.cell_cl4), m, lower, upper)
df_Quan$type <- as.factor(c(1,2,2,2,1,2,2,2,1,2,2,1,2,2,2)) ## to add color

## plot
tiff("Extended Data Figure 7C.tiff", units="mm", width=80, height=90, res=300)
g <- ggplot(data=df_Quan, aes(x=label, y=m, ymin=lower, ymax=upper, color = type))+
  geom_pointrange(size = 0.2)+
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Log2 expression") +
  theme_bw()+theme(legend.position = "none",axis.text.x = element_text(size = 6.2),plot.margin = margin(0.1, 0.3, 1, -0.1, "cm"),
                   axis.title.x = element_text(size = 12),axis.text.y = element_text(size = 12),
                   panel.grid.major = element_blank(),panel.grid.minor = element_blank())+annotate("text", x=1, y=0.9, label="cluster 4", color = "red", size = 5)
g+scale_color_manual(values=c( "1"="#ff4242", "2"="#1f1fd8"))
dev.off()



###############################################################################################
####################################   Session Info   #########################################
###############################################################################################

sessionInfo()

##########-------------------------------------------------------------------------------------
R version 3.6.3 (2020-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.3 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.2.20.so

locale:
  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8   
[6] LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
  [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] gplots_3.1.1           RColorBrewer_1.1-2     reshape_0.8.8          biomaRt_2.40.5         gridExtra_2.3          DEGreport_1.20.0      
[7] pheatmap_1.0.12        dplyr_1.0.5            plyr_1.8.6             tibble_3.1.1           msigdbr_7.2.1          GSVA_1.32.0           
[13] fgsea_1.10.1           Rcpp_1.0.6             magrittr_2.0.1         survivalAnalysis_0.1.3 survminer_0.4.9        survival_3.2-10       
[19] ggsignif_0.6.1         ggpubr_0.4.0           MCPcounter_1.1.0       curl_4.3               ggbeeswarm_0.6.0       ggplot2_3.3.3         

loaded via a namespace (and not attached):
  [1] readxl_1.3.1                backports_1.2.1             circlize_0.4.12             Hmisc_4.5-0                 fastmatch_1.1-0            
[6] ConsensusClusterPlus_1.48.0 GSEABase_1.46.0             splines_3.6.3               BiocParallel_1.18.1         GenomeInfoDb_1.20.0        
[11] digest_0.6.27               htmltools_0.5.1.1           fansi_0.4.2                 checkmate_2.0.0             memoise_2.0.0              
[16] cluster_2.1.0               openxlsx_4.2.3              limma_3.40.6                ComplexHeatmap_2.0.0        annotate_1.62.0            
[21] Nozzle.R1_1.1-1             matrixStats_0.58.0          prettyunits_1.1.1           jpeg_0.1-8.1                colorspace_2.0-0           
[26] blob_1.2.1                  ggrepel_0.9.1               haven_2.4.0                 xfun_0.22                   crayon_1.4.1               
[31] RCurl_1.98-1.3              graph_1.62.0                genefilter_1.66.0           zoo_1.8-9                   glue_1.4.2                 
[36] gtable_0.3.0                zlibbioc_1.30.0             XVector_0.24.0              GetoptLong_1.0.5            DelayedArray_0.10.0        
[41] car_3.0-10                  shape_1.4.5                 BiocGenerics_0.30.0         abind_1.4-5                 scales_1.1.1               
[46] DBI_1.1.1                   edgeR_3.26.8                rstatix_0.7.0               progress_1.2.2              xtable_1.8-4               
[51] lasso2_1.2-20               htmlTable_2.1.0             tmvnsim_1.0-2               clue_0.3-59                 foreign_0.8-72             
[56] bit_4.0.4                   km.ci_0.5-2                 Formula_1.2-4               stats4_3.6.3                httr_1.4.2                 
[61] htmlwidgets_1.5.3           ellipsis_0.3.1              farver_2.1.0                pkgconfig_2.0.3             XML_3.99-0.3               
[66] nnet_7.3-12                 locfit_1.5-9.4              utf8_1.2.1                  labeling_0.4.2              tidyselect_1.1.0           
[71] rlang_0.4.10                later_1.1.0.1               AnnotationDbi_1.46.1        munsell_0.5.0               cellranger_1.1.0           
[76] tools_3.6.3                 cachem_1.0.4                generics_0.1.0              RSQLite_2.2.7               broom_0.7.6                
[81] ggdendro_0.1.22             stringr_1.4.0               fastmap_1.1.0               yaml_2.2.1                  knitr_1.32                 
[86] bit64_4.0.5                 zip_2.1.1                   caTools_1.18.2              survMisc_0.5.5              purrr_0.3.4                
[91] nlme_3.1-141                mime_0.10                   rstudioapi_0.13             compiler_3.6.3              shinythemes_1.2.0          
[96] beeswarm_0.3.1              png_0.1-7                   geneplotter_1.62.0          stringi_1.5.3               forcats_0.5.1              
[101] lattice_0.20-38             Matrix_1.2-17               psych_2.1.3                 KMsurv_0.1-5                vctrs_0.3.7                
[106] pillar_1.6.0                lifecycle_1.0.0             tidytidbits_0.2.3           GlobalOptions_0.1.2         data.table_1.14.0          
[111] cowplot_1.1.1               bitops_1.0-6                httpuv_1.5.5                GenomicRanges_1.36.1        R6_2.5.0                   
[116] latticeExtra_0.6-29         promises_1.2.0.1            KernSmooth_2.23-15          rio_0.5.26                  vipor_0.4.5                
[121] IRanges_2.18.3              gtools_3.8.2                MASS_7.3-51.4               assertthat_0.2.1            SummarizedExperiment_1.14.1
[126] DESeq2_1.24.0               rjson_0.2.20                withr_2.4.2                 mnormt_2.0.2                S4Vectors_0.22.1           
[131] GenomeInfoDbData_1.2.1      parallel_3.6.3              hms_1.0.0                   rpart_4.1-15                tidyr_1.1.3                
[136] carData_3.0-4               logging_0.10-108            Biobase_2.44.0              shiny_1.6.0                 base64enc_0.1-3

