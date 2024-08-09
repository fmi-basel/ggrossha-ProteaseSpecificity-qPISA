# This file contains the R code used to perform various analyses and figures presented in the manuscript:

# "Massively parallel quantification of substrate turnover defines protease subsite cooperativity
# Rajani Kanth Gudipati, Dimos Gaidatzis, Jan Seebacher, Sandra Muehlhaeusser, Georg Kempf, 
# Simone Cavadini, Daniel Hess, Charlotte Soneson and Helge Großhans"

# This R script produces figures and tables related to the analysis of dpp-4. The analyis of the dpf-3 
# experiment can be performed by replacing the dpp-4 peptide log2 fold changes by the dpf-3 peptide
# log2 fold changes (see line "replace by ")

library(SingleCellExperiment)
library(ggplot2)
library(ggseqlogo)
library(gplots)

# download this R data file from https://www.ebi.ac.uk/pride/, under the accession PXD042089
# it contains all the proteomics data required to perform the analysis below
inFile <- "PD_TMT_L2_210730_210804_JS_2083_TMT16plex_Fr_semiTrypsin_2MC_uniquePep_SPS65_NAc_PeptideGroups_Abundance_MinProb_center.median_limma_einprot0.7.0_sce_publi.rds"

sce <- readRDS(inFile)

# remove peptides with ambigous n terminus flanks
sce <- sce[nchar(do.call(rbind,strsplit(rowData(sce)[,"Annotated.Sequence"],"\\."))[,1]) == 3,]

# only use the first mapping position if peptide maps to multiple proteins
rowData(sce)[,"Positions.in.Master.Proteins"] <- sapply(strsplit(rowData(sce)[,"Positions.in.Master.Proteins"],";"),function(x){x[1]})

# remove peptides with low detection levels
sce <- sce[rowMeans(assay(sce, "log2_Abundance_norm")) > -4.7,]

# remove duplicate sequences according to max intensity
pep_to_ind_L <- split(1:nrow(sce),rowData(sce)[,"Sequence"])
pep_avgIntensity <- rowMeans(assay(sce, "log2_Abundance_norm"))

sce <- sce[sapply(pep_to_ind_L,function(x){x[which.max(pep_avgIntensity[x])]}),]

# grab the peptide intensities
TL <- assay(sce, "log2_Abundance_norm")
rownames(TL) <- rowData(sce)[,"Sequence"]

# extract peptide meta data
aln1 <- do.call(rbind,strsplit(rowData(sce)$Positions.in.Master.Proteins," "))
aln2 <- cbind(rowData(sce)[,"Sequence"],aln1,do.call(rbind,strsplit(gsub("]","",gsub("\\[","",aln1[,2])),"-")))

# define some color schemes
cbbPal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colSchemeAALet <- c("G","A","V","L","I","F","M","P","W","S","T","Y","C","N","Q","D","E","H","K","R")
colSchemeAAGroup <- c(rep("Nonpolar",9),rep("Polar",6),rep("Acidic",2),rep("Basic",3))
colSchemeAACol <- c(rep("black",9),rep(cbbPal[4],6),rep(cbbPal[7],2),rep(cbbPal[3],3))
colSchemeAA = make_col_scheme(chars=colSchemeAALet, groups=colSchemeAAGroup, cols=colSchemeAACol)
bwr.colors <- colorRampPalette(c("blue", "white", "red"))

# extract flanking sequence to determine if a peptide is a trypitc+ (or here called original as is was there before the digestion)
seqAnn3 <- do.call(rbind,strsplit(rowData(sce)[,"Annotated.Sequence"],"\\."))
seqAnn3[,1] <- substr(seqAnn3[,1],2,nchar(seqAnn3[,1])-1)
seqAnn3[,3] <- substr(seqAnn3[,3],2,nchar(seqAnn3[,3])-1)
seqAnn3_dAA <- cbind(seqAnn3[,1],substr(seqAnn3[,2],nchar(seqAnn3[,2]),nchar(seqAnn3[,2])))
# determine original peptides (tryptic or at the termini)
selTryptic <- (seqAnn3_dAA[,1] %in% c("R","K")) & (seqAnn3_dAA[,2] %in% c("R","K"))
selTermini <- (seqAnn3[,1]=="-" & (seqAnn3_dAA[,2] %in% c("R","K"))) | (seqAnn3[,3]=="-" & (seqAnn3_dAA[,1] %in% c("R","K")))

selOriginal <- selTryptic | selTermini
selNotOriginal <- !selOriginal


#pdf("plots/pairsCorrs.pdf")
CORM <- cor(TL);diag(CORM) <- max(CORM[lower.tri(CORM)])
heatmap.2(CORM, dendrogram="none",Colv = NA,Rowv="none",col=colorRampPalette(c("black","white"))(256),trace="none",density.info="none",margins=c(14,14))
#dev.off()


# averge the intensities of the 3 replicates
TLA <- 1/3*(TL[,seq(1,ncol(TL),by=3)]+TL[,seq(2,ncol(TL),by=3)]+TL[,seq(3,ncol(TL),by=3)])
colnames(TLA) <- gsub("_S[0-9]+","",colnames(TLA))

plLims <- range(TLA)

#pdf("plots/dpf3_dpp4_noenzyme_vs_wt.pdf",width=13,height=7,pointsize=16)
par(mfrow=c(1,2))
plot(TLA[,c(3,1)],pch='.',cex=2,col=selNotOriginal+1,xlim=plLims,ylim=plLims,xaxt="n",yaxt="n")
axis(1, at = seq(-8,8,by=2))
axis(2, at = seq(-8,8,by=2))
plot(TLA[,c(4,2)],pch='.',cex=2,col=selNotOriginal+1,xlim=plLims,ylim=plLims,xaxt="n",yaxt="n")
axis(1, at = seq(-8,8,by=2))
axis(2, at = seq(-8,8,by=2))
legend(x="topleft",legend=c("Tryptic+","Non-Tryptic+"),col=c("black","red"),bty="n", pch = "*",)
#dev.off()

# extract the first new peptide positions
MTseqsStart <- do.call(rbind,lapply(strsplit(rowData(sce)[,"Sequence"],""),function(x){x[1:6]}))
colnames(MTseqsStart) <- paste0("P",1:ncol(MTseqsStart))

# subset the tryptic+ seqences
MTseqsStartOrig <- MTseqsStart[selOriginal,]

# subset the tryptic+ peptides
TLAOrig <- TLA[selOriginal,]

# calculate the log2 peptide fold change and threshold negative numbers to zero
# this is the quantity that the linear models try to predict
# replace by pmin(TLAOrig[,1]-TLAOrig[,3],0) to perform all the analyses for dpf3 instead of dpp4
TLAOrigD <- pmin(TLAOrig[,2]-TLAOrig[,4],0)


# determine all possible amino acids
allAA <- sort(unique(as.vector(MTseqsStartOrig)))


# assess the performance of various linear models. Note that the naming scheme of the variables here 
# does not use proper nomenclature due to the requirement of special characters (P2,P1,P1',P2',P3',P4'). 
# These are renamed further down
fm_r_P1 <- summary(lm(TLAOrigD ~ P1,data=data.frame(MTseqsStartOrig)))$adj.r.squared
fm_r_P2 <- summary(lm(TLAOrigD ~ P2,data=data.frame(MTseqsStartOrig)))$adj.r.squared
fm_r_P3 <- summary(lm(TLAOrigD ~ P3,data=data.frame(MTseqsStartOrig)))$adj.r.squared
fm_r_P4 <- summary(lm(TLAOrigD ~ P4,data=data.frame(MTseqsStartOrig)))$adj.r.squared
fm_r_P5 <- summary(lm(TLAOrigD ~ P5,data=data.frame(MTseqsStartOrig)))$adj.r.squared
fm_r_P6 <- summary(lm(TLAOrigD ~ P6,data=data.frame(MTseqsStartOrig)))$adj.r.squared
fm_r_P123 <- summary(lm(TLAOrigD ~ P1 + P2 + P3,data=data.frame(MTseqsStartOrig)))$adj.r.squared
fm_r_P123_I12 <- summary(lm(TLAOrigD ~ P1 + P2 + P3 + P1:P2,data=data.frame(MTseqsStartOrig)))$adj.r.squared
fm_r_P123_I23 <- summary(lm(TLAOrigD ~ P1 + P2 + P3 + P2:P3,data=data.frame(MTseqsStartOrig)))$adj.r.squared
fm_r_P123_I13 <- summary(lm(TLAOrigD ~ P1 + P2 + P3 + P1:P3,data=data.frame(MTseqsStartOrig)))$adj.r.squared
fm_r_P123_I11I23 <- summary(lm(TLAOrigD ~ P1 + P2 + P3 + P1:P2 + P2:P3,data=data.frame(MTseqsStartOrig)))$adj.r.squared

fm_P123456 <- lm(TLAOrigD ~ P1 + P2 + P3 + P4 + P5 + P6,data=data.frame(MTseqsStartOrig))
fm_r_P123456 <- summary(fm_P123456)$adj.r.squared
fm_P123456_CM <- rbind(A=0,matrix(fm_P123456$coefficients[-1],ncol=6))
rownames(fm_P123456_CM)[-1] <- substr(names(fm_P123456$coefficients[-1]),3,3)[1:19]
fm_P123456_CMN <- t(t(fm_P123456_CM)-apply(fm_P123456_CM,2,median)) # median normalize


#pdf("plots/dpp4_aa_model_6positions_seqlogo.pdf",width=8,height=6)
l6p <- ggseqlogo(-fm_P123456_CMN, method='custom',col_scheme=colSchemeAA) + ylab('Contribution') + theme_bw()
l6p <- l6p + theme(text=element_text(size=20),panel.spacing = unit(0.3,"cm"), panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank()) 
l6p + scale_x_continuous(breaks=c(1,2,3,4,5,6),labels= c("2","1","1'","2'","3'","4'"),position="bottom")
#dev.off()

# collect all the R^2 values from the models to show in a figure and use proper nomenclature
all_rs <- c(fm_r_P1,fm_r_P2,fm_r_P3,fm_r_P4,fm_r_P5,fm_r_P6,fm_r_P123,fm_r_P123456,fm_r_P123_I13,fm_r_P123_I12,fm_r_P123_I23,fm_r_P123_I11I23)
names(all_rs) <- c("P2","P1","P1'","P2'","P3'","P4'","P2 + P1 + P1'","P2 + P1 + P1' + P2' + P3' + P4'",
                "P2 + P1 + P1' + P2:P1'","P2 + P1 + P1' + P2:P1","P2 + P1 + P1' + P1:P1'","P2 + P1 + P1' + P2:P1 + P1:P1'")

#pdf("plots/dpp4_aa_models_performance.pdf",width=10,height=5.5,pointsize=16)
par(mar=c(5, 14, 1.5, 1.5))
bpCol <- c("darkgray","red")[rev(((1:length(all_rs)) %in% grep(":",names(all_rs)))+1)]
bpPos <- barplot(100*rev(all_rs),las=2,horiz=TRUE,xlab="Explained variance (adjusted, %)",xlim=c(0,max(100*all_rs)+15),col=bpCol)
text(100*rev(all_rs),bpPos-0.1,paste0(round(100*rev(all_rs),1),"%"),pos=4)
#dev.off()

# put together the final model
fullModel_mm <- model.matrix(TLAOrigD ~ P1 + P2 + P3 + P1:P2 + P2:P3,data=data.frame(MTseqsStartOrig)) # define the model matrix
fullModel_mms <- fullModel_mm[,colSums(fullModel_mm ) >= 10] # remove variables that are supported with less than 10 data points
fullModel <- lm(TLAOrigD ~ 0 + fullModel_mms) # fit the model

# extract the model coefficients
fullModel_coeffs <- fullModel$coefficients
names(fullModel_coeffs) <- gsub("fullModel_mms","",names(fullModel_coeffs))

fullModel_coeffs_P1 <- c(A=0,fullModel_coeffs[grep("^P1.$",names(fullModel_coeffs))])
names(fullModel_coeffs_P1) <- gsub("P1","",names(fullModel_coeffs_P1))

fullModel_coeffs_P2 <- c(A=0,fullModel_coeffs[grep("^P2.$",names(fullModel_coeffs))])
names(fullModel_coeffs_P2) <- gsub("P2","",names(fullModel_coeffs_P2))
fullModel_coeffs_P2P <- fullModel_coeffs_P2 + fullModel_coeffs[1]

fullModel_coeffs_P3 <- c(A=0,fullModel_coeffs[grep("^P3.$",names(fullModel_coeffs))])
names(fullModel_coeffs_P3) <- gsub("P3","",names(fullModel_coeffs_P3))

fullModel_coeffs_P1P2 <- fullModel_coeffs[grep("P1.:P2",names(fullModel_coeffs))]
names(fullModel_coeffs_P1P2) <- gsub("P2","",gsub("P1","",names(fullModel_coeffs_P1P2)))

fullModel_coeffs_P2P3 <- fullModel_coeffs[grep("P2.:P3",names(fullModel_coeffs))]
names(fullModel_coeffs_P2P3) <- gsub("P3","",gsub("P2","",names(fullModel_coeffs_P2P3)))

# reorganzie the coefficients such that the model can be visualized as a collection of sequence logos
fullModel_coeffs_P1P2_M <- matrix(NA,nrow=length(allAA),ncol=length(allAA),dimnames=list(allAA,allAA))
fullModel_coeffs_P1P2_M[do.call(rbind,strsplit(names(fullModel_coeffs_P1P2),":"))] <- fullModel_coeffs_P1P2
fullModel_coeffs_P1P2_M[1,] <- 0
fullModel_coeffs_P1P2_M[,1] <- 0

fullModel_coeffs_P2P3_M <- matrix(NA,nrow=length(allAA),ncol=length(allAA),dimnames=list(allAA,allAA))
fullModel_coeffs_P2P3_M[do.call(rbind,strsplit(names(fullModel_coeffs_P2P3),":"))] <- fullModel_coeffs_P2P3
fullModel_coeffs_P2P3_M[1,] <- 0
fullModel_coeffs_P2P3_M[,1] <- 0

fullModel_coeffs_P1P2_MP <- fullModel_coeffs_P1P2_M + fullModel_coeffs_P1
fullModel_coeffs_P2P3_MP <- t(t(fullModel_coeffs_P2P3_M) + fullModel_coeffs_P3) 

fullModel_coeffs_P1P2_MP_cm <- colMeans(fullModel_coeffs_P1P2_MP,na.rm=TRUE)
fullModel_coeffs_P2P3_MP_cm <- rowMeans(fullModel_coeffs_P2P3_MP,na.rm=TRUE)

fullModel_coeffs_P1P2_MPN <- t(t(fullModel_coeffs_P1P2_MP) - fullModel_coeffs_P1P2_MP_cm)
fullModel_coeffs_P2P3_MPN <- fullModel_coeffs_P2P3_MP - fullModel_coeffs_P2P3_MP_cm
fullModel_coeffs_P2PN <- fullModel_coeffs_P2P + fullModel_coeffs_P1P2_MP_cm + fullModel_coeffs_P2P3_MP_cm

ipL <- lapply(allAA,function(c){-cbind(P1=fullModel_coeffs_P1P2_MPN[,c],P2=fullModel_coeffs_P2PN*as.numeric(allAA %in% c),P3=fullModel_coeffs_P2P3_MPN[c,])})
names(ipL) <- allAA

ipM <- t(sapply(1:length(ipL),function(i){M=ipL[[i]];v1=M[,1];names(v1)=paste0("P1:",names(v1));v3=M[,3];names(v3)=paste0("P3:",names(v3));c("P2"=M[names(ipL)[i],2],v1,v3)}))
rownames(ipM) <- names(ipL)

# rename parameters (P1->P2 P2->P1 P3->P1')
ipMR <- ipM
colnames(ipMR) <- gsub("XPL","P2",gsub("P3","P1'",gsub("P2","P1",gsub("P1","XPL",colnames(ipM)))))

# save the model parameters
#write.table(round(ipMR,4),"dpp4_modelParams.txt",sep="\t",quote=FALSE)

# calculate the model prediction manually and check that these are the same as the fitting result from the liner model itself (but with opposite sign)
fullModel_pred <- apply(MTseqsStartOrig,1,function(x){ipL[[x[2]]][x[1],1]+ipL[[x[2]]][x[2],2]+ipL[[x[2]]][x[3],3]})
#plot(fullModel$fitted.values,fullModel_pred)

realVsPred <- as.matrix(data.frame("dpp4_experiment"=TLAOrigD,"dpp4_prediction"=fullModel_pred))

# save a table that compares the real log2 peptide fold change to the predicted one
#write.table(realVsPred,"dpp4_realVsPredicted.txt",sep="\t",quote=FALSE)

# extract the data to visualize in the enzyme specificity sequence logo
# sort the P1 amino acids by the maximum absolute score that each P1 position could achieve
ipL_S <- ipL[names(sort(sapply(ipL,function(x){sum(apply(abs(x),2,max,na.rm=TRUE),na.rm=TRUE)}),decreasing=TRUE))[1:10]]
names(ipL_S) <- paste("X ",names(ipL_S)," X")

#pdf("plots/dpp4_aa_model_2interactions_seqlogo.pdf",width=12,height=6)
gp1 <- ggseqlogo(ipL_S, method='custom',ncol=length(ipL_S),col_scheme=colSchemeAA) + ylab('Contribution') + theme_bw()
gp1 <- gp1 + theme(text=element_text(size=20),panel.spacing = unit(0.3,"cm"), panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank()) 
gp1 + scale_x_continuous(breaks=c(1,2,3),labels= c("2","1","1'"),position="top")
#dev.off()


# do a simple seq logo based on enrichment cutoff
sseqlogoLen <- 8
TLAOrigDS <- TLAOrigD[nchar(names(TLAOrigD)) >= sseqlogoLen]
TLAOrigDS_sel <- TLAOrigDS < -log2(10)

#pdf("plots/dpp4_10foldEnr_seqlogo.pdf",width=6,height=3.5)
TLAOrigDS_seqs <- substr(names(TLAOrigDS)[TLAOrigDS_sel],1,sseqlogoLen)
ggseqlogo(TLAOrigDS_seqs)
#dev.off()

# background amino acid distribution
#pdf("plots/AminoAcid_dpp4_bg_seqlogo.pdf",width=6,height=3.5)
ggseqlogo(substr(names(TLAOrigD)[nchar(names(TLAOrigD)) >= sseqlogoLen],1,8)) + ylim(0,3)
#dev.off()
