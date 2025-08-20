################################################################################
# Libraries
################################################################################

library(parallel)
library(gplots)
library(ggplot2)

################################################################################
# Functions
################################################################################

panel.cor.pearson <- function(x, y, digits = 2, prefix = "", cex.cor,  ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use="pairwise.complete.obs", method="pearson")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * abs(r))
}
panel.hist40 <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE, breaks=40)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
panel.smooth <- function(...){
  smoothScatter(..., nrpoints = 0, colramp = colorRampPalette(c("white", blues9[2:9])), add = TRUE)
}


################################################################################
# Load data
################################################################################

load("003-lfc_limma_AA_per_batch.Rdata")
he = lfc_he_aa
rm(lfc_he_aa)
load("004-mutation_effects.Rdata")
global_model = read.delim("005-thermodynamic_models/001-global_model/002-output/mochi_project/task_1/predictions/predicted_phenotypes_all.txt")
load("005-thermodynamic_models/001-global_model/001-dummy_conversion.Rdata")
load("005-thermodynamic_models/002-individual_models/003-Mochi_output.Rdata")
dl_pred = read.delim("000-data/TableS2_deep_learning_validation_set_predictions.txt")
load("006-synthetic_bZIPs_binding_scores.Rdata")
synth_bZ_pred = read.csv("000-data/binder_optimization.csv") # file is named results in the output of the deep learning github


all_hep = 1:5
all_pos = letters[1:7]
all_mut = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
all_hp = paste0(rep(all_hep,each=7),all_pos)

#colorPalette
colPal = colorRampPalette(c("cyan4","darkred"))
colPalDiff = colorRampPalette(c("cyan4","grey90", "darkred"))
colPalDiff2 = colorRampPalette(c("cyan4","grey99", "darkred"))

################################################################################
# Process data
################################################################################

# shorten wt_prey name
he$wt_prey = do.call("rbind",strsplit(he$full_id_prey,"_"))[,1]

# combine binding scores and mutation effects
he$dLFC = dlfc$logFC[match(paste(he$full_id_bait,he$full_id_prey),paste(dlfc$full_id_bait,dlfc$full_id_prey))]
he$dLFC.pval = dlfc$P.Value[match(paste(he$full_id_bait,he$full_id_prey),paste(dlfc$full_id_bait,dlfc$full_id_prey))]
rm(dlfc)

# add predictions from the global model
# first, convert the dummy sequence to that of the wild-type partner
global_model = cbind(conv[,1:2], global_model)
rm(conv)
he$global_model_pred = global_model$mean[match(paste(he$full_id_bait,he$full_id_prey), paste(global_model$full_id_bait,global_model$full_id_prey))]

# add predictions from individual models, as well as ddG_mut, ddG_bZ and dG (ddG_mut for the wild type baits correspond to dG_ref)
mo2 = do.call("rbind",lapply(mo,function(x){x$predictions}))

tmp = mo2[match(paste(he$wt_bait,he$id_bait,he$wt_prey), paste(mo2$wt_bait,mo2$id_bait,mo2$wt_prey)), c("mean","ddG_mut","ddG_bZ","dG")]
names(tmp)[1] = "ind_model_pred"
he = cbind(he,tmp)
rm(tmp)

# compute projected ddG
proj = do.call("rbind",lapply(mo,function(x){
  
  pr = x$predictions
  lin = x$linear_weights
  
  proj = log(lin[1]/(pr$fitness - lin[2]) - 1)
  
  pr$proj_dG = proj
  proj_dG_wt = proj[pr$id_bait == "0xx"]
  names(proj_dG_wt) = pr$wt_prey[pr$id_bait == "0xx"]
  pr$proj_dG_wt = proj_dG_wt[pr$wt_prey]
  
  pr$ddG_spe = pr$proj_dG - pr$dG
  pr$ddG_tot = pr$ddG_mut + pr$ddG_spe
  
  pr
  
}))

he = cbind(he, proj[match(paste(he$wt_bait,he$id_bait,he$wt_prey), paste(proj$wt_bait,proj$id_bait,proj$wt_prey)), c("proj_dG","proj_dG_wt","ddG_spe","ddG_tot")])

rm(mo2, proj, global_model)
gc()


# we keep a copy before filtering, because we need to filter first before computing mutation effect FDRs, and then add the FDRs back to the full table before writing it
he2 = he




# filter out low coverage variants and non-expected variants
he = he[he$i1 > 10 & he$i2 > 10 & he$i3 > 10 & he$o1 > 0 & he$o2 > 0 & he$o3 > 0 & he$expected_bait, ]

# filter out ATF7 as a bait since most variants are lost
he = he[he$wt_bait != "ATF7",]

# wild types
wt = he[he$id_bait == "0xx",]


# compute FDR for mutation effect. Do not consider WT since there's no mutation effect by definition
he$fdr = NA
he$fdr[he$id_bait != "0xx"] = p.adjust(he$dLFC.pval[he$id_bait!="0xx"],method="fdr") # compute FDR only for mutants since it doesn't make sense for WT interactions

he2$fdr = he$fdr[match(paste(he2$full_id_bait,he2$full_id_prey),paste(he$full_id_bait,he$full_id_prey))]
he2 = he2[,c(1:33,43,34:42)]
# save table S1
write.table(he2, file="000-data/TableS1_all_unfiltered_data.txt", sep="\t", quote=F, row.names=F)
rm(he2)


# vector of the 53 mutated bZ
all_baits = sort(unique(he$wt_bait))
# vector of all 54 WT partners
all_preys = sort(unique(he$wt_prey))

# all WT pairs
all_wt = as.data.frame(expand.grid(wt_bait = all_baits, wt_prey = all_preys))
all_wt$logFC = wt$logFC[match(paste(all_wt$wt_bait,all_wt$wt_prey), paste(wt$wt_bait,wt$wt_prey))]

he$wt_logFC = all_wt$logFC[match(paste(he$wt_bait,he$wt_prey), paste(all_wt$wt_bait,all_wt$wt_prey))]


all_pairs = as.data.frame(as.matrix(expand.grid(all_baits,all_preys,all_hep,all_pos,all_mut)))
names(all_pairs) = c("wt_bait","wt_prey","hep","pos","mut")

all_pairs$logFC = he$logFC[match(do.call("paste0",all_pairs), do.call("paste0",he[,c("wt_bait","wt_prey","id_bait")]))]
all_pairs$dLFC = he$dLFC[match(do.call("paste0",all_pairs[,1:5]), do.call("paste0",he[,c("wt_bait","wt_prey","id_bait")]))]
all_pairs = all_pairs[order(all_pairs$wt_bait,all_pairs$wt_prey,all_pairs$hep,all_pairs$pos,all_pairs$mut),]



tmp = as.vector(as.matrix(synth_bZ_pred[,6:59]))
synth_bZ_pred2 = data.frame(
  synth_bZ = rep(synth_bZ_pred$Name, 54),
  originate_from = rep(synth_bZ_pred$wt, 54),
  wt_partner = rep(names(synth_bZ_pred)[6:59],each=nrow(synth_bZ_pred)),
  pred_logFC = tmp
  )

lfc = merge(lfc, synth_bZ_pred2, by.x=c("bait","prey"), by.y=c("synth_bZ","wt_partner"), all=T)
lfc = lfc[,c(1,22,2:21,23)]

write.table(lfc,file="000-data/TableS3_synthetic_bZIPs.txt",quote=F,sep="\t",row.names=F)

################################################################################
# Fig S1 - pairs plot input and output
################################################################################

pdf("007-figures/FigS1_pairs_plot.pdf")
pairs(log10(he[,c("i1","i2","i3","o1","o2","o3")]), pch=".", upper.panel = panel.cor.pearson, diag.panel = panel.hist40, lower.panel = panel.smooth)
dev.off()

################################################################################
# Fig S2 - WT heatmap
################################################################################

m = t(matrix(all_wt$logFC, nrow=53)) #columns correspond to baits and rows to preys

rownames(m) = all_preys
colnames(m) = all_baits


co1 = cor(m,use="pairwise.complete.obs")
co2 = cor(t(m),use="pairwise.complete.obs")
hc1 = hclust(as.dist(1-co1))
hc2 = hclust(as.dist(1-co2))

m = m[hc2$order, hc1$order]

ra = range(m, na.rm=T)
ra[1] = floor(ra[1])
ra[2] = ceiling(ra[2])

pdf("007-figures/FigS2_wt_heatmap.pdf",width=8,height=8)
par(mar=c(4,4,4,4))
image(1:nrow(m), 1:ncol(m), m, col=colPal(ra[2]-ra[1]), zlim=ra, axes=F, xlab="", ylab="",main="")
axis(2,at = 1:ncol(m), colnames(m), cex.axis=0.5,las=2)
axis(1,at = 1:nrow(m), rownames(m), cex.axis=0.5,las=2)
box()

image(seq(ra[1],ra[2],1),1,matrix(seq(ra[1],ra[2],1)), col=colPal(ra[2]-ra[1]), zlim=ra, axes=F, xlab="binding score", ylab="",main="color key")
axis(1,at = seq(ra[1],ra[2],1))
box()

plot(as.dendrogram(hc1), main="baits")
plot(as.dendrogram(hc2), main="preys")

dev.off()





################################################################################
# Fig 1D - heatmap^2
################################################################################

mb = matrix(all_pairs$dLFC, ncol=53)
all_pairs = all_pairs[order(all_pairs$wt_prey),]
mp = matrix(all_pairs$dLFC, ncol=54)
all_pairs = all_pairs[order(all_pairs$wt_bait),]
colnames(mb) = all_baits
colnames(mp) = all_preys

hc1 = hclust(as.dist(1-cor(mb,use="pairwise.complete.obs")))
hc2 = hclust(as.dist(1-cor(mp,use="pairwise.complete.obs")))


pdf("007-figures/Fig1C_heatmap.pdf",width=54,height=53)
par(mfrow=c(53,54),mar=c(0,0,0,0),oma=c(3,3,3,3), xpd=TRUE)

ra = ceiling(max(abs(all_pairs$dLFC), na.rm=T))

for(i in hc1$labels[hc1$order]){
  for(j in hc2$labels[hc2$order]){
    cat("\r",i,"\t",j,"\t\t\t")
    m = t(matrix(all_pairs$dLFC[all_pairs$wt_bait == i & all_pairs$wt_prey == j], ncol=35))
    image(1:nrow(m), 1:ncol(m), m, col=colPalDiff(2*ra), zlim=c(-ra,ra), axes=F, xlab="", ylab="",main="")
    box()
    
  }
}
dev.off()


png("007-figures/Fig1D_heatmap.png",width=5400,height=5300)
par(mfrow=c(53,54),mar=c(0,0,0,0),oma=c(3,3,3,3), xpd=TRUE)

ra = ceiling(max(abs(all_pairs$dLFC), na.rm=T))

for(i in hc1$labels[hc1$order]){
  for(j in hc2$labels[hc2$order]){
    cat("\r",i,"\t",j,"\t\t\t")
    m = t(matrix(all_pairs$dLFC[all_pairs$wt_bait == i & all_pairs$wt_prey == j], ncol=35))
    image(1:nrow(m), 1:ncol(m), m, col=colPalDiff(2*ra), zlim=c(-ra,ra), axes=F, xlab="", ylab="",main="")
    box()
    
  }
}
dev.off()


pdf("007-figures/Fig1D_zoom_BACH1_MAFF.pdf", width=14, height=8)
i ="BACH1"
j = "MAFF"
m = t(matrix(all_pairs$dLFC[all_pairs$wt_bait == i & all_pairs$wt_prey == j], ncol=35))
image(1:nrow(m), 1:ncol(m), m, col=colPalDiff(2*ra), zlim=c(-ra,ra), axes=F, xlab="Position", ylab="Substitution",main="BACH1-MAFF")
box()
axis(1, at = 1:35, labels = paste0(rep(all_hep, each=7), rep(all_pos, 5)))
axis(2, at = 1:20, labels = all_mut, las=2)
dev.off()

pdf("007-figures/Fig1D_dendro_and_colorKey.pdf")
image(seq(-ra,ra,1),1,matrix(seq(-ra,ra,1)), col=colPalDiff(2*ra), zlim=c(-ra,ra), axes=F, xlab="mutation effect", ylab="",main="color key")
axis(1,at = seq(-ra,ra,1))
box()
plot(as.dendrogram(hc1), main="baits")
plot(as.dendrogram(hc2), main="preys")
dev.off()


################################################################################
# Fig S3 - distribution mutation effects
################################################################################

pdf("007-figures/FigS3_distri_mutation_effects.pdf", height=16, width=8)
par(mfrow=c(2,1))
xx=hist(he$dLFC,breaks=100, xlab="Mutation effect")
abline(v=c(-2,2))
xx$counts = log10(xx$counts/sum(xx$counts))+6
plot(xx,axes=F, xlab="Mutation effect", ylab="Frequency (log10)",ylim=c(0,6))
abline(v=c(-2,2))
axis(1)
axis(2, at=0:6, labels=-6:0)
dev.off()



# stats on significant mutation effects
nrow(he[!is.na(he$fdr),])
nrow(he[!is.na(he$fdr) & he$fdr < 0.05 & abs(he$dLFC) > 2,])
nrow(he[!is.na(he$fdr) & he$fdr < 0.05 & abs(he$dLFC) > 2,]) / nrow(he[!is.na(he$fdr),])
nrow(he[!is.na(he$fdr) & he$dLFC > 2 & he$fdr < 0.05,])
nrow(he[!is.na(he$fdr) & he$dLFC > 2 & he$fdr < 0.05,]) / nrow(he[!is.na(he$fdr) & he$fdr < 0.05 & abs(he$dLFC) > 2,]) 
nrow(he[!is.na(he$fdr) & he$dLFC > 2 & he$fdr < 0.05,]) / nrow(he[!is.na(he$fdr),])
nrow(he[!is.na(he$fdr) & he$dLFC < -2 & he$fdr < 0.05,])
nrow(he[!is.na(he$fdr) & he$dLFC < -2 & he$fdr < 0.05,]) / nrow(he[!is.na(he$fdr) & he$fdr < 0.05 & abs(he$dLFC) > 2,]) 
nrow(he[!is.na(he$fdr) & he$dLFC < -2 & he$fdr < 0.05,]) / nrow(he[!is.na(he$fdr),])



################################################################################
# Fig S4 - wild-type binding score vs No of significant interaction
################################################################################

he_int = split(he,list(he$wt_bait, he$wt_prey), drop=T)

# matrix with number of significant positive and negative mutation effects for each pair of bZIP
he_int2 = do.call("rbind",lapply(he_int,function(x){
  
  c(x$wt_logFC[1], nrow(x[!is.na(x$fdr) & x$fdr < 0.05 & x$dLFC > 2,]), nrow(x[!is.na(x$fdr) & x$fdr < 0.05 & x$dLFC < -2,]))
  
}))

pdf("007-figures/FigS4_wt_binding_score_vs_No_detrimental_mut.pdf")
plot(he_int2[,1],he_int2[,3],log="y", xlab="wild type binding score", ylab="No. of significantly detrimental mutations (log10)")
dev.off()


################################################################################
# Fig S5 - specificity rewiring
################################################################################

# add bait and prey names
tmp = do.call("rbind",strsplit(rownames(he_int2),"\\."))
he_int2 = data.frame(tmp,he_int2)

he$wt_int_decreased = he_int2[match(paste(he$wt_bait,he$wt_prey),paste(he_int2[,1],he_int2[,2])),5]

# how many mmutation create new interactions
mut_new_int = he[!is.na(he$fdr) & he$fdr < 0.05 & he$dLFC > 2 & he$wt_int_decreased == 0,]
nrow(mut_new_int)
nrow(mut_new_int[!duplicated(paste(mut_new_int$wt_bait, mut_new_int$wt_prey)),])

# For each variant, identify specificity rewiring mutations 
variant_list = split(he[!is.na(he$wt_logFC),], he$full_id_bait[!is.na(he$wt_logFC)])
spe_rewiring = unlist(lapply(variant_list, function(x){
  nrow(x[x$wt_int_decreased > 15 & x$wt_logFC > 1,]) > 0 & # the corresponding WT has at least one true WT interaction
  nrow(x[x$dLFC > 2 & x$fdr < 0.05 & x$logFC > 1 & x$wt_logFC < 1 & !is.na(x$fdr),]) > 0 & # at least one of the non-interacting partner is increased by the mutation
  nrow(x[x$wt_logFC > 1 & x$logFC > 1,]) == 0 # none of the partner interact with both the mutant and the WT
  
}))

# how many re there?
sum(spe_rewiring)


x = variant_list[["CREB5_N42L"]]

pdf("007-figures/FigS5_CREB5_3aL_spe_rewiring.pdf")
x$wt.se = he$se.logFC[match(paste(x$wt_bait, x$wt_prey), paste(he$wt_bait, he$wt_prey)[he$id_bait=="0xx"])]
plot(x$wt_logFC, x$logFC, xlab=paste(x$wt_bait[1], "wild type"), ylab = paste(x$wt_bait[1], x$id_bait[1]), xlim=c(min(x$wt_logFC - x$wt.se),max(x$wt_logFC + x$wt.se)), ylim=c(min(x$logFC-x$se.logFC),max(x$logFC+x$se.logFC)))
segments(x$wt_logFC - x$wt.se, x$logFC, x$wt_logFC + x$wt.se, x$logFC)
segments(x$wt_logFC, x$logFC - x$se.logFC, x$wt_logFC, x$logFC + x$se.logFC)
dev.off()

# which bZIPs are more involved
tmp = table(do.call("rbind",strsplit(names(variant_list[spe_rewiring]),"_"))[,1])
length(tmp)
mean(tmp)


################################################################################
# Fig S6 - Global model
################################################################################

pdf("007-figures/FigS6_global_model_performance.pdf")
par(mar=c(4,4,4,4))
ggplot(he, aes(x=global_model_pred, y=logFC) ) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis", trans="log", breaks=10^(0:6)) +
  theme_bw() +
  geom_abline(slope=1,intercept=0, linewidth = 0.2)
dev.off()

1 - sum((he$logFC-he$global_model_pred)^2,na.rm=T) / sum((he$logFC - mean(he$logFC))^2,na.rm=T)

################################################################################
# Fig 2C - distribution of individual Mochi models' R2
################################################################################

r = unlist(mclapply(all_baits, function(x){
  
  x = he[he$wt_bait == x,]
  1 - sum((x$logFC-x$ind_model_pred)^2,na.rm=T) / sum((x$logFC - mean(x$logFC))^2,na.rm=T)
  
},mc.cores=length(all_baits)))

pdf("007-figures/Fig2C_distri_R2.pdf")
hist(r, breaks=20)
dev.off()


fos = he[he$wt_bait == "FOS",]
1 - sum((fos$logFC-fos$ind_model_pred)^2,na.rm=T) / sum((fos$logFC - mean(fos$logFC))^2,na.rm=T)
fos = fos[fos$wt_prey %in% c("JUN","JUNB","JUND"),]
1 - sum((fos$logFC-fos$ind_model_pred)^2,na.rm=T) / sum((fos$logFC - mean(fos$logFC))^2,na.rm=T)




################################################################################
# Fig S7 - average binding score vs Mochi R2
################################################################################

av_binding=aggregate(he$logFC, list(he$wt_bait), mean, na.rm=T)

pdf("007-figures/FigS7_avLFC_vs_MochiR2.pdf")
plot(av_binding$x,r,xlab="Average binding score", ylab="MoCHI R2")
dev.off()

cor(av_binding$x,r, method="spearman")

################################################################################
# Fig 2D - pred vs measured binding scores, overall
################################################################################

r2 = 1 - sum((he$logFC-he$ind_model_pred)[!is.na(he$global_model_pred)]^2,na.rm=T) / sum((he$logFC - mean(he$logFC))[!is.na(he$global_model_pred)]^2,na.rm=T)

pdf("007-figures/Fig2D_pred_vs_obs.pdf", width=12.5, height=7)
par(mar=c(4,4,4,4))
ggplot(he[!is.na(he$global_model_pred),], aes(x=ind_model_pred, y=logFC) ) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis", trans="log", breaks=10^(0:6)) +
  theme_bw() +
  geom_abline(slope=1,intercept=0, linewidth = 0.2)
dev.off()


################################################################################
# Fig 2E; S8 - Correlation between ddG_mut
################################################################################

# unique mutations
umut = unique(he$std_id_bait[!is.na(he$ddG_mut)])

ddG_perfect_match = do.call("cbind",lapply(all_baits,function(x){
  
  x = he[he$wt_bait == x,]
  x = x[!duplicated(x$std_id_bait),]
  rownames(x) = x$std_id_bait
  x[umut,"ddG_mut"]
  
}))

colnames(ddG_perfect_match) = all_baits
rownames(ddG_perfect_match) = umut

corel = cor(ddG_perfect_match, use="pairwise.complete.obs")
hc_perfect_match = hclust(as.dist(1-corel))

ddG_perfect_match = ddG_perfect_match[,hc_perfect_match$order]
corel = corel[hc_perfect_match$order,hc_perfect_match$order]

range(corel[upper.tri(corel)])
mean(corel[upper.tri(corel)])
median(corel[upper.tri(corel)])

# now remove bZIPs for which the R2 of the Mochi model is < 0.4
bad_bZ = all_baits[r<0.4]

corel2 = corel[!(rownames(corel) %in% bad_bZ),!(colnames(corel) %in% bad_bZ)]
range(corel2[upper.tri(corel2)])
mean(corel2[upper.tri(corel2)])
median(corel2[upper.tri(corel2)])

hc = hclust(as.dist(1-corel2))
corel2 = corel2[hc$order,hc$order]

pdf("007-figures/FigS8_distri_correl_ddGmut.pdf")
hist(corel[upper.tri(corel)],xlab="Pearson correlation coefficient",main="all bZIPs")
hist(corel2[upper.tri(corel2)],xlab="Pearson correlation coefficient",main="without low interaction bZIPs")
dev.off()

ra = max(abs(corel2), na.rm=T)
ra = ceiling(ra*100)/100

pdf("007-figures/Fig2E-heatmap_global_ddGmut.pdf",height=8, width=8)
par(mar=c(6,6,6,6))

image(1:nrow(corel2), 1:ncol(corel2),corel2, col=colPalDiff(ra*200), zlim = c(-ra,ra), , axes=F, xlab="", ylab="")
box()
axis(1,at=1:ncol(corel2),labels = colnames(corel2),las=2, cex.axis=0.8)
axis(2,at=1:ncol(corel2),labels = colnames(corel2),las=2, cex.axis=0.8)

image(seq(-ra,ra,0.01),1,matrix(seq(-ra,ra,0.01)), col=colPalDiff(200*ra), zlim=c(-ra,ra), axes=F, xlab="mutation effect", ylab="",main="color key")
axis(1,at = seq(-1,1,0.1))
box()

plot(as.dendrogram(hc))
dev.off()



################################################################################
# Fig 2F - Example specificity plot for BATF3
################################################################################

co = 1:8
names(co) = c(letters[1:7],"x")

batf3 = mo[["BATF3"]]

lin = batf3$linear_weights
batf3 = batf3$predictions

batf3 = batf3[batf3$i1 > 10 & batf3$i2 > 10 & batf3$i3 > 10 & batf3$o1 > 0 & batf3$o2 > 0 & batf3$o3 > 0 & batf3$expected_bait,]

xl = range(batf3$dG,na.rm=T)
yl = range(batf3$fitness,na.rm=T)
s = seq(xl[1],xl[2],0.01)


pdf("007-figures/Fig2F_spePlots_BATF3.pdf",height=8,width=8)
par(mfrow = c(2,1),mar=c(4,4,4,4))

x = batf3[batf3$wt_prey == "BACH1",]
plot(x$dG, x$fitness, xlim=xl, ylim=yl, type="n", xlab="dG", ylab="binding score",main="BACH1")
text(x$dG, x$fitness, paste0(x$pos,x$mut_aa), col=co[substr(x$id_bait,2,2)], cex = 0.5)

lines(s,lin[1]/(1+exp(s))+lin[2])

abline(v = x$dG[x$id_bait == "0xx"], lty=3)
abline(h = x$fitness[x$id_bait == "0xx"],lty=3)


x = batf3[batf3$wt_prey == "JUN",]
plot(x$dG, x$fitness, xlim=xl, ylim=yl, type="n", xlab="dG", ylab="binding score",main="BACH1")
text(x$dG, x$fitness, paste0(x$pos,x$mut_aa), col=co[substr(x$id_bait,2,2)], cex = 0.5)

lines(s,lin[1]/(1+exp(s))+lin[2])

abline(v = x$dG[x$id_bait == "0xx"], lty=3)
abline(h = x$fitness[x$id_bait == "0xx"],lty=3)


dev.off()



################################################################################
# Fig 3; Fig S9 - Specificity matrices
################################################################################

um = paste0(rep(all_hep,each=140),rep(all_pos,each=20),all_mut)
mut_prey = as.data.frame(as.matrix(expand.grid(all_preys,um)))
names(mut_prey) = c("wt_prey","id_bait")
mut_prey = mut_prey[order(mut_prey$wt_prey, mut_prey$id_bait),]

zl = ceiling(max(abs(he$logFC-he$ind_model_pred), na.rm=T))

he$spe = he$logFC - he$ind_model_pred


pdf("007-figures/Fig3A_S9_specificity_maps.pdf",width=16,height=8)
par(mar=c(4,4,4,4))
lapply(all_baits, function(x){
  
  x = he[he$wt_bait == x,]
  x = x[x$id_bait != "0xx",]
  
  a = mut_prey
  a$spe = x$spe[match(paste(a$wt_prey, a$id_bait), paste(x$wt_prey, x$id_bait))]
  
  m = matrix(a$spe, ncol=length(all_preys))
  colnames(m) = all_preys
  
  m = m[,hc2$labels[hc2$order]] # order based on clustering of Fig 1D
  
  image(1:nrow(m),1:ncol(m),m, axes=F, col=colPalDiff2(zl*2), zlim=c(-zl,zl),xlab="",ylab="", main=x$wt_bait[1])
  axis(1,at=seq(0.5,700.5,20), labels = NA)
  axis(1,at=seq(10.5,690.5,20), tick=F, labels = all_hp)
  axis(2, at=1:ncol(m), labels=colnames(m), las=2, cex.axis=0.5)
  box()

})
dev.off()

pdf("007-figures/Fig3A_legend.pdf")
tmp = matrix(-zl:zl)
image(-zl:zl,1:ncol(tmp),tmp,col=colPalDiff2(zl*2), zlim=c(-zl,zl), xlab="specific mutation effect",ylab="",axes=F)
axis(1, at=seq(-10,10,2), labels=seq(-10,10,2))
dev.off()


pdf("007-figures/Fig3B_distri_spe.pdf")
hist(he$spe,xlab="specific mutation effect",main="")
abline(v=c(-6,6))
dev.off()


nrow(he[he$spe > 6 & !is.na(he$spe),])
nrow(he[he$spe < -6 & !is.na(he$spe),])

he$pos = he$spe > 6
he$neg = he$spe < -6

ret = aggregate(he[,c("pos","neg")], list(he$wt_bait, substr(he$id_bait,1,2)), sum, na.rm=T)
ret = ret[ret[,2] != "0x",]

# we need to add position 5g for ATF1 and CREB1 to be able to make a matrix, since these are one aa shorter
ret = rbind(ret, c("ATF1", "5g", 0, 0),c("CREB1", "5g", 0, 0))
ret = ret[order(ret[,1], ret[,2]),]
ret$pos = as.numeric(ret$pos)
ret$neg = as.numeric(ret$neg)

# matrix for positions with at least one positive specific effect
m = matrix(ret$pos > 0, ncol=53)
colnames(m) = all_baits
m = m[,hc1$labels[hc1$order]] # order based on clustering of Fig 1D


# matrix for positions with at least one negative specific effect
m2 = matrix(ret$neg > 0, ncol=53)
colnames(m2) = unique(ret[,1])
m2 = m2[,hc1$labels[hc1$order]]


pdf("007-figures/Fig3C_specificity_matrix.pdf", width=8,height=8)
par(mar=c(4,4,4,4))
m3 = m + (m2+1)*2 # allows encoding pos, neg, dual
image(1:nrow(m3),1:ncol(m3),m3, axes=F,xlab="",ylab="", col=c("white","blue","red","violet"))
axis(1,at=seq(0.5,35.5,1), labels = NA)
axis(1,at=1:35, tick=F, labels = all_hp, cex.axis=0.5)
axis(2, at=1:ncol(m3), labels=colnames(m3), las=2, cex.axis=0.5)
box()
dev.off()


pos_det_per_pos = rowSums(m2)
names(pos_det_per_pos) = all_hp
pos_det_per_heptad_pos = aggregate(pos_det_per_pos, list(substr(names(pos_det_per_pos),2,2)),sum)
pos_det_per_heptad = aggregate(pos_det_per_pos, list(substr(names(pos_det_per_pos),1,1)),sum)


neg_det_per_pos = rowSums(m)
names(neg_det_per_pos) = all_hp
neg_det_per_heptad_pos = aggregate(neg_det_per_pos, list(substr(names(neg_det_per_pos),2,2)),sum)
neg_det_per_heptad = aggregate(neg_det_per_pos, list(substr(names(neg_det_per_pos),1,1)),sum)


pdf("007-figures/Fig3D_barplot_heptad.pdf")
par(mfrow=c(2,1))
barplot(pos_det_per_heptad$x, names.arg = pos_det_per_heptad$Group.1, col="red", main="positive determinants")
barplot(neg_det_per_heptad$x, names.arg = neg_det_per_heptad$Group.1, col="blue", main="negative determinants")
dev.off()

pdf("007-figures/Fig3E_barplot_heptad_position.pdf")
par(mfrow=c(2,1))
barplot(pos_det_per_heptad_pos$x, names.arg = pos_det_per_heptad_pos$Group.1, col="red", main="positive determinants")
barplot(neg_det_per_heptad_pos$x, names.arg = neg_det_per_heptad_pos$Group.1, col="blue", main="negative determinants")
dev.off()



################################################################################
# Fig 4A - Comparisons of total ddG
################################################################################

restricted_he = he[he$proj_dG > -3 & he$proj_dG < 3 & he$proj_dG_wt > -3 & he$proj_dG_wt < 3 & !is.na(he$ddG_tot) & he$id_bait != "0xx",]
uid = unique(paste(restricted_he$std_id_bait, restricted_he$wt_prey)[!is.na(restricted_he$ddG_tot)])

# build a matrix of ddG_tot for all bZIP
# we do this only for informative bZIPs, i.e. those with a Mochi R2 > 0.4
proj2 = do.call("cbind",lapply(all_baits[!(all_baits %in% bad_bZ)],function(x){
  
  x = restricted_he[restricted_he$wt_bait == x,]
  x$ddG_tot[match(uid,paste(x$std_id_bait,x$wt_prey))]
  
}))
rownames(proj2) = uid
colnames(proj2) = all_baits[!(all_baits %in% bad_bZ)]

co = cor(proj2,use="pairwise.complete.obs")

hc = hclust(as.dist(1-co))
co = co[hc$order, hc$order]


pdf("007-figures/Fig4A_C_correlation_MAFF_G_BATF2_3.pdf")
tmp = proj2[,c("MAFF","MAFG")]
tmp = tmp[!is.na(rowSums(tmp)),]
plot(tmp[,1],tmp[,2],pch=".",xlab="Total ddG MAFF", ylab="Total ddG MAFG")
text(-2,4,paste("R =", sprintf("%.3f",co["MAFF","MAFG"])))

tmp = proj2[,c("BATF2","BATF3")]
tmp = tmp[!is.na(rowSums(tmp)),]
plot(tmp[,1],tmp[,2],pch=".",xlab="Total ddG BATF2", ylab="Total ddG BATF3")
text(-2,4,paste("R =", sprintf("%.3f",co["BATF2","BATF3"])))

tmp = proj2[,c("MAFF","BATF2")]
tmp = tmp[!is.na(rowSums(tmp)),]
plot(tmp[,1],tmp[,2],pch=".",xlab="Total ddG MAFF", ylab="Total ddG BATF2")
text(-2,4,paste("R =", sprintf("%.3f",co["MAFF","BATF2"])))
dev.off()


pdf("007-figures/Fig4D_correlation_matrix_tot_ddG.pdf",width=14,height=14)
par(mar=c(6,6,6,6))

image(1:nrow(co),1:ncol(co),co, axes=F, xlab="", ylab="",col=colPalDiff(100),zlim=c(-1,1))
box()
axis(1,at=1:ncol(co),labels = colnames(co),las=2)
axis(2,at=1:ncol(co),labels = colnames(co),las=2)


plot(as.dendrogram(hc))

image(as.matrix(seq(-1,1,0.02)),col=colPalDiff(100),zlim=c(-1,1),axes=F)
box()
dev.off()



################################################################################
# Fig 5B-D - CNN models
################################################################################

random = split(dl_pred, dl_pred$fold_random)
id = split(dl_pred, dl_pred$fold_id)
aa = split(dl_pred, dl_pred$fold_aa)
pos = split(dl_pred, dl_pred$fold_pos)
bZ = split(dl_pred, dl_pred$fold_bZ)

r2_random = unlist(lapply(random,function(x){
  1 - sum((x$norm_logFC_random-x$valid_pred_random)^2) / sum((x$norm_logFC_random - mean(x$norm_logFC_random))^2)
}))
r2_id = unlist(lapply(id,function(x){
  1 - sum((x$norm_logFC_id-x$valid_pred_id)^2) / sum((x$norm_logFC_id - mean(x$norm_logFC_id))^2)
}))
r2_aa = unlist(lapply(aa,function(x){
  1 - sum((x$norm_logFC_aa-x$valid_pred_aa)^2) / sum((x$norm_logFC_aa - mean(x$norm_logFC_aa))^2)
}))
r2_pos = unlist(lapply(pos,function(x){
  1 - sum((x$norm_logFC_pos-x$valid_pred_pos)^2) / sum((x$norm_logFC_pos - mean(x$norm_logFC_pos))^2)
}))
r2_bZ = unlist(lapply(bZ,function(x){
  1 - sum((x$norm_logFC_bZ-x$valid_pred_bZ)^2) / sum((x$norm_logFC_bZ - mean(x$norm_logFC_bZ))^2)
}))

r2 = cbind(r2_random, r2_id, r2_aa, r2_pos, r2_bZ)
means = colMeans(r2)
sds = apply(r2,2,sd)

pdf("007-figures/Fig5B_predicted_vs_measured_plot_random_split.pdf",height=8, width=9)
par(mar=c(4,4,4,4))
ggplot(dl_pred, aes(x=valid_pred_random, y=norm_logFC_random) ) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis", trans="log", breaks=10^(0:6)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()

# overall R2 for the random split
1 - sum((dl_pred$norm_logFC_random-dl_pred$valid_pred_random)^2) / sum((dl_pred$norm_logFC_random - mean(dl_pred$norm_logFC_random))^2)

pdf("007-figures/Fig5C_barplot_CNN.pdf")
barplot2(means[5:1], horiz=T, ci.l = means[5:1]-sds[5:1]/sqrt(10)*1.96, ci.u = means[5:1]+sds[5:1]/sqrt(10)*1.96, plot.ci=T, col="#441853")
dev.off()


r2_cnn_aa = 1 - sum((dl_pred$norm_logFC_aa_singleFold-dl_pred$valid_pred_cnn_aa_singleFold)^2, na.rm=T) / sum((dl_pred$norm_logFC_aa_singleFold - mean(dl_pred$norm_logFC_aa_singleFold, na.rm=T))^2, na.rm=T)
r2_esm_aa = 1 - sum((dl_pred$norm_logFC_aa_singleFold-dl_pred$valid_pred_esm_aa_singleFold)^2, na.rm=T) / sum((dl_pred$norm_logFC_aa_singleFold - mean(dl_pred$norm_logFC_aa_singleFold, na.rm=T))^2, na.rm=T)
r2_esm_bZ = 1 - sum((dl_pred$norm_logFC_esm_bZ_singleFold-dl_pred$valid_pred_esm_bZ_singleFold)^2, na.rm=T) / sum((dl_pred$norm_logFC_esm_bZ_singleFold - mean(dl_pred$norm_logFC_esm_bZ_singleFold, na.rm=T))^2, na.rm=T)


pdf("007-figures/Fig5D_barplot_CNN_ESM.pdf")
barplot2(c(r2_esm_bZ,r2_esm_aa,r2_cnn_aa), names.arg = c("ESM bZ", "ESM aa", "CNN aa"),horiz=T, col="#441853")
dev.off()



################################################################################
# Fig 5E - synthetic bZIPs
################################################################################


b = sort(unique(lfc$bait)[grep("targeting",unique(lfc$bait))])
p = unique(lfc$prey)
all_synth = as.data.frame(as.matrix(expand.grid(b,p)))
names(all_synth) = c("bait","prey")
all_synth = merge(all_synth, lfc, by = names(all_synth), all.x=T)
m_obs = t(matrix(all_synth$logFC, nrow=54)) # rows are baits
m_pred = t(matrix(all_synth$pred_logFC, nrow=54)) # rows are baits

targets = do.call("rbind",strsplit(unique(all_synth$bait),"_"))[,1]


pdf("007-figures/Fig5E_heatmap_synth_bZ.pdf", width=nrow(m_obs)/5, height = ncol(m_obs)/5)
par(mar=c(4,4,4,4))

ra = range(m_pred,na.rm=T)
ra[1] = floor(ra[1]*10)/10
ra[2] = ceiling(ra[2]*10)/10
image(1:nrow(m_pred), 1:ncol(m_pred), m_pred, col=colPal((ra[2]-ra[1])*10), axes=F, xlab="", ylab="")
axis(1, at = seq(0.5,nrow(m_pred)+0.5,10), labels = NA)
axis(1, at = seq(5.5,nrow(m_pred)-4.5,10), tick = F, labels = unique(targets))
axis(2, at = 1:ncol(m_pred), tick=F, unique(all_synth$prey),las=2, cex.axis=0.7)
box()

m = matrix(seq(ra[1]*10,ra[2]*10,1))
image(seq(ra[1],ra[2],0.1),1,matrix(seq(ra[1]*10,ra[2]*10,1)), col=colPal((ra[2]-ra[1])*10), axes=F, xlab="measured binding score", ylab="",main="color key")
axis(1,at = seq(ra[1],ra[2],2.5), labels=seq(ra[1],ra[2],2.5))
box()

ra = range(m_obs,na.rm=T)
ra[1] = floor(ra[1]*10)/10
ra[2] = ceiling(ra[2]*10)/10
image(1:nrow(m_obs), 1:ncol(m_obs), m_obs, col=colPal((ra[2]-ra[1])*10), axes=F, xlab="", ylab="")
axis(1, at = seq(0.5,nrow(m_obs)+0.5,10), labels = NA)
axis(1, at = seq(5.5,nrow(m_obs)-4.5,10), tick = F, labels = unique(targets))
axis(2, at = 1:ncol(m_obs), tick=F, unique(all_synth$prey),las=2, cex.axis=0.7)
box()

m = matrix(seq(ra[1]*10,ra[2]*10,1))
image(seq(ra[1],ra[2],0.1),1,matrix(seq(ra[1]*10,ra[2]*10,1)), col=colPal((ra[2]-ra[1])*10), axes=F, xlab="measured binding score", ylab="",main="color key")
axis(1,at = seq(-2.5,ra[2],2.5), labels=seq(-2.5,ra[2],2.5))
box()

dev.off()
