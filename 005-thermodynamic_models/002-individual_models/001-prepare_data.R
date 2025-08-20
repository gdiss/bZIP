library(stringr)
library(parallel)

load("../../003-lfc_limma_AA_per_batch.Rdata")

bZIP <- sort(unique(lfc_he_aa$wt_bait))
bZIP <- data.frame(bZIP = bZIP, dummy = rep(NA, length(bZIP)))

lfc_he_aa = lfc_he_aa[lfc_he_aa$wt_bait != "ATF7",]


tmp = NULL
for(i in 1:(nrow(bZIP)-1)){
  x = paste(rep("A", i-1), collapse = "")
  y = paste(rep("A", 53-i), collapse = "")
  tmp = c(tmp, paste0(x, "C", y))
}

w = which(bZIP$bZIP == "DDIT3")
bZIP$dummy[1:(w-1)] = tmp[1:(w-1)]
bZIP$dummy[(w+1):54] = tmp[w:53]
bZIP$dummy[w] = paste(rep("A", 53), collapse = "")


var_seq = lfc_he_aa$finalBaitAA
str_sub(var_seq,  -1) = ""

lfc_he_aa$var_seq = var_seq
lfc_he_aa = lfc_he_aa[!grepl("\\*",lfc_he_aa$var_seq),]

tmp = do.call("rbind",strsplit(lfc_he_aa$full_id_prey,"_"))[,1]
lfc_he_aa$prey = bZIP$dummy[match(tmp, bZIP$bZIP)]
lfc_he_aa$aa_seq = paste0(lfc_he_aa$var_seq, lfc_he_aa$prey)

d = split(lfc_he_aa, lfc_he_aa$wt_bait)

lapply(d,function(x){
  
  cat(x$wt_bait[1],"\n")
  ref = which(x$id_bait == "0xx" & x$full_id_prey == "DDIT3_x0x")
  n = stringdist::stringdist(x$aa_seq[ref], x$aa_seq, method = "hamming")
  
  all_variants = data.frame(nt_seq = NA,
             aa_seq = x$aa_seq,
             Nham_nt = NA,
             Nham_aa = n,
             WT = FALSE,
             fitness = x$logFC,
             sigma = x$se.logFC
             )
  all_variants$WT[ref] = TRUE
  
  save(all_variants, file = paste0("001-input_tables/",x$wt_bait[1],".RData"))
})


lapply(1:length(d),function(x){
  txt = paste0("phenotype\ttransformation\ttrait	file\nBinding\tTwoStateFractionFolded\tBinding\t001-input_tables/",names(d)[x],".RData")
  out = paste0("001-input_tables/",names(d)[x],"_design.txt")
  write(txt, file = out)
  
  output_dir = paste0("002-output/",bZIP$bZIP[i])
  dir.create(output_dir, showWarnings = F)
})

