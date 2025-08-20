library(stringr)
library(parallel)

load("../../003-lfc_limma_AA_per_batch.Rdata")

bZIP <- sort(unique(lfc_he_aa$wt_bait))
bZIP <- data.frame(bZIP = bZIP, dummy = rep(NA, length(bZIP)))

lfc_he_aa = lfc_he_aa[lfc_he_aa$wt_bait != "ATF7",]


# make dummy sequence to encode wild-type partners
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

ref = which(lfc_he_aa$id_bait == "0xx" & lfc_he_aa$wt_bait == "JDP2" & lfc_he_aa$full_id_prey == "DDIT3_x0x")
n = stringdist::stringdist(lfc_he_aa$aa_seq[ref], lfc_he_aa$aa_seq, method = "hamming")

all_variants = data.frame(aa_seq = lfc_he_aa$aa_seq,
                          Nham_aa = n,
                          WT = FALSE,
                          fitness = lfc_he_aa$logFC,
                          sigma =lfc_he_aa$se.logFC
)
all_variants$WT[ref] = TRUE

save(all_variants, file = "001-Mochi_input.RData")

conv = lfc_he_aa[,c("full_id_bait","full_id_prey","aa_seq")]
save(conv, file="001-dummy_conversion.Rdata")
