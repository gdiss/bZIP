library(stringr)
library(parallel)


bZ = read.delim("../../000-data/TableS6_wt_bZIPs_aa_seq.txt")
bZ = bZ[order(bZ$Gene),]

# table for vairant id mapping
load("../../003-lfc_limma_AA_per_batch.Rdata")
var_seq = lfc_he_aa$finalBaitAA
str_sub(var_seq,  -1) = ""
lfc_he_aa$finalBaitAA_2 = var_seq
map_var_id = lfc_he_aa[!duplicated(lfc_he_aa$finalBaitAA_2),]
lfc_he_aa$wt_prey = do.call("rbind", strsplit(lfc_he_aa$full_id_prey,"_"))[,1]


# conversion between dummy seq and wt id
tmp = NULL
for(i in 1:(nrow(bZ)-1)){
  x = paste(rep("A", i-1), collapse = "")
  y = paste(rep("A", 53-i), collapse = "")
  tmp = c(tmp, paste0(x, "C", y))
}
w = which(bZ$Gene == "DDIT3")
bZ$dummy = NA
bZ$dummy[1:(w-1)] = tmp[1:(w-1)]
bZ$dummy[(w+1):54] = tmp[w:53]
bZ$dummy[w] = paste(rep("A", 53), collapse = "")



tmp_conv = 28:62
names(tmp_conv) = paste0(rep(1:5,each=7),letters[1:7])

xx = bZ$Gene[bZ$Gene != "ATF7"]
mo = mclapply(xx, function(x){
  
  pr = read.delim(paste0("002-output/",x,"/mochi_project/task_1/predictions/predicted_phenotypes_all.txt"))
  ddG = read.delim(paste0("002-output/",x,"/mochi_project/task_1/weights/weights_Binding.txt"))
  lin = read.delim(paste0("002-output/",x,"/mochi_project/task_1/weights/linears_weights_Binding.txt"))
  
  wt_seq = str_sub(pr$aa_seq, -53, -1)
  var_seq = str_sub(pr$aa_seq, 1, -54)
  
  pr = cbind(pr, map_var_id[match(var_seq,map_var_id$finalBaitAA_2),c("wt_bait","id_bait","expected_bait")])
  pr$wt_prey = bZ$Gene[match(wt_seq, bZ$dummy)]
  pos = as.numeric(tmp_conv[substr(pr$id_bait,1,2)])
  pos[!pr$expected_bait] = NA
  pr$pos = pos
  pr$wt_aa = sapply(pos,function(i){substr(bZ$bZIP_aa[bZ$Gene == x],i,i)})
  pr$mut_aa = substr(pr$id_bait,3,3)
  pr$mut_aa[!pr$expected_bait] = NA
  
  tmp = map_var_id[!duplicated(map_var_id$std_id_bait),]
  
  ddG$Pos[ddG$id_ref == "WT"] = 0
  ddG = ddG[order(ddG$Pos),]
  ddG_bZ = ddG[(nrow(ddG)-52):nrow(ddG),]
  ddG_bZ$id = bZ$Gene[bZ$Gene != "DDIT3"]
  names(ddG_bZ)[1] = "wt_prey"
  ddG_bZ = rbind(ddG_bZ,ddG_bZ[1,])
  ddG_bZ[nrow(ddG_bZ),1] = "DDIT3"
  ddG_bZ[nrow(ddG_bZ),2:15] = NA
  ddG_bZ[nrow(ddG_bZ),16:18] = 0
  ddG_bZ[nrow(ddG_bZ),20:22] = NA
  
  
  ddG = ddG[1:(nrow(ddG)-53),]
  ddG$id = tmp$id_bait[match(ddG$id_ref, tmp$std_id_bait)]
  ddG$id[is.na(ddG$id)] = ""
  ddG$id[ddG$id_ref == "WT"] = "0xx"
  ddG$expected_bait = tmp$expected_bait[match(ddG$id_ref, tmp$std_id_bait)]
  ddG$expected_bait[ddG$id_ref == "WT"] = TRUE
  ddG$expected_bait[is.na(ddG$expected_bait)] = FALSE
  
  k = match(paste(pr$wt_bait, pr$id_bait, pr$wt_prey), paste(lfc_he_aa$wt_bait, lfc_he_aa$id_bait, lfc_he_aa$wt_prey))
  pr$i1 = lfc_he_aa$i1[k]
  pr$i2 = lfc_he_aa$i2[k]
  pr$i3 = lfc_he_aa$i3[k]
  pr$o1 = lfc_he_aa$o1[k]
  pr$o2 = lfc_he_aa$o2[k]
  pr$o3 = lfc_he_aa$o3[k]
  pr$ddG_mut = ddG$mean[match(pr$id_bait, ddG$id)]
  pr$ddG_mut[pr$id_bait == "0xx"] = 0
  pr$ddG_bZ = ddG_bZ$mean[match(pr$wt_prey, ddG_bZ$wt_prey)]
  pr$dG = ddG$mean[ddG$id == "0xx"] + pr$ddG_bZ + pr$ddG_mut
  
  lin = c(mean(lin$kernel),mean(lin$bias))
  
  list(predictions=pr,ddG_mut=ddG,ddG_bZ=ddG_bZ,linear_weights=lin)
  
},mc.cores=24)
names(mo) = xx


save(mo, file="003-Mochi_output.Rdata")



