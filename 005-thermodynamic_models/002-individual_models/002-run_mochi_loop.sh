#!/bin/bash

module use /tungstenfs/groups/gbioinfo/Appz/easybuild/modules/all
module use /tungstenfs/groups/gbioinfo/Appz/modules
module load MoCHI/20240409

xxx_values=("ATF1"    "ATF2"    "ATF3"    "ATF4"    "ATF5"    "ATF6"    "ATF6B"   "BACH1"   "BACH2"  "BATF"    "BATF2"   "BATF3"   "CEBPA"   "CEBPB"   "CEBPD"   "CEBPE"   "CEBPG"   "CREB1"  "CREB3"   "CREB3L1" "CREB3L2" "CREB3L3" "CREB3L4" "CREB5"   "CREBL2"  "CREBRF"  "CREBZF" "CREM"    "DBP"     "DDIT3"   "FOS"     "FOSB"    "FOSL1"   "FOSL2"   "HLF"     "JDP2"   "JUN"     "JUNB"    "JUND"    "MAF"     "MAFA"    "MAFB"    "MAFF"    "MAFG"    "MAFK"   "NFE2"    "NFE2L1"  "NFE2L2"  "NFE2L3"  "NFIL3"   "NRL"     "TEF"     "XBP1")


# Loop through each XXX value
for xxx in "${xxx_values[@]}"; do
  output_directory="002-output/$xxx/"
  design_file="001-input_tables/${xxx}_design.txt"
  stdout="002-output/$xxx/stdout.txt"
  stderr="002-output/$xxx/stderr.txt"

  echo $xxx
  
  run_mochi.py --output_directory "$output_directory" --model_design "$design_file" 1>$stdout 2>$stderr
done