#!/usr/bin/env bash
# Run from the directory containing SparK.py and the ENCODE example inputs.
set -euo pipefail

python SparK.py \
-pr chr12:6520512-6640512 \
-cf HepG2_H3K27AC_1_ENCFF495QSO.bigWig.bdg.bw HepG2_H3K27AC_2_ENCFF348RLL.bigWig.bdg.bw HepG2_H3K4me3_1_ENCFF699DRO.bigWig.bdg.bw HepG2_H3K4me3_2_ENCFF400FYO.bigWig.bdg.bw \
-gtf gencode.v34.chr_patch_hapl_scaff.annotation.gtf \
-gl H3K27AC H3K4me3 H3K27AC-2 H3K4me3-2 \
-dg GAPDH IFFO1 NOP2 CHD4 LPAR5 \
-o 1_bw

python SparK.py \
-pr chr12:6520512-6640512 \
-cf HepG2_H3K27AC_1_ENCFF495QSO.bigWig.bdg.bw HepG2_H3K27AC_2_ENCFF348RLL.bigWig.bdg.bw HepG2_H3K4me3_1_ENCFF699DRO.bigWig.bdg.bw HepG2_H3K4me3_2_ENCFF400FYO.bigWig.bdg.bw \
-gtf gencode.v34.chr_patch_hapl_scaff.annotation.gtf \
-gl H3K27AC H3K4me3 H3K27AC-2 H3K4me3-2 \
-dg GAPDH IFFO1 NOP2 CHD4 LPAR5 \
-f 000000 \
-o 1_custom_color_bw

python SparK.py \
-pr chr12:6527512-6550512 \
-tf K562_H3K27AC_1_ENCFF779QTH.bigWig.bdg.bw K562_H3K27AC_2_ENCFF945XHA.bigWig.bdg.bw K562_H3K4me3_1_ENCFF804OLI.bigWig.bdg.bw K562_H3K4me3_2_ENCFF352VRB.bigWig.bdg.bw \
-cf HepG2_H3K27AC_1_ENCFF495QSO.bigWig.bdg.bw HepG2_H3K27AC_2_ENCFF348RLL.bigWig.bdg.bw HepG2_H3K4me3_1_ENCFF699DRO.bigWig.bdg.bw HepG2_H3K4me3_2_ENCFF400FYO.bigWig.bdg.bw \
-tg 1 1 2 2 \
-cg 1 1 2 2 \
-gl H3K27AC H3K4me3 \
-l HepG2_cells K562_cells \
-gtf gencode.v34.chr_patch_hapl_scaff.annotation.gtf \
-ps averages \
-dg GAPDH \
-o 2_averages_bw

python SparK.py \
-pr chr12:6527512-6550512 \
-tf K562_H3K27AC_1_ENCFF779QTH.bigWig.bdg.bw K562_H3K27AC_2_ENCFF945XHA.bigWig.bdg.bw K562_H3K4me3_1_ENCFF804OLI.bigWig.bdg.bw K562_H3K4me3_2_ENCFF352VRB.bigWig.bdg.bw \
-cf HepG2_H3K27AC_1_ENCFF495QSO.bigWig.bdg.bw HepG2_H3K27AC_2_ENCFF348RLL.bigWig.bdg.bw HepG2_H3K4me3_1_ENCFF699DRO.bigWig.bdg.bw HepG2_H3K4me3_2_ENCFF400FYO.bigWig.bdg.bw \
-tg 1 1 2 2 \
-cg 1 1 2 2 \
-gl H3K27AC H3K4me3 \
-l HepG2_cells K562_cells \
-gtf gencode.v34.chr_patch_hapl_scaff.annotation.gtf \
-ps averages \
-f CE358B 005CFF \
-bed bedtest.txt bedtest.txt \
-bedcol EF1414 00B7FF \
-o 3_bw

python SparK.py \
-pr chr12:6527512-6550512 \
-tf K562_H3K27AC_1_ENCFF779QTH.bigWig.bdg.bw K562_H3K27AC_2_ENCFF945XHA.bigWig.bdg.bw K562_H3K4me3_1_ENCFF804OLI.bigWig.bdg.bw K562_H3K4me3_2_ENCFF352VRB.bigWig.bdg.bw \
-cf HepG2_H3K27AC_1_ENCFF495QSO.bigWig.bdg.bw HepG2_H3K27AC_2_ENCFF348RLL.bigWig.bdg.bw HepG2_H3K4me3_1_ENCFF699DRO.bigWig.bdg.bw HepG2_H3K4me3_2_ENCFF400FYO.bigWig.bdg.bw \
-tg 1 1 2 2 \
-cg 1 1 2 2 \
-gl H3K27AC H3K4me3 \
-l HepG2_cells K562_cells \
-gtf gencode.v34.chr_patch_hapl_scaff.annotation.gtf \
-f all_grey \
-sp yes \
-o 4_bw

python SparK.py \
-pr chr12:6533612-6539012 \
-tf K562_H3K27AC_1_ENCFF779QTH.bigWig.bdg.bw K562_H3K27AC_2_ENCFF945XHA.bigWig.bdg.bw K562_RNAseq_plus_1_ENCFF682BLJ.bigWig.bdg.bw K562_RNAseq_plus_2_ENCFF745PIR.bigWig.bdg.bw \
-cf HepG2_H3K27AC_1_ENCFF495QSO.bigWig.bdg.bw HepG2_H3K27AC_2_ENCFF348RLL.bigWig.bdg.bw HepG2_RNAseq_plus_1_ENCFF576COO.bigWig.bdg.bw HepG2_RNAseq_plus_2_ENCFF599HGT.bigWig.bdg.bw \
-tg 1 1 2 2 \
-cg 1 1 2 2 \
-gl H3K27AC RNA-seq \
-l HepG2_cells K562_cells \
-gtf gencode.v34.chr_patch_hapl_scaff.annotation.gtf \
-bed bedtest.txt bedtest.txt \
-bedcol EF1414 00B7FF \
-dg GAPDH \
-o 5_bw

python SparK.py \
-pr chr12:6540512-6545512 \
-tf K562_H3K27AC_1_ENCFF779QTH.bigWig.bdg.bw K562_H3K27AC_2_ENCFF945XHA.bigWig.bdg.bw K562_H3K4me3_1_ENCFF804OLI.bigWig.bdg.bw K562_H3K4me3_2_ENCFF352VRB.bigWig.bdg.bw \
-cf HepG2_H3K27AC_1_ENCFF495QSO.bigWig.bdg.bw HepG2_H3K27AC_2_ENCFF348RLL.bigWig.bdg.bw HepG2_H3K4me3_1_ENCFF699DRO.bigWig.bdg.bw HepG2_H3K4me3_2_ENCFF400FYO.bigWig.bdg.bw \
-tg 1 1 2 2 \
-cg 1 1 2 2 \
-gl H3K27AC H3K4me3 \
-l HepG2_cells K562_cells \
-gtf gencode.v34.chr_patch_hapl_scaff.annotation.gtf \
-pt STD \
-o 6_bw

python SparK.py \
-pr chr12:6533612-6539012 \
-tf K562_H3K27AC_1_ENCFF779QTH.bigWig.bdg.bw K562_H3K27AC_2_ENCFF945XHA.bigWig.bdg.bw K562_H3K4me3_1_ENCFF804OLI.bigWig.bdg.bw K562_H3K4me3_2_ENCFF352VRB.bigWig.bdg.bw \
-cf HepG2_H3K27AC_1_ENCFF495QSO.bigWig.bdg.bw HepG2_H3K27AC_2_ENCFF348RLL.bigWig.bdg.bw HepG2_H3K4me3_1_ENCFF699DRO.bigWig.bdg.bw HepG2_H3K4me3_2_ENCFF400FYO.bigWig.bdg.bw \
-tg 1 1 2 2 \
-cg 1 1 2 2 \
-gl H3K27AC H3K4me3 \
-l HepG2_cells K562_cells \
-gtf gencode.v34.chr_patch_hapl_scaff.annotation.gtf \
-o 7_bw

python SparK.py \
-pr chr12:6527512-6550512 \
-tf K562_H3K27AC_1_ENCFF779QTH.bigWig.bdg.bw K562_H3K27AC_2_ENCFF945XHA.bigWig.bdg.bw K562_H3K4me3_1_ENCFF804OLI.bigWig.bdg.bw K562_H3K4me3_2_ENCFF352VRB.bigWig.bdg.bw \
-cf HepG2_H3K27AC_1_ENCFF495QSO.bigWig.bdg.bw HepG2_H3K27AC_2_ENCFF348RLL.bigWig.bdg.bw HepG2_H3K4me3_1_ENCFF699DRO.bigWig.bdg.bw HepG2_H3K4me3_2_ENCFF400FYO.bigWig.bdg.bw \
-tg 1 1 2 2 \
-cg 1 1 2 2 \
-gl H3K27AC H3K4me3 \
-l HepG2_cells K562_cells \
-gtf gencode.v34.chr_patch_hapl_scaff.annotation.gtf \
-ps averages \
-f CE358B 005CFF \
-bed bedtest.txt bedtest.txt \
-bedcol EF1414 00B7FF \
-bedlab custom1 custom2 \
-o 8_custom_bed_labels_bw

python SparK.py \
-pr chr12:6533612-6539012 \
-cf HepG2_H3K27AC_1_ENCFF495QSO.bigWig.bdg.bw \
-o 9_bw

python SparK.py \
-pr chr12:6533612-6539012 \
-tf K562_H3K27AC_1_ENCFF779QTH.bigWig.bdg.bw K562_H3K27AC_2_ENCFF945XHA.bigWig.bdg.bw K562_H3K4me3_1_ENCFF804OLI.bigWig.bdg.bw K562_H3K4me3_2_ENCFF352VRB.bigWig.bdg.bw \
-cf HepG2_H3K27AC_1_ENCFF495QSO.bigWig.bdg.bw HepG2_H3K27AC_2_ENCFF348RLL.bigWig.bdg.bw HepG2_H3K4me3_1_ENCFF699DRO.bigWig.bdg.bw HepG2_H3K4me3_2_ENCFF400FYO.bigWig.bdg.bw \
-tg 1 1 2 2 \
-cg 1 1 2 2 \
-gl H3K27AC H3K4me3 \
-l HepG2_cells K562_cells \
-gtf gencode.v34.chr_patch_hapl_scaff.annotation.gtf \
-pt sine \
-o 10_sine_bw

python SparK.py \
-pr chr12:6520512-6640512 \
-cf HepG2_H3K27AC_1_ENCFF495QSO.bigWig.bdg.bw HepG2_H3K27AC_2_ENCFF348RLL.bigWig.bdg.bw HepG2_H3K4me3_1_ENCFF699DRO.bigWig.bdg.bw HepG2_H3K4me3_2_ENCFF400FYO.bigWig.bdg.bw \
-gtf gencode.v34.chr_patch_hapl_scaff.annotation.gtf \
-gl H3K27AC H3K4me3 H3K27AC-2 H3K4me3-2 \
-dg GAPDH IFFO1 NOP2 CHD4 LPAR5 \
-cs 60 150 47 50 \
-o 11_1_custom_y_axis_bw

