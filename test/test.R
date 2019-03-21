#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)!=1) {
  stop("one argument must be supplied (linux/crossplatform).n", call.=FALSE)
} 

library(MAGMA.Celltyping)
setwd('/Volumes/botweiser/microglia_iPSC/MCT_test/')
gwas_1 = 'AD_sumstats_Jansen_test.txt'
gwas_2 = 'body_BMIz_test.txt'
gwas_3 = 'body_HEIGHTz_test.txt'
gwas_4 = 'CLOZUK_SCZ_test.txt'
gwas_5 = 'IGAP_test.txt'
gwas_6 = 'PGC_SCZ2_test.txt'
genome_ref_dir = "/Volumes/botweiser/reference/bed/g1000_eur"
genome_ref_path = sprintf("%s/g1000_eur", genome_ref_dir)


if (args[1]=='linux') {
  col_headers = format_sumstats_for_magma_linux(gwas_1)
  col_headers = format_sumstats_for_magma_linux(gwas_2)
  col_headers = format_sumstats_for_magma_linux(gwas_3)
  col_headers = format_sumstats_for_magma_linux(gwas_4)
  col_headers = format_sumstats_for_magma_linux(gwas_5)
  col_headers = format_sumstats_for_magma_linux(gwas_6)
} else if (args[1]=='crossplatform') {
  col_headers = format_sumstats_for_magma_crossplatform(gwas_1)
  col_headers = format_sumstats_for_magma_crossplatform(gwas_2)
  col_headers = format_sumstats_for_magma_crossplatform(gwas_3)
  col_headers = format_sumstats_for_magma_crossplatform(gwas_4)
  col_headers = format_sumstats_for_magma_crossplatform(gwas_5)
  col_headers = format_sumstats_for_magma_crossplatform(gwas_6)
} else {stop("Incorrect argument specified.n", call.=FALSE)}

