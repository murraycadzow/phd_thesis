panel <- read.delim(paste0('~/data/NZ_coreExome_1kgp/nz_1kgp.panel'), stringsAsFactors = FALSE) %>% filter(!pop %in% c('NAD','WPN','POL','EPN'))
markers <- read.table('~/data/NZ_coreExome_1kgp/nz_1kg_markers.txt', stringsAsFactors = FALSE, header = FALSE)
names(markers) <- c("chrom","chrom_start","marker","ref","alt")
gwas_cat <- read.delim('~/data/gwas_catalog//gwas_catalog_v1.0.1-associations_e89_r2017-06-19.tsv', header=TRUE, stringsAsFactors = FALSE, sep='\t')


# load the significant data
## frequency based
lower_sig_stats <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_lower_sig_stat_genes_filtered_ns3_resized-31-8-2017.RDS')
lower_sig_stats <- bind_rows(lapply(names(lower_sig_stats), function(y){
  lower_sig_stats[[y]] <- bind_rows(
    lapply(names(lower_sig_stats[[y]]), function(x){
      lower_sig_stats[[y]][[x]] %>% mutate(pop = x, statname =y, chrom = as.character(chrom)) 
    }
    ))
})) %>% filter(!pop %in% c("NAD","EPN","WPN","POL"))
lower_sig_stats_binned <-readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_lower_sig_stat_filtered_ns3_resized_binned-2-11-2017.RDS') %>% mutate(chrom = paste0('chr',chrom)) %>% select('statname' = stat, everything())



upper_sig_stats <-readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_upper_sig_stat_genes_filtered_ns3_resized-31-8-2017.RDS')
upper_sig_stats <- bind_rows(lapply(names(upper_sig_stats), function(y){
  upper_sig_stats[[y]] <- bind_rows(
    lapply(names(upper_sig_stats[[y]]), function(x){
      upper_sig_stats[[y]][[x]] %>% mutate(pop = x, statname =y, chrom = as.character(chrom)) 
    }
    ))
})) %>% filter(!pop %in% c("NAD","EPN","WPN","POL"))

upper_sig_stats_binned <-readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_upper_sig_stat_filtered_ns3_resized_binned-2-11-2017.RDS') %>% mutate(chrom = paste0('chr',chrom)) %>% select('statname' = stat, everything())

## permutation data
# reorder the lower and upper columns so that the numbers make sense in the tables
fdr_low <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_chr_stats_popgenome_long_ns3_resized_binned_fdr_low.2-11-2017.RDS') %>% select(stat, pop, lower_per, 'lowerCI' = upper, 'FDR' = med, 'upperCI' = lower, lowest_quant, super_pop)

fdr_high <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_chr_stats_popgenome_long_ns3_resized_binned_fdr_high.2-11-2017.RDS') %>% select(stat, pop, upper_per, 'lowerCI' = lower, 'FDR' = med, 'upperCI' = upper, highest_quant, super_pop)

perm_median <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/permutations_median_thresholds.2-11-2017.RDS')
perm_ci2.5 <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/permutations_ci2.5_thresholds.2-11-2017.RDS')
perm_ci97.5 <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/permutations_ci97.5_thresholds.2-11-2017.RDS')


## load global summary info
global_summary <- readRDS("~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_summary_popgenome_filtered_3ns_resized.31-8-2017.RDS") %>% filter(!pop %in% c("WPN","EPN","POL","NAD"))
global_super_summary <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_summary_superpop_popgenome_filtered_3ns_resized.6-11-2017.RDS')


## load ihs and nsl data
ihs_clus_regions <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/ihs_clus_regions-15-12-2017.RDS')
ihs_clus <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/sig_ihs_clus-15-12-2017.RDS')
sig_ihs_val <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/sig_ihs_values_with_genes_14-11-2017.RDS') %>% filter(abs(statvalue) > 2.6 )

nsl_clus_regions <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/nsl_clus_regions-15-12-2017.RDS')
nsl_clus <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/sig_nsl_clus-15-12-2017.RDS')
sig_nsl_val <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/sig_nsl_values_with_genes_14-11-2017.RDS')  %>% filter(abs(statvalue) > 2.6 )


#nsl_clus <-readRDS('~/data/NZ_coreExome_1kgp/haplotype/nsl_clus_regions-14-7-2017.RDS')

#100kb window 10kb slide pairwise fst
# FST can't be negative. This keeps NAs
pol_fst <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/windowed_poly_chr_fst_popgenome_ns3.14-11-2017.RDS') %>% mutate_at(.vars = vars(-contains('chrom')), function(x){ifelse(x > 0, x, 0)})
pol_fst_99 <- pol_fst %>%  gather('pop', "fst", -contains('chrom')) %>% group_by(pop) %>% filter(fst > quantile(fst, 0.99, na.rm = TRUE))
#source('~/Git_repos/bookdown_thesis/scripts/Annotate_genes.R')
#pol_fst_99genes <- txdb_gene_annotate(pol_fst_99 %>% mutate(chrom = paste0('chr',chrom)) %>% dplyr::select(contains("chrom")) %>% GenomicRanges::GRanges())
pol_fst_99genes <- readRDS("~/data/NZ_coreExome_1kgp/100kbWindow_intra/windowed_poly_chr_fst_popgenome_ns3_99per_genes.25-1-2018.RDS") 

# load xpehh
xpehh <- read_delim('~/data/NZ_coreExome_1kgp/haplotype/sig_xpehh_29-10-2017.csv', delim = ',') %>% filter(!pop1 %in% c("POL","EPN","WPN", "NAD"), !pop2 %in% c("POL","EPN","WPN", "NAD")) %>% filter(abs(xpehh_value) > 2.6)

# find all of the genes from POL that have something to suggest they may have been selected
table(unique(lower_sig_stats$SYMBOL) %in% (lower_sig_stats %>% filter(statname %in% c('Tajima.D','Fay.Wu.H','Zeng.E', 'Fu.Li.F', 'Fu.Li.D'), pop %in% c('CIM','NZM','SAM','TON')) %>% .[['SYMBOL']] %>% unique()))

load('~/data/gwas_catalog/diseaseGR-25-7-2017.RData')# brings in objects called {gc_urate_gout,urate_gout,kd,metsyn,obesity,t2d}_GR
load('~/data/gwas_catalog/zhang_and_extra_18-11-2017.RData') # brings in objects called malaria_GR, zhang_immune_GR, neurological_GR, psychiatric_GR

# made from Thesis/selectionpipeline/consec_regions.R
consec_lower_regions <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_lower_sig_consec_regions_annotated-15-12-2017.RDS') %>% filter(!pop %in% c('WPN','EPN','NAD'))
# this is incorrect
consec_upper_regions <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_upper_sig_consec_regions_annotated-15-12-2017.RDS') %>% filter(!pop %in% c('WPN','EPN','NAD'))


plot_gene <- function(p, g){
  wh <- genesymbol[g]
  wh <- range(wh, ignore.strand = TRUE)
  wh_df <- data.frame(wh)
  stats_plot<- bind_rows(sig_ihs_val, sig_nsl_val) %>% filter(pop == p, genename == g)%>% mutate(pop = p, statname = case_when(statid == 20 ~ "nSL", statid == 26 ~ "iHS"), chrom = paste0('chr',chrom)) %>%  select(SYMBOL = genename, chrom, chrom_start, chrom_end, statname) %>%  bind_rows(., lower_sig_stats %>% filter(pop == p, SYMBOL == g) %>% select(contains("chrom"), statname, SYMBOL) ) %>% GenomicRanges::GRanges() %>% ggbio::autoplot(., aes(fill = statname, colour = statname)) + theme(legend.position = 'bottom')
  #ggbio::fixed(stats_plot) <- TRUE
  
  marker_plot <- markers %>% mutate(chrom = paste0('chr', chrom), chrom_end = chrom_start +1) %>% filter(chrom %in% wh_df$seqnames, chrom_start >= wh_df$start -5000, chrom_start <= wh_df$end + 5000) %>% GRanges() %>%  ggbio::autoplot(. ) + theme(legend.position = 'none')
  
  gene_plot <- ggbio::autoplot(Homo.sapiens::Homo.sapiens, which = wh, gap.geom = "chevron")
  #ggbio::fixed(gene_plot) <- TRUE
  hasAxis(gene_plot) <- TRUE
  ggbio::tracks(marker_plot, gene_plot, stats_plot ) + ggtitle(g)
}

#load the pathways tables from panther
panther <- list()
for (file in list.files(path = 'data/Panther/')){
  panther[[file]] <- read.delim(paste0('data/Panther/',file), header = FALSE, stringsAsFactors = FALSE) %>% mutate(file  = file, pop = substring(file, 1, 3))
}

# load the kegg pathways analysis
kegg <- list()
list.files('~/Git_repos/bookdown_thesis/data/enrichr_kegg2016/', include.dirs = FALSE, pattern = '*.txt')
for(file in list.files('~/Git_repos/bookdown_thesis/data/enrichr_kegg2016/', include.dirs = FALSE, pattern = '*.txt')){
  kegg[[file]] <- read.delim(paste0('~/Git_repos/bookdown_thesis/data/enrichr_kegg2016/', file), header = TRUE, stringsAsFactors = FALSE)
  
}

kegg_sig <-list()
kegg_sig <- lapply(names(kegg), function(n){kegg[[n]] %>% filter(Adjusted.P.value < 0.05) %>% arrange(Adjusted.P.value) %>% mutate(filename = n)} )
names(kegg_sig) <- names(kegg)
kegg_sig<- bind_rows(kegg_sig) %>% mutate(stat = lapply(filename, function(fn){strsplit(fn, '_')[[1]][2]})) %>% mutate(stat = substring(stat, 1, nchar(stat)-4), pop = substring(filename, 1,3)) %>% select(-filename)

# for use with fst data where the pop is "POP1_POP2"
# will return the populations so that the pops are aligned to the specified pop2
popswitch <- function(p, pop2){
  p1 <- substring(p,1,3)
  p2 <- substring(p,5,8)
  if(p2 == pop2){
    return(paste0(p1,p2))
  } else{
    return(paste0(p2,p1))
  }
}


# for a given population pull out the genes that had haplotypic significant results AND intersected a significant window from the intra-population statistics
extract_sig <- function(popname){
  bind_rows(sig_ihs_val, sig_nsl_val) %>% filter(pop == popname) %>% select(pop, 'SYMBOL' = genename) %>% inner_join(., bind_rows(., lower_sig_stats %>% filter(pop == popname) %>% select(pop, SYMBOL)), by = c("SYMBOL","pop")) %>% filter(pop == popname, !is.na(SYMBOL)) %>% .[['SYMBOL']] %>% unique()
}

# for a given super population pull out the genes that intersected a significant window from the intra-population statistics
# only uses genes that have both haplotypic AND SFS
extract_sig_super <- function(superpopname){
  bind_rows(sig_ihs_val, sig_nsl_val) %>% 
    select(pop, superpop, 'SYMBOL' = genename) %>% 
    filter(superpop == superpopname) %>% 
    inner_join(., lower_sig_stats %>%
                 left_join(., panel %>% select(pop, "superpop" = super_pop) %>% distinct(), by = 'pop')%>% 
                 select(superpop, SYMBOL, pop) %>% filter(superpop == superpopname), by = c('pop','superpop', "SYMBOL")) %>% 
    filter(!is.na(SYMBOL)) %>% .[['SYMBOL']] %>% unique()
}


#Calculate the genome wide DHS regions using a threshold of at least 60 cell lines having a score > 400
dhs_clusters <- readr::read_delim('~/data/UCSC_DNase1HypersensitivityClusters/hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed.gz', delim = '\t', col_names = c("chrom","chromStart","chromEnd","name","score","sourceCount","sourceIds","sourceScores")) %>% data.frame()
dhs_regions <- dhs_clusters %>% filter(sourceCount > 60 & score > 400) %>% GRanges()


# for working out where variants were
# snpeff_effect <- function(trans){
#   strsplit(trans, '|', fixed = TRUE)[[1]][[2]]
# }
# 
# 
# extract_snp_effect <-function(info_col){
#   map(info_col , function(x){strsplit(x, split = ';')[[1]]}) %>% map_chr(., function(x){x[grep(x, pattern = 'ANN')]}) %>%  #take info column and get the ANN field
#     map(., function(x){strsplit(x, ',') %>% unlist()}) %>% # split the ANN field into transcripts
#     map_chr(., function(x){lapply(x, snpeff_effect) %>% unlist() %>% unique() %>% paste(collapse = ',')}) # pull out the effect from each transcript and find the unique values
# }
# 
# snpeff <- read.table('~/data/NZ_coreExome_1kgp/snpEff_Annotated/NZ_1KGP.chr_1-22.marker_effects.txt', header = FALSE, stringsAsFactors = FALSE)
# snpeff <- snpeff %>% mutate(effect = pull(., V8) %>% extract_snp_effect()) %>% select(-V8)
# snpeff$l <- strsplit(snpeff$effect, split = ',')
# snpeff <- bind_rows(lapply(seq_along(snpeff[,1]), function(x){eff = unlist(snpeff[x,10]); data.frame(snpeff[rep(x, length(eff)),1:8], eff, stringsAsFactors = FALSE)}))
#saveRDS(snpeff, file = '~/data/NZ_coreExome_1kgp/snpEff_Annotated/snpeff_terms.RDS')
snpeff <- readRDS(file = '~/data/NZ_coreExome_1kgp/snpEff_Annotated/snpeff_terms.RDS')


genelists <- data.frame(pheno = c(rep("obesity", length(obesity_GR$SYMBOL)), rep("urate", length( urate_goutGR$SYMBOL)), rep('t2d', length(t2d_GR$SYMBOL)), rep("kd", length(kd_GR$SYMBOL)), rep("metsyn", length(metsyn_GR$SYMBOL))),  genename = c(obesity_GR$SYMBOL, urate_goutGR$SYMBOL, t2d_GR$SYMBOL, kd_GR$SYMBOL, metsyn_GR$SYMBOL), stringsAsFactors = FALSE) %>% distinct() %>%  mutate(present = 1) %>% spread(pheno, present, 0)


pickrell2009 <- list(
  biaka = 'COX5B,ACTR1B,ZAP70,KIAA1641,MAN2A1,HIST1H4G,HIST1H3F,HIST1H2BH,HIST1H3G,HIST1H2BI,HIST1H4H,BTN3A2,BTN2A2,BTN3A1,BTN2A3,BTN3A3,BTN2A1,BTN1A1,HMGN4,RBMS3,OR10Q1,OR10W1,OR5B17,OR5B3,OR5B2,OR5B12,OR5B21,CPNE4,TBX18,DCAMKL1,SOHLH2,AFF3,SLC28A3,PRTG,NEDD4,SUHW4,TCF12,TBC1D15,TPH2,TCF12,RTN4,FLJ42562,SGEF,KIAA0350,ITGA6,PDK1,CREB5,NSUN6,ARL5B,NEK10,SLC4A7,SCRN1,FKBP14,PLEKHA8,WIPF3,C7orf10,CSNK1G3,SETBP1,PHACTR3,THADA,PLEKHH2,CTNNA3,FLJ37357,SOX6,GDPD4,PAK1,MYO7A,RIOK3,C18orf8,NPC1,ANKRD29,C18orf45,C3orf48,EFHB,RAB5A,PCAF,SGOL1,PFKP,PITRM1,TRIO,EML4,COX7A2L,DNER,PID1,SPOCK1,EPAS1,ATP6V1E2,RHOQ,LOC388965,INPP4A,MGC26733,CNGA3,OPRM1,DKFZP586P0123,PPME1,P4HA3,PGM2L1,KCNE3,B2M,TRIM69,C15orf43,NRXN1,SOD2,WTAP,ACAT2,TCP1,MRPL18,PNLDC1,MAS1,GNA12,CARD11,TTC23,LRRC28,DMN,MAP2,PCDHGA1,PCDHGA2,PCDHGA3,PCDHGB1,PCDHGA4,PCDHGB2,PCDHGA5,PCDHGB3,PCDHGA6,PCDHGA7,PCDHGB4,PCDHGA8,PCDHGB5,PCDHGA9,PCDHGB6,PDCD7,CLPX,CILP,PARP16,PUNC,LRRN1,FER,EPB41L2,TGFBI,SMAD5,TRPC7,FUT8,C14orf145,TSHR,TSCOT,ZFP37,POMC,DNMT3A,SYN3,TIMP3,C8orf34,NPHP4,KCNAB2,DSCAM,EPB41L4B,C9orf4,C9orf5,CTNNAL1,SCLT1,C4orf33,CSS3,COL21A1,GMDS,TAAR2,TAAR5,TAAR1,VNN1,VNN3,VNN2,C6orf192,RPS12' %>% strsplit(., split = ',') %>% unlist() %>% unique(),
  bantu = 'RYR1,SIPA1L3,DPF1,PPP1R14A,SPINT2,LOC541469,C19orf33,YIF1B,KCNK6,C19orf15,PSMD8,GGN,SPRED3,RASGRP4,FAM98C,FDX1,ARHGAP20,PLXDC2,SBF2,ANP32B,C9orf156,FOXE1,HEMGN,HS3ST2,USP31,KIAA1627,SEC24D,SYNPO2,NDST3,PRSS12,MAN2A1,SUMF1,ITPR1,PPP1R12B,SYT2,JARID1B,ETFB,CLDND2,NKG7,LIM2,SIGLEC12,SIGLEC6,ZNF175,SIGLEC5,SIGLEC10,SIGLEC8,LOC729767,KCNJ3,PREX1,ARFGEF2,SYT1,MCM5,RASD2,MB,LOC284912,APOL6,APOL5,SIM1,SUSD1,ROD1,HSDL2,GAN,CMIP,C12orf26,TMTC2,TRPV1,CTNS,TRPV3,CARKL,TAX1BP3,TMEM93,P2RX5,ITGAE,GSG2,ITGB8,LOC440131,PTHLH,FSTL5,NCAM1,RNF32,LMBR1,NOM1,LINGO2,SETBP1,NUBPL,MAMDC4,LCN10,LCN6,LCN8,UNQ2541,TMEM141,C9orf86,PHPT1,EDF1,TRAF2,FBXW5,C8G,PTGDS,FLJ45224,C9orf142,CLIC3,ABCA2,C9orf139,FUT7,KIAA1984,LCN1,OR52L1,OR56A1,OR56B4,OR52B2,OR52W1,C11orf42,C11orf56,CNGA4,CCKBR,OR56A4,CSMD3,ANKRD19,BICD2,ZNF484,FGD3,SAMD3,KIAA1913,COL8A1,FLJ44076,C3orf26,FILIP1L,RNF165,LOXHD1,RIMS1,EPB41L4B,C9orf4,C9orf5,CTNNAL1,PCDH17,FLJ40296,TSCOT,ZFP37,CCDC66,C3orf63,ARHGEF3,CCDC102B,DLG2,APOL3,APOL4,APOL2,APOL1,MYH9,RAB38,FYCO1,CXCR6,XCR1,CCR1,PPFIA2,HS3ST4,FHIT,HIST1H4G,HIST1H3F,HIST1H2BH,HIST1H3G,HIST1H2BI,HIST1H4H,BTN3A2,BTN2A2,BTN3A1,BTN2A3,BTN3A3,BTN2A1,BTN1A1,HMGN4,AXIIR,ZNF131,MGC42105,CCT8,C21orf7,BACH1,MAS1L,OR2H2,GABBR1,MOG,HLA−F,INPP4A,MGC26733,CNGA3,MGC40405,CDK6,DKFZP564O0523,PEX1,RFXDC2,TEX9,SEMA3E,AGK,FLJ40852,SSBP1,TAS2R3,TAS2R4,TAS2R5,LOC136242,KIAA1147,RGS6,TNN,KIAA0040,TNR,SCLT1,C4orf33,PKN2,PDE9A,WDR4,NDUFV3,WDR40C,C10orf59,LIPL1,LIPF,CPA6,DEPDC2,LMBRD1,COL19A1,TCERG1,GPR151,PPP2R2B,CDC27,MYL4,ITGB3,C17orf57,VPS13C,C21orf29,KRTAP10−1,KRTAP10−2,KRTAP10−3,KRTAP10−4,KRTAP10−5,KRTAP10−6,KRTAP10−7,KRTAP10−8,KRTAP10−9,KRTAP10−10,KRTAP10−11,KRTAP12−4,KR,JMJD1A,VPS24,RNF103,RMND5A,CTNNA3,ISL2,ZNF291,ETFA' %>% strsplit(., split = ',') %>% unlist() %>% unique(),
  europe = 'CMYA3,RGS9,SLC44A1,TMEM30B,PRKCH,HIF1A,DCP2,MCC,RFX3,FLJ35880,NGLY1,OXSM,PDE11A,AGPS,TTC30A,TTC30B,PCMTD1,C8ORFK32,COL22A1,KCND2,TSPAN12,HLA−A,HCG9,ZNRD1,RNF39,TRIM31,TRIM40,TRIM15,TRIM26,FLJ45422,TRIM39,RPP21,PPP1R11,LAMA3,C18orf17,GTF3C3,UNQ6411,PGAP1,ANKRD44,TCERG1,GPR151,PPP2R2B,COL22A1,CDH12,ACVR1C,ACVR1,UPP2,GLCE,PAQR5,PREP,TYRP1,C9orf150,NRG3,HMGCR,COL4A3BP,POLK,MPHOSPH6,SNX25,LRP2BP,ANKRD37,C4orf20,CCDC110,RHOJ,GPHB5,PPP2R5E,DENND1A,FGF1,ARHGAP26,PDZRN4,ELMO1,EDNRA,AK3L1,AK3L2,JAK1,PTPRM,ATF6,OLFML2B,NOS1AP,C1orf111,SH2D1B,GRIK1,CCDC53,NUP37,PMCH,C12orf48,COL6A3,PRLH,RAB17,MLPH,NRG3,FAM109A,SH2B3,ATXN2,CUTL2,SMYD3,NRSN1,DCDC2,TMEM66,LEPROTL1,DCTN6,NTRK2,FERD3L,TMEM100,PCTP,MAF,AUH,NFIL3,SLC8A1,AGBL1,FBXO2,FBXO44,FBXO6,EPB41L4B,C9orf4,C9orf5,CTNNAL1,APBA2,NELL1,CDH10,SYNE2,KLHL15,EIF2S3,ZFX,IER3,PRR3,HLA−E,GNL1,ABCF1,PPP1R10,MRPS18B,C6orf134,C6orf136,DHX16,NRM,MDC1,TUBB,FLOT1,DKFZP686A01247,PHOX2B,TMEM33,RYBP,COL8A1,TUBA3C,USP15,MON2,WDR16,STX8,USP43,LCN2,CIZ1,DNM1,GOLGA2,TRUB2,COQ4,SLC27A4,URM1,CEECAM1,C9orf16,PRB4,PRB1,PRB2,NCOR2,FAM101A,C1GALT1,COL28A1,WWOX,MAPK8IP3,IFT140,CRAMP1L,HN1L,NME3,MRPS34,EME2,SPSB3,NUBP2,IGFALS,FAHD1,C16orf73,HAGH,DRG1,EIF4ENIF1,SFI1,PISD,C22orf30,GRIK2,RPS6KA2,RNASET2,FGFR1OP,COMMD10,SEMA6A,ZPLD1,TIMP2,LGALS3BP,CANT1,C1QTNF1,FLJ21865,DNAH7,STK17B,HECW2' %>% strsplit(., split = ',') %>% unlist() %>% unique(),
  mideast = 'SLC44A1,ERBB4,KCNH5,RHOJ,GPHB5,PPP2R5E,CDH18,C8ORFK32,COL22A1,RGS9,APBA2,KCND2,TSPAN12,NGLY1,OXSM,CCDC102B,NRG3,C20orf19,XRN2,NKX2−2,PDZRN4,CRB2,DENND1A,PRKCH,HIF1A,BACH1,GRIK1,FLJ20699,PPARA,LOC150383,PKDREJ,GTSE1,TRMU,CELSR1,SMYD3,RAVER2,JAK1,COL6A3,PRLH,RAB17,MLPH,TRPS1,COL22A1,TMEM132B,TCERG1,GPR151,PPP2R2B,PRDM10,C11orf37,APLP2,ST14,ZBTB44,SCML4,FLJ10159,ACVR1C,ACVR1,UPP2,FGF1,ARHGAP26,COMMD3,FLJ35880,TSLP,TMPRSS7,PLCXD2,PHLDB2,ABHD10,TAGLN3,CYP26B1,BBS9,GLCE,PAQR5,KIF23,RPLP1,LOC729416,MAPK8IP3,IFT140,CRAMP1L,HN1L,NME3,MRPS34,EME2,SPSB3,NUBP2,IGFALS,FAHD1,C16orf73,HAGH,WWOX,ANKS1B,FLJ46363,CLEC12A,CLEC1B,CLEC12B,CLEC9A,CLEC1A,CLEC7A,GABARAPL1,KLRD1,KLRK1,C12orf59,OLR1,GPR64,MLC1,BRD1,ZBED4,ALG12,CRELD2,PIM3,FLJ41993,LOC164714,FOXA2,ERN1,TEX2,PECAM1,PRICKLE2,SYNE2,KITLG,PARD6B,BCAS4,ADNP,DPM1,MOCS3,PCMTD1,CBR4,SH3RF1,NEK1,DLG2,C14orf124,TGM1,CHMP4A,MDP−1,NEDD8,GMPR2,TINF2,RABGGTA,DHRS1,C14orf21,LTB4R2,LTB4R,ADCY4,NFATC4,KIAA0323,CMA1,RIPK3,CIDEB,CBLN3,MYT1L,TMEM100,PCTP,SMYD2,PTPN14,NRG3,OR2T8,OR2L13,OR2L8,OR2AK2,OR2L2,OR2L3,OR2M5,OR2M2,OR2M3,CDH12,EDNRA,POT1,GPR37,LOC401398,LYPLA1,MRPL15,DOK5,MYOD1,SERGEF,KCNC1,NEB,ARL5A,CACNB4,CDH10,HLA−A,HCG9,ZNRD1,RNF39,TRIM31,TRIM40,TRIM15,TRIM26,FLJ45422,TRIM39,RPP21,PPP1R11,FAM109A,SH2B3,ATXN2,CUTL2,ST6GALNAC1,MXRA7,SFRS2,MFSD11,MGAT5B,PTDSR,LOC124512,STK32B,COMMD10,SEMA6A,CBLN4,CAND1,FAM83B,HTR1F,CGGBP1,GRM8,PHF17,SCLT1,C4orf33,DCP2,MCC,SMAD6,PTPRM,HMGCR,COL4A3BP,POLK,CHRNA6,THAP1,RNF170,CHRNB3,SLC24A5,MYEF2' %>% strsplit(., split = ',') %>% unlist() %>% unique(),
  sasia = 'SLC24A5,MYEF2,SLC12A1,DUT,TMEM100,PCTP,KCND2,TSPAN12,KCNH5,RHOJ,GPHB5,PPP2R5E,C8orf15,MTMR9,AMAC1L2,XKR6,FLJ20699,PPARA,LOC150383,PKDREJ,GTSE1,TRMU,CELSR1,NRG3,NGLY1,OXSM,C6orf148,KCNQ5,C8ORFK32,COL22A1,CRB2,DENND1A,PRR3,RNF39,TRIM31,TRIM40,TRIM15,TRIM26,FLJ45422,TRIM39,RPP21,HLA−E,GNL1,ABCF1,FBXW2,PSMD5,PHF19,TRAF1,C5,CEP110,RAB14,FLJ35880,CXCR7,LGALS8,HEATR1,ACTN2,MTR,SLC44A1,N6AMT1,ZNF294,C20orf133,ACVR1C,ACVR1,UPP2,GPR55,CAB39,ITM2C,PSMD1,LCN2,CIZ1,DNM1,GOLGA2,TRUB2,COQ4,SLC27A4,URM1,CEECAM1,C9orf16,LIN28B,BVES,POPDC3,PREP,APBA2,TCERG1,GPR151,PPP2R2B,RGS9,SLC25A36,SPSB4,ACPL2,DKFZP686A01247,PHOX2B,TMEM33,CCDC102B,DCP2,MCC,GPBAR1,PNKD,ARPC2,TMBIM1,MGC50811,SLC11A1,CTDSP1,VIL1,USP37,AAMP,TMTC2,NUDT18,DOK2,XPO7,NPM2,FGF17,EPB49,RAI16,HR,USP25,C21orf34,KIAA1324L,DMTF1,RARB,TOP2B,ACSL3,KCNE4,ANKS1B,PLEKHG3,SPTB,CHURC1,GPX2,RAB15,FNTB,MAX,PDE11A,AGPS,TTC30A,TTC30B,NTRK2,PCMTD1,EFCBP1,TMEM55A,OTUD6B,PRDM10,C11orf37,APLP2,ST14,ZBTB44,ROR2,SPTLC1,PARD6B,BCAS4,ADNP,DPM1,MOCS3,RP1L1,C8orf74,SOX7,PINX1,XKR6,CHRNA6,THAP1,RNF170,CHRNB3,IPMK,ZCD1,UBE2D1,TFAM,CCT8,C21orf7,BACH1,GRIK1,ELMO1,HTR3A,ZBTB16,ISX,C12orf64,GALNT17,CBR4,SH3RF1,SYNE2,SMAD6,NDUFB6,TAF1L,LOC401498,DNAH7,STK17B,HECW2,PLCB4,C20orf103,PAK7,PLCZ1,CAPZA3,COL22A1,FOXA2,TSPYL5,UBE2N,MRPL42,SOCS2,CRADD,OR2B3,OR5U1,OR5V1,OR12D3,LOC651503,ZNF783,PDIA4,ZNF786,ZNF425,ZNF398,ZNF212,ZNF282,KITLG,CDH13,HSBP1,EFCBP2,MBTPS1,OSGIN1,MLYCD,LOC146167,ADRBK2,NAALADL2,FTMT,DIAPH3' %>% strsplit(., split = ',') %>% unlist() %>% unique(),
  easia = 'RANBP2,FLJ32745,EDAR,ERBB4,HEY2,NCOA7,SLC5A7,SULT1C3,SULT1C1,SULT1C2,GCC2,PCDH15,EXOC6,CYP26C1,CYP26A1,FER1L3,CEP55,APC,SRP19,DCP2,MCC,REEP5,ERC1,MRPL27,EME1,LRRC59,FLJ20920,CHAD,RSAD1,EPN3,SPATA20,CACNA1G,ABCC3,ANKRD40,CROP,XYLT2,MYCBPAP,CBR4,SH3RF1,NEK1,NRG1,MAP3K7IP1,MGAT3,SMCR7L,ATF4,RPS19BP1,CACNA1I,FER,PRKG1,KIF11,HHEX,EXOC6,LOC284757,CRB2,DENND1A,XYLT1,GRM3,KIAA1324L,ITGA6,PDK1,CYP1B1,MGC34824,ARL6IP2,LEPR,PDE4B,CXXC4,DAPK2,FAM96A,SNX1,SNX22,PPIB,CSNK1G1,SDK1,CUGBP2,ADAMTS18,SLC44A5,ACADM,RABGGTB,MSH4,EPB41L4B,C9orf4,C9orf5,CTNNAL1,COL11A1,DYNLT1,VIL2,LOC202459,SYTL3,FOXP1,FLJ34931,CLIP4,ALK,APOL3,APOL4,APOL2,APOL1,MYH9,TXN2,FOXRED2,EIF3S7,USP38,GAB1,ARHGEF7,C13orf16,CADPS,C10orf38,ITGA8,GRM8,ZNF800,GCC1,ARF5,FSCN3,PAX4,NAALADL2,ZNF277P,IFRD1,FLJ39575,NEK10,SLC4A7,UBE2D2,CXXC5,PSD2,NRG2,ZPBP,SERPINI1,GOLPH4,C9orf93,IRF2BP2,TNNI3K,C1orf173,PDE11A,C6orf10,TNXB,CREBL1,FKBPL,PRRT1,PPT2,EGFL8,AGPAT1,RNF5,AGER,GPSM3,NOTCH4,LOC401252,PBX2,TNFSF18,TNFSF4,FAM44B,LIN28B,BVES,POPDC3,PREP,ZC3HAV1L,ZC3HAV1,TM4SF19,FLJ25996,ZDHHC19,OSTalpha,PCYT1A,MGC33212,C8orf15,MTMR9,AMAC1L2,XKR6,ASCC1,C10orf104,DDIT4,DNAJB12,CBARA1,BOC,WDR52,CCDC52,THADA,PLEKHH2' %>% strsplit(., split = ',') %>% unlist() %>% unique(),
  america = 'PRM1,KIAA0350,SOCS1,TNP2,PRM3,PRM2,C16orf75,LOC400499,GPR115,GPR111,CD2AP,C1QDC1,TSPAN11,HSF2,SERINC1,PKIB,HDAC3,C5orf16,FCHSD1,CENTD3,PCDH1,DIAPH1,TMPRSS11F,CHR415SYT,TMPRSS11B,FUT10,RBM13,C8orf41,RNF122,DUSP26,C16orf68,BCAS1,CYP24A1,FLJ16641,NRCAM,PNPLA8,THAP5,DNAJB9,FCGR2A,HSPA6,FCGR3A,FCGR2C,FCGR3B,FCGR2B,FCRLA,FCRLB,DUSP12,ATF6,OLFML2B,C21orf87,BRWD1,HMGN1,WRB,C21orf13,SH3BGR,LOC728776,PKIB,FABP7,SMPDL3A,MUC17,TRIM56,SERPINE1,AP1S1,VGF,FLJ39237,MOGAT3,PLOD3,ZNHIT1,RIMS1,RAD54B,KIAA1429,RBM35A,DPY19L4,ZNF659,FAM83B,PRMT8,EFCAB4B,PARP11,NEK7,SSBP2,B3GALT5,PCP4,LOC150084,GALNT10,SAP30L,HAND1,CCND1,FLJ42258,ORAOV1,FGF19,TSSC1,TTC15,LHX1,AATF,ACACA,C22orf34,BRD1,ZBED4,SPTBN1,RTN4,FLJ42562,KIAA0409,ILK,TAF10,DCHS1,MRPL17,OR2AG2,OR2AG1,OR6A2,OR10A5,OR10A2,TPP1,KCNC2,PLXNA4B,SMC6,FLJ40869,KCNS3,MCMDC1,ASF1A,C6orf60,SSH1,DAO,SVOP,USP30,ALKBH2,UNG,WDR25,KIAA1446,C14orf70,FBXL17,LRRC55,AGTRL1,TNKS1BP1,BHLHB3,SSPN,ITPR2,ITGB6,RBMS1,FLJ21908,P11,RAPGEF3,RNF144,RSAD2,OPRM1,PIP3−E,FBXO47,PLXDC1,CACNB1,RPL19,STAC2,FHOD3,C18orf10,KIAA1328,ZNF438,KRT28,KRT40,KRT25,KRT26,KRT27,KRT10,KRT12,KRT20,KRT23,KRT39,KRTAP3−3,KRTAP3−2,KRTAP3−1,KRTAP1−5,KRTAP1−3,TMEM99,GPR113,SELI,C2orf39,OTOF,HADHB,THRB,VCX−C,KAL1,NOSTRIN,SPBC25,G6PC2,ABCB11,DHRS9,GPR116,GPR110,AQP9,CENPA,KCNK3,DPYSL5,MAPRE3,C2orf18,MAN1A1,GBP6,LRRC8B,APOBEC1,GDF3,DPPA3,CLEC4C,NANOG,SLC13A1,PTPRA,FAM113A,VPS16,GNRH2,MRPS26,OXT,AVP,UBOX5,CD302,LY75,PLA2R1,KIAA1370,ONECUT1,C6orf107,TAF11,ANKS1A,TCP11,SLC16A10,KIAA1919,REV3L,LARP2,PGRMC2,CYP39A1,TDRD6,PLA2G7,SLC25A27,ADAM10,SLTM,FAM63B' %>% strsplit(., split = ',') %>% unlist() %>% unique(),
  oceania = 'FAM77D,SUCLG2,GRB14,COBLL1,FLJ39822,SMOC1,SLC8A3,OR5U1,OR5V1,OR12D3,OR12D2,OR11A1,OR10C1,OR2H1,MAS1L,C3orf57,LOC131149,CRADD,SCN2A,SCN3A,FAM130A2,MTBP,SNTB1,TSSC1,TTC15,ADI1,RNASEH1,RPS7,COLEC11,RGS7BP,P18SRP,C20orf196,CHGB,CGI−09,MCM8,CRLS1,C20orf75,C20orf42,EDEM3,FAM129A,CPEB2,FLJ39743,IGF1R,C13orf21,CTGF,ENPP1,WDR72,C9orf126,PPP6C,RABEPK,HSPA5,DUPD1,DUSP13,SAMD8,VDAC2,MYST4,ANKRD10,ARHGEF7,WWOX,LRRC40,SFRS11,LRRC7,PPM1L,B3GALNT1,NMD3,SGCZ,KIF6,LRCH1,ESD,FAM77D,C4orf16,TIFA,ALPK1,ING3,FLJ21986,WNT16,FAM3C,IGF2R,PNLDC1,MAS1,PDE11A,AGPS,TTC30A,TTC30B,CCDC46,ACTN1,WDR22,ADAM11,DBF4B,HIGD1B,EFTUD2,CCDC103,GFAP,C1QL1,GJA7,LOC146909,ERBB4,FGL1,PCM1,ASAH1,CRB2,STRBP,DENND1A,LACE1,FOXO3A,LRRC18,C10orf72,XYLB,ENDOGL1,ACVR2B,SCN5A,ERBB4,SDCCAG10,ADAMTS6,LRRC3B,TNRC9,FGF14,C1orf21,ECM1,ADAMTSL4,C1orf138,MCL1,ENSA,GOLPH3L,HORMAD1,CTSS,CTSK,ARNT,RBM43,NMI,RASSF4,C10orf10,C10orf25,ZNF22,FRMPD4,ANKRD15,DMRT1,DMRT3,DMRT2,NFAM1,SERHL,CGI−96,FEM1C,TICAM2,TMED7,PTPRT,PHC3,PRKCI,SKIL,CLDN11,RFTN2,MARS2,BOLL,LHFPL3,ORC5L,FERD3L,C11orf69,CD59,FBXO3,LMO2,LCP2,KCNIP1,KCNMB1,ASTN1,FAM5B,ARFGAP3,PACSIN2,TTLL1,BIK,IKZF2,ZNF239,ZNF485,ZNF32,COL21A1,DST,ADAMTS20,COQ5,COX6A1,TRIAP1,15E1.2,SFRS9,DYNLL1,RNF10,POP5,CABP1,KIAA0152,ACADS,MGC5139,ANKRD12,TWSG1,PCSK2,BFSP1,ADK,PDE4D,STK31' %>% strsplit(., split = ',') %>% unlist() %>% unique()
)
