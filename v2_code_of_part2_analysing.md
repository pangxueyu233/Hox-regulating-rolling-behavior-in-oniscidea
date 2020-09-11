

## Figure2 making

```R
Summary <- read.csv("/mnt/data/user_data/xiangyu/workshop/roly_poly/alignment.csv")
Summary$Map_ratio <- 100*(Summary$Mapped_reads/Summary$Total_reads)
Summary$Forward_strand_ratio <- 100*(Summary$Forward_strand/Summary$Total_reads)
Summary$Both_pairs_mapped_ratio <- 100*(Summary$Both_pairs_mapped/Summary$Total_reads)
Summary$Reverse_strand_ratio <- 100*(Summary$Reverse_strand/Summary$Total_reads)
library(ggpubr)
p1 <- ggbarplot(Summary, x = "data_type", y = "Map_ratio", legend = "none", fill="data_type",
       add = c("mean_se", "jitter"),title="Map_ratio (%)")
p2 <- ggbarplot(Summary, x = "data_type", y = "Forward_strand_ratio", legend = "none", fill="data_type",
       add = c("mean_se", "jitter"),title="Forward_strand_ratio (%)")
p3 <- ggbarplot(Summary, x = "data_type", y = "Reverse_strand_ratio", legend = "none", fill="data_type",
       add = c("mean_se", "jitter"),title="Reverse_strand_ratio (%)")
p4 <- ggbarplot(Summary, x = "data_type", y = "Both_pairs_mapped_ratio", legend = "none", fill="data_type",
       add = c("mean_se", "jitter"),title="Both_pairs_mapped_ratio (%)")
aa <- plot_grid(p1,p2,p3,p4,nrow=1)
ggsave("/mnt/data/user_data/xiangyu/workshop/roly_poly/figure_making/map_ratio_sum.svg", plot=aa,width = 14, height = 5,dpi=1080)
```

![image-20200810105009657](v2_code_of_part2_analysing.assets/image-20200810105009657.png)

```R
library(ggplot2)
library(ggmsa)
f <- "/mnt/data/user_data/xiangyu/workshop/roly_poly/similarity/final_out/hox/sub_fa/all_merge_hoxc10.fa"
ggmsa(f, font = NULL, color = "Chemistry_NT", seq_name = TRUE)
library(phangorn)
library(ggtree)
tipseq <- read.phyDat(f, format = "fasta")
bears <- read.phyDat(f,format="fasta",type="DNA")
dm <- dist.ml(bears)
treeUPGMA <- upgma(dm)
rootedtree <- root(treeUPGMA, outgroup='hoxc10_ref')
library(ggtree)
p <- ggtree(rootedtree, branch.length='none') + geom_tiplab()
aa <- msaplot(p, f, offset=1.5) + labs(title="Hoxc10")
ggsave("/mnt/data/user_data/xiangyu/workshop/roly_poly/figure_making/hoxc10_tree.svg", plot=aa,width = 10, height = 5,dpi=1080)

```



![image-20200810105810672](v2_code_of_part2_analysing.assets/image-20200810105810672.png)

```R
library(ggplot2)
library(ggmsa)
f <- "/mnt/data/user_data/xiangyu/workshop/roly_poly/similarity/final_out/hox/sub_fa/all_merge_hoxd11.fa"
ggmsa(f, font = NULL, color = "Chemistry_NT", seq_name = TRUE)
library(phangorn)
library(ggtree)
tipseq <- read.phyDat(f, format = "fasta")
bears <- read.phyDat(f,format="fasta",type="DNA")
dm <- dist.ml(bears)
treeUPGMA <- upgma(dm)
rootedtree <- root(treeUPGMA, outgroup='hoxd11_ref')
library(ggtree)
p <- ggtree(rootedtree, branch.length='none') + geom_tiplab()
aa <- msaplot(p, f, offset=1.5) + labs(title="Hoxd11")
ggsave("/mnt/data/user_data/xiangyu/workshop/roly_poly/figure_making/Hoxd11_tree.svg", plot=aa,width = 10, height = 5,dpi=1080)
```

![image-20200810105851557](v2_code_of_part2_analysing.assets/image-20200810105851557.png)

~~~R
library(ggplot2)
library(ggmsa)
library(Biostrings)
f <- "/mnt/data/user_data/xiangyu/workshop/roly_poly/AA_SIM/sub_fa/all_merge_Hoxc9_AA.fa"
library(Biostrings)
x <- readAAStringSet(f)
d <- as.dist(stringDist(x, method = "hamming")/width(x)[1])
library(ape)
tree <- bionj(d)
library(ggtree)
p <- ggtree(tree, branch.length='none') + geom_tiplab()
data <- tidy_msa(x)
p <- p + geom_facet(geom = geom_msa, data = data,  panel = 'msa',
               font = NULL, color = "Chemistry_AA") +
    xlim_tree(1) + labs(title="Hoxc9")
~~~

<img src="code_of_part2_analysing.assets/image-20200830145041883.png" alt="image-20200830145041883" style="zoom:50%;" />

~~~R
library(ggplot2)
library(ggmsa)
library(Biostrings)
f <- "/mnt/data/user_data/xiangyu/workshop/roly_poly/AA_SIM/sub_fa/all_merge_Hoxc10_AA.fa"
library(Biostrings)
x <- readAAStringSet(f)
d <- as.dist(stringDist(x, method = "hamming")/width(x)[1])
library(ape)
tree <- bionj(d)
library(ggtree)
p <- ggtree(tree, branch.length='none') + geom_tiplab()
data <- tidy_msa(x)
p <- p + geom_facet(geom = geom_msa, data = data,  panel = 'msa',
               font = NULL, color = "Chemistry_AA") +
    xlim_tree(1) + labs(title="Hoxc10")
~~~

<img src="code_of_part2_analysing.assets/image-20200830145115984.png" alt="image-20200830145115984" style="zoom:50%;" />



## Figure6 making

```R
counts_res_1 <- read.csv(file = "/mnt/data/user_data/xiangyu/workshop/roly_poly/bwa_files/mm10/All_data_normalised_and_res.csv")
counts_res_1 <- na.omit(counts_res_1)
counts_res_1 <- subset(counts_res_1,padj < 0.05)
counts_res_1 <- subset(counts_res_1,abs(log2FoldChange) > 0.5)
counts_res_1 <- counts_res_1[order(counts_res_1$log2FoldChange,decreasing=TRUE),]
counts <- counts_res_1[,c(1:7)]
rownames(counts) <- counts$X
counts <- counts[,-1]
sampletable <- data.frame(bam_files=c("curl1.sort.bam","curl2.sort.bam",
  "curl3.sort.bam","strait1.sort.bam",
  "strait2.sort.bam","strait3.sort.bam"),
  sample_name=c("curl1","curl2","curl3","strait1","strait2","strait3"),
  group=c("curl","curl","curl","strait","strait","strait")
  )
rownames(sampletable) <- sampletable$sample_name
Roly_Poly_seurat <- CreateSeuratObject(counts, project='Roly_Poly',meta.data=sampletable)
Idents(Roly_Poly_seurat) <- Roly_Poly_seurat$group
mark_gene <- grep("Hox",rownames(counts),value=TRUE)
library(RColorBrewer)
pdf(file="/mnt/data/user_data/xiangyu/workshop/roly_poly/figure_making/mm10_DGE.pdf",width = 5, height = 5)
aa <- XY_heatmap(seurat_obj=Roly_Poly_seurat,group="group",genes=rownames(counts),all_num=FALSE,new_names=NULL,labels_rot=90,
  assay_sel="RNA",color=colorRampPalette(c("navy", "white", "firebrick3"))(50)[50:1],min_and_max_cut=2,show_row_names=FALSE,mark_gene=mark_gene,label_size=0)
dev.off()
```

![image-20200810110049340](v2_code_of_part2_analysing.assets/image-20200810110049340.png)

```R
GOdownres_1_all <- readRDS("/mnt/data/user_data/xiangyu/workshop/roly_poly/bwa_files/mm10/Curl_DN_GO_BP.rds")
GOupres_1_all <- readRDS("/mnt/data/user_data/xiangyu/workshop/roly_poly/bwa_files/mm10/Curl_UP_GO_BP.rds")
ff <- barplot(GOupres_1_all,showCategory=5)
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
p1 <- ff1 + labs(title="Curl vs Strait UP")

ff <- barplot(GOdownres_1_all,showCategory=5)
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
p2 <- ff1+ labs(title="Curl vs Strait DN")
aa <- cowplot::plot_grid(p1,p2)
aa
```

![image-20200812201411988](v2_code_of_part2_analysing.assets/image-20200812201411988.png)

```R
library(ggplot2)
library(ggmsa)
f <- "/mnt/data/user_data/xiangyu/workshop/roly_poly/similarity/final_out/hox/sub_fa/all_merge_hoxc9.fa"
ggmsa(f, font = NULL, color = "Chemistry_NT", seq_name = TRUE)
library(phangorn)
library(ggtree)
tipseq <- read.phyDat(f, format = "fasta")
bears <- read.phyDat(f,format="fasta",type="DNA")
dm <- dist.ml(bears)
treeUPGMA <- upgma(dm)
rootedtree <- root(treeUPGMA, outgroup='hoxc9_ref')
library(ggtree)
p <- ggtree(rootedtree, branch.length='none') + geom_tiplab()
aa <- msaplot(p, f, offset=1.5) + labs(title="Hoxc9")
ggsave("/mnt/data/user_data/xiangyu/workshop/roly_poly/figure_making/hoxc9_tree.svg", plot=aa,width = 10, height = 5,dpi=1080)

```





## sessionInfo

```R
sessionInfo()
```

```R
R version 3.6.3 (2020-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.6 LTS

Matrix products: default
BLAS:   /usr/lib/libblas/libblas.so.3.6.0
LAPACK: /usr/lib/lapack/liblapack.so.3.6.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=zh_CN.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=zh_CN.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=zh_CN.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=zh_CN.UTF-8 LC_IDENTIFICATION=C

attached base packages:
 [1] splines   stats4    parallel  stats     graphics  grDevices utils
 [8] datasets  methods   base

other attached packages:
 [1] ggtree_2.0.2                ape_5.3
 [3] ggmsa_0.0.4                 gdtools_0.2.2
 [5] ggpubr_0.2.5                magrittr_1.5
 [7] BuenColors_0.5.5            MASS_7.3-51.5
 [9] scales_1.1.0                future.apply_1.5.0
[11] future_1.17.0               tidyr_1.0.2
[13] nichenetr_0.1.0             iTALK_0.1.0
[15] stringr_1.4.0               data.table_1.12.8
[17] DESeq2_1.26.0               org.Mm.eg.db_3.10.0
[19] PoiClaClu_1.0.2.1           RColorBrewer_1.1-2
[21] pheatmap_1.0.12             GenomicAlignments_1.22.1
[23] SummarizedExperiment_1.16.1 DelayedArray_0.12.2
[25] BiocParallel_1.20.1         matrixStats_0.56.0
[27] GenomicFeatures_1.38.2      Rsamtools_2.2.3
[29] Biostrings_2.54.0           XVector_0.26.0
[31] GenomicRanges_1.38.0        GenomeInfoDb_1.22.1
[33] trqwe_0.1                   velocyto.R_0.6
[35] cowplot_1.0.0               pathview_1.26.0
[37] org.Hs.eg.db_3.10.0         topGO_2.38.1
[39] SparseM_1.78                GO.db_3.10.0
[41] AnnotationDbi_1.48.0        IRanges_2.20.2
[43] S4Vectors_0.24.3            graph_1.64.0
[45] clusterProfiler_3.14.3      DOSE_3.12.0
[47] plyr_1.8.6                  monocle_2.14.0
[49] DDRTree_0.1.5               VGAM_1.1-2
[51] ggplot2_3.2.1               Biobase_2.46.0
[53] BiocGenerics_0.32.0         irlba_2.3.3
[55] densityClust_0.3            Rtsne_0.15
[57] gplots_3.0.3                proxy_0.4-23
[59] Matrix_1.2-18               Seurat_3.1.5
[61] dplyr_0.8.5

loaded via a namespace (and not attached):
  [1] minpack.lm_1.2-1
  [2] graphlayouts_0.6.0
  [3] pbapply_1.4-2
  [4] lattice_0.20-40
  [5] haven_2.2.0
  [6] vctrs_0.2.4
  [7] V8_3.0.2
  [8] fastICA_1.2-2
  [9] mgcv_1.8-31
 [10] blob_1.2.1
 [11] survival_3.1-11
 [12] prodlim_2019.11.13
 [13] DBI_1.1.0
 [14] SingleCellExperiment_1.8.0
 [15] rappdirs_0.3.1
 [16] uwot_0.1.8
 [17] jpeg_0.1-8.1
 [18] zlibbioc_1.32.0
 [19] MatrixModels_0.4-1
 [20] pcaMethods_1.78.0
 [21] ChIPseeker_1.25.0
 [22] htmlwidgets_1.5.1
 [23] mvtnorm_1.1-0
 [24] GlobalOptions_0.1.1
 [25] miscTools_0.6-26
 [26] laeken_0.5.1
 [27] leiden_0.3.3
 [28] DEoptimR_1.0-8
 [29] tidygraph_1.1.2
 [30] Rcpp_1.0.4.6
 [31] readr_1.3.1
 [32] KernSmooth_2.23-16
 [33] gdata_2.18.0
 [34] limma_3.42.2
 [35] Hmisc_4.4-0
 [36] ShortRead_1.44.3
 [37] RSpectra_0.16-0
 [38] fastmatch_1.1-0
 [39] ranger_0.12.1
 [40] digest_0.6.25
 [41] png_0.1-7
 [42] qlcMatrix_0.9.7
 [43] sctransform_0.2.1
 [44] ggraph_2.0.2
 [45] pkgconfig_2.0.3
 [46] docopt_0.6.1
 [47] gower_0.2.1
 [48] iterators_1.0.12
 [49] URD_1.1.0
 [50] reticulate_1.15
 [51] network_1.16.0
 [52] circlize_0.4.8
 [53] modeltools_0.2-23
 [54] xfun_0.13
 [55] zoo_1.8-7
 [56] tidyselect_1.0.0
 [57] reshape2_1.4.4
 [58] purrr_0.3.4
 [59] ica_1.0-2
 [60] viridisLite_0.3.0
 [61] rtracklayer_1.46.0
 [62] rlang_0.4.5
 [63] hexbin_1.28.1
 [64] glue_1.4.0
 [65] ensembldb_2.10.2
 [66] ggseqlogo_0.1
 [67] lava_1.6.7
 [68] europepmc_0.3
 [69] ggsignif_0.6.0
 [70] recipes_0.1.10
 [71] labeling_0.3
 [72] chipseq_1.36.0
 [73] DEsingle_1.6.0
 [74] gggenes_0.4.0
 [75] class_7.3-15
 [76] preprocessCore_1.48.0
 [77] RMTstat_0.3
 [78] Sierra_0.2.0
 [79] DO.db_2.9
 [80] annotate_1.64.0
 [81] jsonlite_1.6.1
 [82] systemfonts_0.2.0
 [83] bit_1.1-15.2
 [84] gridExtra_2.3
 [85] gmodels_2.18.1
 [86] stringi_1.4.6
 [87] bitops_1.0-6
 [88] RSQLite_2.2.0
 [89] randomForest_4.6-14
 [90] KEGGgraph_1.46.0
 [91] rstudioapi_0.11
 [92] nlme_3.1-147
 [93] qvalue_2.18.0
 [94] locfit_1.5-9.4
 [95] VariantAnnotation_1.32.0
 [96] listenv_0.8.0
 [97] ggthemes_4.2.0
 [98] gridGraphics_0.5-0
 [99] dbplyr_1.4.3
[100] TTR_0.23-6
[101] readxl_1.3.1
[102] lifecycle_0.2.0
[103] ggfittext_0.8.1
[104] timeDate_3043.102
[105] munsell_0.5.0
[106] cellranger_1.1.0
[107] hwriter_1.3.2
[108] visNetwork_2.0.9
[109] caTools_1.18.0
[110] codetools_0.2-16
[111] lmtest_0.9-37
[112] htmlTable_1.13.3
[113] triebeard_0.3.0
[114] lsei_1.2-0
[115] xtable_1.8-4
[116] ROCR_1.0-7
[117] diptest_0.75-7
[118] BiocManager_1.30.10
[119] Signac_0.2.5
[120] scatterplot3d_0.3-41
[121] abind_1.4-5
[122] farver_2.0.3
[123] FNN_1.1.3
[124] RANN_2.6.1
[125] askpass_1.1
[126] biovizBase_1.34.1
[127] sparsesvd_0.2
[128] RcppAnnoy_0.0.16
[129] patchwork_1.0.0.9000
[130] tibble_3.0.1
[131] dichromat_2.0-0
[132] cluster_2.1.0
[133] tidytree_0.3.3
[134] ellipsis_0.3.0
[135] prettyunits_1.1.1
[136] lubridate_1.7.8
[137] ggridges_0.5.2
[138] igraph_1.2.5
[139] RcppEigen_0.3.3.7.0
[140] fgsea_1.12.0
[141] slam_0.1-47
[142] destiny_3.0.1
[143] VIM_5.1.1
[144] htmltools_0.4.0
[145] BiocFileCache_1.10.2
[146] plotly_4.9.2.1
[147] XML_3.99-0.3
[148] ModelMetrics_1.2.2.2
[149] e1071_1.7-3
[150] foreign_0.8-76
[151] withr_2.2.0
[152] fitdistrplus_1.0-14
[153] randomcoloR_1.1.0.1
[154] bit64_0.9-7
[155] foreach_1.5.0
[156] ProtGenerics_1.18.0
[157] robustbase_0.93-6
[158] scde_2.14.0
[159] combinat_0.0-8
[160] GOSemSim_2.13.1
[161] rsvd_1.0.3
[162] memoise_1.1.0
[163] forcats_0.5.0
[164] rio_0.5.16
[165] geneplotter_1.64.0
[166] gamlss_5.1-6
[167] gamlss.data_5.1-4
[168] caret_6.0-86
[169] curl_4.3
[170] DiagrammeR_1.0.5
[171] fdrtool_1.2.15
[172] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
[173] urltools_1.7.3
[174] xts_0.12-0
[175] acepack_1.4.1
[176] edgeR_3.28.1
[177] checkmate_2.0.0
[178] Gviz_1.30.3
[179] npsurv_0.4-0
[180] maxLik_1.3-8
[181] rjson_0.2.20
[182] openxlsx_4.1.4
[183] ggrepel_0.8.2
[184] distillery_1.0-7
[185] Lmoments_1.3-1
[186] tools_3.6.3
[187] sandwich_2.5-1
[188] soGGi_1.18.0
[189] RCurl_1.98-1.2
[190] car_3.0-7
[191] ggplotify_0.0.5
[192] xml2_1.3.2
[193] httr_1.4.1
[194] assertthat_0.2.1
[195] AnnotationFilter_1.10.0
[196] boot_1.3-24
[197] globals_0.12.5
[198] R6_2.4.1
[199] nnet_7.3-13
[200] RcppHNSW_0.2.0
[201] treeio_1.10.0
[202] progress_1.2.2
[203] genefilter_1.68.0
[204] KEGGREST_1.26.1
[205] gtools_3.8.2
[206] shape_1.4.4
[207] Rook_1.1-1
[208] carData_3.0-3
[209] colorspace_1.4-1
[210] generics_0.0.2
[211] base64enc_0.1-3
[212] smoother_1.1
[213] pillar_1.4.3
[214] Rgraphviz_2.30.0
[215] tweenr_1.0.1
[216] sp_1.4-1
[217] ggplot.multistats_1.0.0
[218] HSMMSingleCell_1.6.0
[219] rvcheck_0.1.8
[220] GenomeInfoDbData_1.2.2
[221] extRemes_2.0-11
[222] gtable_0.3.0
[223] bdsmatrix_1.3-4
[224] zip_2.0.4
[225] knitr_1.28
[226] RcppArmadillo_0.9.860.2.0
[227] latticeExtra_0.6-29
[228] gamlss.dist_5.1-6
[229] biomaRt_2.42.1
[230] Cairo_1.5-12
[231] pscl_1.5.5
[232] flexmix_2.3-15
[233] quantreg_5.55
[234] vcd_1.4-7
[235] BSgenome_1.54.0
[236] openssl_1.4.1
[237] backports_1.1.6
[238] plotrix_3.7-8
[239] ipred_0.9-9
[240] enrichplot_1.6.1
[241] brew_1.0-6
[242] hms_0.5.3
[243] ggforce_0.3.1
[244] polyclip_1.10-0
[245] grid_3.6.3
[246] numDeriv_2016.8-1.1
[247] bbmle_1.0.23.1
[248] lazyeval_0.2.2
[249] Formula_1.2-3
[250] tsne_0.1-3
[251] crayon_1.3.4
[252] MAST_1.12.0
[253] pROC_1.16.2
[254] svglite_1.2.3
[255] viridis_0.5.1
[256] rpart_4.1-15
[257] compiler_3.6.3
```

