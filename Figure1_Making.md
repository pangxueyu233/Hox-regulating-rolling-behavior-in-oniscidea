# Figure1 mkaing and Code of analysing

## 1. Loading required packages

```R
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(Matrix)
  library(proxy)
  library(gplots)
  library(Rtsne)
  library(densityClust)
  library(irlba)
  library(monocle)
  library(plyr)
  library(DOSE)
  library(clusterProfiler)
  library(topGO)
  library(AnnotationDbi)
  library(cowplot)
  library(ggplot2)
  library(trqwe)
  library(Rsamtools)
  library(GenomicFeatures)
  library(GenomicAlignments)
  library(BiocParallel)
  library(pheatmap)
  library(RColorBrewer)
  library(PoiClaClu)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(DESeq2)
  library(stringr)
  library(tidyr)
})
source("/mnt/data/user_data/xiangyu/programme/R_PACKAGES/my_code/MyBestFunction_scRNA.R")
library(future)
library(future.apply)
options(future.globals.maxSize = 300 * 1024^3)
plan("multiprocess", workers = 15)
plan()
library(scales)
library(BuenColors)			
```

## Figure1 making

~~~R
AA <- data.frame(type=c("Roll1", "Roll2", "Roll3", "Roll4", "Roll5", "Straight1", "Straight2", "Straight3", "Straight4", "Straight5"),
  length=c(0.9, 1.1, 0.7, 0.6, 0.5, 0.8, 0.8, 0.8, 0.6, 0.6),
  new_name=c("Roll", "Roll", "Roll", "Roll", "Roll", "Straight", "Straight", "Straight", "Straight", "Straight"))
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
 data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}

df3 <- data_summary(AA, varname="length", 
                    groupnames=c("new_name"))
df3$new_name=as.factor(df3$new_name)
p <- ggplot(df3,aes(x=new_name, y=length, fill=new_name))
p <- p + geom_bar(stat="identity", position=position_dodge())
p <- p + geom_errorbar(aes(ymin=length-sd, ymax=length+sd), width=.2,
                 position=position_dodge(.9))
p
~~~

![image-20200813180410925](Figure1_Making.assets/image-20200813180410925.png)