# load gene x count matrix
library(dplyr)
library(data.table)
data_path <- '/home/rstudio/result/res20230521/Swanda2023/cpm/count_mat_1.csv.gz'
#data_path <- '/home/rstudio/result/res20230521/Swanda2023/cpm/count_mat_2.csv.gz'

dt = fread(data_path,header=TRUE) %>% as.matrix(rownames=1)
smpls <- colnames(dt)
grp <- unlist(lapply(smpls, function(x){strsplit(x , "_[12]")[[1]][1]}))
genes <- rownames(dt)

library(edgeR)
d <- DGEList(
  counts=dt,
  genes=genes,
  group=grp)

# filering
d_full <- d
idx_keep <- rowSums(cpm(d)>100) >= 1
d <- d[idx_keep,]
dim(d)
d$samples$lib.size <- colSums(d$counts)

# normalize data by TMM
d <- calcNormFactors(d,method = 'TMM')
nc <- cpm(d, normalized.lib.sizes=TRUE) %>% as.data.table
genes_now <- d$genes %>% as.vector
row.names(nc) <- genes_now$genes
fwrite(
  nc,
  '/home/rstudio/result/res20230521/Swanda2023/cpm/count_mat_TMM_1.csv',
  row.names = TRUE)

# fix the dispersion value
d <- estimateDisp(d)

pairs <- list(
  c("mRNA_NT","mRNA_-AA"),
  c("mRNA_NT","mRNA_-Cys")
)
#pairs <- list(
#  c("NT","CTNS_KD"),
#  c("NT","CTNS_KD_RESCUE"),
#  c("NT","Cysteamine")
#)

for (p in pairs){
  et <- exactTest(d,pair=p)
  fwrite(
    et$table,
    paste0(
      '/home/rstudio/result/res20230521/Swanda2023/cpm/exact_test_',
      p[1],'_',p[2],
      '.csv'),
    row.names = TRUE)
}